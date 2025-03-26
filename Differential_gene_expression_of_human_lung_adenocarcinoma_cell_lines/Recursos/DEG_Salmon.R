# Verificar e instalar paquetes necesarios
dependencies <- c("BiocManager", "tximport", "tidyr", "dplyr", "org.Hs.eg.db", 
                  "edgeR", "ggplot2", "biomaRt", "statmod", "limma", "rtracklayer", "AnnotationDbi")
missing_packages <- dependencies[!(dependencies %in% installed.packages()[,"Package"])]
if(length(missing_packages)) install.packages(missing_packages)

# Cargar librerías
lapply(dependencies, library, character.only = TRUE)

# Definir nombres de muestras y grupos experimentales
samples <- c("H1975_1", "H1975_2", "HCC827_1", "HCC827_2")
group <- factor(gsub("_.*", "", samples))

# Buscar archivos de cuantificación de Salmon
salmon_files <- list.files(path = "./", pattern = "_quant.sf", full.names = TRUE)

# Importar archivo de anotación GTF y crear mapa de transcritos a genes
gtf <- rtracklayer::import("gencode.v47.annotation.gtf")
gtf_df <- as.data.frame(gtf)
tx2gene_map <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id)

# Cargar datos de cuantificación de Salmon
tx_data <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene_map)

# Extraer matriz de conteos y asignar nombres de columnas
counts <- tx_data[["counts"]]
colnames(counts) <- samples
rownames(counts) <- sub("\\.\\d+$", "", rownames(counts))

# Anotación de genes utilizando org.Hs.eg.db
gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(counts),
                                          columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                                          keytype = "ENSEMBL")

# Limpiar y ordenar anotaciones
gene_annotations <- gene_annotations[!is.na(gene_annotations$SYMBOL), ]
gene_annotations <- gene_annotations[!duplicated(gene_annotations$SYMBOL), ]
counts <- counts[rownames(counts) %in% gene_annotations$ENSEMBL, ]
gene_annotations <- gene_annotations[match(rownames(counts), gene_annotations$ENSEMBL), ]
rownames(counts) <- gene_annotations$ENTREZID

# Calcular Counts Per Million (CPM)
cpm <- edgeR::cpm(counts)
log2_cpm <- edgeR::cpm(counts, log=TRUE)

# Convertir matriz a tibble para manipulación con tidyverse
log2_cpm_df <- as_tibble(log2_cpm, rownames = "geneID")
log2_cpm_pivot <- pivot_longer(log2_cpm_df, cols = -1, names_to = "samples", values_to = "expression")

# Visualización de la distribución de expresión con Violin Plot
ggplot(log2_cpm_pivot, aes(x = samples, y = expression, fill = samples)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y = "log2 expression", x = "Sample", title = "Log2 Counts per Million (CPM)", subtitle = "Unfiltered, non-normalized") +
  theme_bw()

# Crear objeto DGEList para análisis de expresión diferencial
dge_list <- DGEList(counts)
dge_list$samples$group <- group
dge_list$annotations <- gene_annotations

# Filtrar genes con baja expresión
dge_list <- dge_list[filterByExpr(dge_list), keep.lib.sizes = FALSE]
log2.cpm.filtered <- edgeR::cpm(counts.DGEList, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")

# Gráfico de expresión después de filtrado de baja expresión
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized") +
  theme_bw()

# Normalización usando método TMM
dge_list <- calcNormFactors(dge_list, method = "TMM")
log2_cpm_filtered <- cpm(dge_list, log=TRUE)
log2_cpm_filtered_df <- as_tibble(log2_cpm_filtered, rownames = "geneID")
log2_cpm_filtered_pivot <- pivot_longer(log2_cpm_filtered_df, cols = -1, names_to = "samples", values_to = "expression")

# Gráfico de expresión después de método TM
ggplot(log2_cpm_filtered_pivot, aes(x = samples, y = expression, fill = samples)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y = "log2 expression", x = "Sample", title = "Log2 Counts per Million (CPM)", subtitle = "Filtered and normalized") +
  theme_bw()

# Multidimensional Scaling Plot
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(dge_list, col = colors[group], pch = pch[group])

# Análisis de expresión diferencial
design <- model.matrix(~ 0 + group, data = dge_list$samples)
colnames(design) <- levels(group)

dispersion <- estimateDisp(dge_list, design, robust = TRUE)
fit <- glmQLFit(dge_list, design, robust = TRUE)
contrast <- makeContrasts(HCC827 - H1975, levels = design)
result <- glmQLFTest(fit, contrast = contrast)
is.de1 <- decideTests(result, adjust.method = "BH", p.value = 0.05, lfc=1)
summary(is.de1)

# Obtener resultados ajustados
adjusted_results <- topTags(result, n = Inf)
data_results <- adjusted_results$table

# Histogramas de valores p
hist(data_results$PValue, breaks = 50, col = 'grey', main = 'Distribución de valores p', xlab = 'P-value')
