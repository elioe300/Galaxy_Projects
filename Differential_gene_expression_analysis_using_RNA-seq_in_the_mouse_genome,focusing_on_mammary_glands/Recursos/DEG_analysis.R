# Instalación y carga de paquetes necesarios
paquetes <- c("BiocManager", "edgeR", "org.Mm.eg.db", "ggplot2", "biomaRt", "statmod", "limma")

# Instalar paquetes que falten
a.instalar <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(a.instalar)) install.packages(a.instalar)

library(edgeR)
library(org.Mm.eg.db)
library(ggplot2)
library(biomaRt)
library(statmod)
library(limma)

# Asegurar que la conexión a internet no cause errores por tiempo de espera
options(timeout = 300)

# Definir metadatos de las muestras
metadata <- data.frame(
  SampleName = c("MCL1.DI", "MCL1.DJ", "MCL1.DK", "MCL1.DL"),
  group = c("basal.pregnant", "basal.pregnant", "basal.lactate", "basal.lactate")
)

# Cargar datos de conteo de genes desde archivos de featureCounts
MCL1.DI <- read.delim("./Galaxy271-[featureCounts on data 215 and data 269_ Counts].tabular", col.names = "MCL1.DI")
MCL1.DJ <- read.delim("./Galaxy273-[featureCounts on data 215 and data 270_ Counts].tabular", col.names = "MCL1.DJ")
MCL1.DK <- read.delim("./Galaxy236-[featureCounts on data 215 and data 43_ Counts].tabular", col.names = "MCL1.DK")
MCL1.DL <- read.delim("./Galaxy234-[featureCounts on data 215 and data 42_ Counts].tabular", col.names = "MCL1.DL")

# Combinar los datos en una única tabla de conteos
counts <- cbind.data.frame(MCL1.DI, MCL1.DJ, MCL1.DK, MCL1.DL, row.names = row.names.data.frame(MCL1.DI))

# Verificar que todas las tablas de conteos tienen las mismas filas (mismos genes)
all_identical <- all(sapply(list(MCL1.DI, MCL1.DJ, MCL1.DK, MCL1.DL), rownames) == rownames(MCL1.DI))

# Conectar a la base de datos Ensembl para obtener anotaciones de genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")

# Obtener las anotaciones de los genes basadas en Ensembl Gene ID
annotations <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = rownames(MCL1.DI),
  mart = ensembl
)

# Filtrar anotaciones para conservar solo genes con IDs válidos en Entrez
annotations <- annotations[!is.na(annotations$entrezgene_id), ]

# Eliminar duplicados en las anotaciones
annotations <- annotations[!duplicated(annotations$entrezgene_id), ]

# Filtrar la tabla de conteos para conservar solo genes presentes en las anotaciones
counts <- counts[rownames(counts) %in% annotations$ensembl_gene_id, ]

# Asegurar que las anotaciones estén en el mismo orden que los rownames de counts
annotations <- annotations[match(rownames(counts), annotations$ensembl_gene_id), ]

# Reemplazar rownames de counts con IDs de Entrez
rownames(counts) <- annotations$entrezgene_id

# Limpiar la descripción de los genes eliminando información innecesaria
annotations$description <- gsub("\\[Source:.*", "", annotations$description)

# Crear objeto DGEList para análisis de expresión diferencial
counts.DGEList <- DGEList(counts)

# Asignar los grupos experimentales a las muestras
counts.DGEList$samples$group <- factor(metadata$group)

# Asociar las anotaciones con el objeto DGEList
counts.DGEList$annotations <- annotations

# Filtrar genes con baja expresión
keep <- filterByExpr(counts.DGEList)
counts.DGEList <- counts.DGEList[keep, keep.lib.sizes = FALSE]

# Normalización de los datos
counts.DGEList <- calcNormFactors(counts.DGEList)

# Ver información de las muestras
counts.DGEList$samples

# Crear matriz de diseño para el modelo lineal de análisis de expresión diferencial
design <- model.matrix(~ 0 + counts.DGEList$samples$group)
colnames(design) <- levels(counts.DGEList$samples$group)

# Estimación de la dispersión
dispersion <- estimateDisp(counts.DGEList, design, robust = TRUE)

# Ajustar modelo lineal generalizado para análisis de expresión diferencial
fit <- glmQLFit(counts.DGEList, design, robust = TRUE)

# Definir contraste para comparar grupos (basal pregnant vs basal lactate)
B.LvsP <- makeContrasts(basal.pregnant - basal.lactate, levels = design)
res <- glmQLFTest(fit, contrast = B.LvsP)

# Obtener resultados ajustados por FDR
res_corrected <- topTags(res, n = Inf)

# Determinar genes diferencialmente expresados (DEGs) usando umbral de p-valor ajustado y logFC
is.de1 <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(is.de1)

# Generar tabla de resultados con estado de regulación diferencial
data <- res_corrected$table
data$DE <- "NO"
data$DE[data$logFC > 2 & data$FDR < 0.05] <- "up"
data$DE[data$logFC < -2 & data$FDR < 0.05] <- "down"

# Crear gráfico de dispersión (Volcano plot) para visualizar genes diferencialmente expresados
ggplot(data, aes(x = logFC, y = -log10(FDR), col = DE)) +
  geom_point(size = 0.2) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10(FDR)")

# Análisis de ontología génica (Gene Ontology, GO) para genes diferencialmente expresados
go <- goana(res, species = "Mm")
# Mostrar las 5 categorías más representativas
topGO(go, number = 5)
