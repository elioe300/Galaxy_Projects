# Instalación y carga de paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)

# Para visualización
to_install <- c("pheatmap", "DOSE", "enrichplot", "ggupset")
install.packages(to_install)

# Funciones ------------------------------------------------
## Función: Convertir una matriz de adyacencia a lista ####
matrix_to_list <- function(pws){
  pathways_list <- list()
  for (pathway in colnames(pws)) {
    pathways_list[[pathway]] <- rownames(pws)[as.logical(pws[, pathway])]
  }
  return(pathways_list)
}

# Obtener los genes presentes en el conjunto de datos
genes_in_dataset <- rownames(res_corrected)

# Leer archivo .gmt con anotaciones de rutas metabólicas
gmt_file <- "./m5.go.bp.v2024.1.Mm.entrez.gmt"
pathways_data <- read.gmt(gmt_file)

# Filtrar rutas metabólicas que contienen genes presentes en el dataset
filtered_pathways <- pathways_data[pathways_data$gene %in% genes_in_dataset, ] 


# Filtrar genes significativos
significant_genes <- data[data$DE != 'NO', ]

# Separar genes en listas de regulados al alza y a la baja
deg_list <- split(significant_genes, significant_genes$DE)

# Definir umbrales para filtrado
padj_threshold <- 0.05   # Umbral de significancia para p-ajustada
gene_count_threshold <- 5 # Número mínimo de genes en una ruta

# Aplicar clusterProfiler en cada subconjunto de DEG
enrichment_results <- lapply(names(deg_list), function(category) {
  enricher(gene = rownames(deg_list[[category]]), TERM2GENE = filtered_pathways)
})
names(enrichment_results) <- names(deg_list)

# Convertir los resultados a un dataframe
results_df <- lapply(names(enrichment_results), function(category) {
  rbind(enrichment_results[[category]]@result)
})
names(results_df) <- names(enrichment_results)
results_df <- do.call(rbind, results_df)

# Agregar columnas adicionales con transformación de datos
results_df <- results_df %>% mutate(
  minuslog10padj = -log10(p.adjust),
  pathway_name = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(results_df))
)

# Filtrar rutas metabólicas significativas
selected_pathways <- unique(results_df$ID[results_df$p.adjust < padj_threshold & results_df$Count > gene_count_threshold])
final_results <- results_df[results_df$ID %in% selected_pathways, ]

# Guardar el resultado
print('Saving clusterprofiler results')
filename = "GO.RDS"
write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)