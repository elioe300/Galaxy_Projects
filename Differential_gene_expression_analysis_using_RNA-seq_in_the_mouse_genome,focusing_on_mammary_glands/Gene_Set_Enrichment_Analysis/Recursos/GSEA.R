# Carga de librerías necesarias
library("fgsea")

# Funciones ===================================================
## Función: Convertir matriz de adyacencia a lista -------------------------
matrix_to_list <- function(pws){
  pathways_list <- list()
  for (pathway in colnames(pws)) {
    pathways_list[[pathway]] <- rownames(pws)[as.logical(pws[, pathway])]
  }
  return(pathways_list)
}

## Función: Preparar archivo GMT --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # Leer archivo GMT
  gmt <- gmtPathways(gmt_file)
  all_genes <- unique(unlist(gmt))
  
  # Convertir GMT a matriz binaria con genes como filas y términos como columnas
  mat <- matrix(NA, dimnames = list(all_genes, names(gmt)),
                nrow = length(all_genes), ncol = length(gmt))
  for (i in 1:ncol(mat)){
    mat[,i] <- as.numeric(all_genes %in% gmt[[i]])
  }
  
  # Filtrar solo genes presentes en los datos de expresión
  genes_filtered <- intersect(genes_in_data, all_genes)
  mat <- mat[genes_filtered, colnames(mat)[which(colSums(mat[genes_filtered,]) > 5)]]
  
  # Convertir matriz de nuevo a lista de conjuntos de genes
  final_list <- matrix_to_list(mat)
  
  if (savefile) {
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print("Conversión exitosa de .gmt a lista de conjuntos de genes!")
  return(final_list)
}

# Análisis ====================================================
## 1. Preparación del fondo de genes -----------------------------------------------
my_genes <- rownames(res_corrected)
bg_genes <- prepare_gmt("./m5.go.bp.v2024.1.Mm.entrez.gmt", my_genes, savefile = FALSE)

# 2. Construcción del ranking de genes ----------------------------------------------
rankings <- sign(res_corrected$table$logFC) * (-log10(res_corrected$table$PValue))
names(rankings) <- rownames(res_corrected)

rankings <- sort(rankings, decreasing = TRUE)
plot(rankings) # Visualizar distribución del ranking

# Manejo de valores extremos en el ranking
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE)

# 3. Ejecución del análisis GSEA ---------------------------------------------------------------
GSEA_results <- fgsea(
  pathways = bg_genes,  # Lista de conjuntos de genes a evaluar
  stats = rankings,
  scoreType = 'std',  # Tipo de ranking (std: valores positivos y negativos)
  minSize = 10,
  maxSize = 500,
  nproc = 1  # Número de núcleos para paralelización
)

# 4. Filtrado y visualización de resultados ------------------------------------------------------
# Obtener las 6 vías más enriquecidas
head(GSEA_results[order(pval), ])

# Definir cuántas vías enriquecidas (positivas y negativas) se incluirán en el análisis.
number_of_top_pathways_up = 10
number_of_top_pathways_down = 10

# Seleccionar las vías enriquecidas más relevantes, tanto positivas como negativas.
topPathwaysUp <- GSEA_results[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEA_results[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]

# Combinar las vías enriquecidas positivas y negativas.
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# Guardar un gráfico con las vías enriquecidas seleccionadas
pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
dev.off()

