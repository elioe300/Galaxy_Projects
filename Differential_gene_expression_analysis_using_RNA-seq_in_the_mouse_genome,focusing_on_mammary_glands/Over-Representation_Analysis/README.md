# Análisis de Sobre-representación en Genes Diferencialmente Expresados (DEG)

## Introducción

En un análisis de expresión diferencial, identificamos genes que presentan cambios significativos en su expresión entre condiciones experimentales. Sin embargo, el listado de DEG por sí solo no proporciona información sobre los procesos biológicos afectados. Para interpretar estos resultados en un contexto funcional, se realiza un análisis de sobre-representación (ORA, por sus siglas en inglés: Over-Representation Analysis).

## ¿Por qué realizar un análisis de sobre-representación?

El ORA permite responder preguntas clave, como:
- ¿Cuáles son las rutas metabólicas o procesos biológicos enriquecidos en los DEG?
- ¿Cuáles son las categorías funcionales más representadas en los genes regulados al alza o a la baja?

Este análisis se basa en comparar el conjunto de DEG con bases de datos de referencia como Gene Ontology, KEGG o Reactome para identificar si algún conjunto de genes está significativamente representado en la lista de genes estudiados.

## Pasos del análisis de sobre-representación

1. **Preprocesamiento de datos:**
   - Obtención de DEG.
```r
# Obtener los genes presentes en el conjunto de datos
genes_in_dataset <- rownames(res_corrected)

```

2. **Filtrado de genes:**
   - Eliminación de genes no significativos según valores de p ajustados, en este caso **FDR > 0.05**.
   - Clasificación en genes regulados al alza y a la baja.
```r
# Filtrar genes significativos
significant_genes <- data[data$DE != 'NO', ]

# Separar genes en listas de regulados al alza y a la baja
deg_list <- split(significant_genes, significant_genes$DE)

# Definir umbrales para filtrado
padj_threshold <- 0.05   # Umbral de significancia para p-ajustada
gene_count_threshold <- 5 # Número mínimo de genes en una ruta
```

3. **Definición del fondo de genes:**
   - Se usa el conjunto total de genes medidos en el experimento como referencia para evitar sesgos.
```r  
# Leer archivo .gmt con anotaciones de rutas metabólicas
gmt_file <- "./m5.go.bp.v2024.1.Mm.entrez.gmt"
pathways_data <- read.gmt(gmt_file)

# Filtrar rutas metabólicas que contienen genes presentes en el dataset
filtered_pathways <- pathways_data[pathways_data$gene %in% genes_in_dataset, ]
```

4. **Análisis de enriquecimiento:**
   - Se emplean herramientas como `clusterProfiler` en R para realizar la prueba de sobre-representación.
   - Se comparan los DEG con bases de datos de vías y procesos biológicos, en este caso con GO.
```r  
enrichment_results <- lapply(names(deg_list), function(category) {
  enricher(gene = rownames(deg_list[[category]]), TERM2GENE = filtered_pathways)
})
names(enrichment_results) <- names(deg_list)
```

5. **Filtrado de resultados:**
   - Se seleccionan solo los términos significativamente enriquecidos según un umbral de p-ajustada (< 0.05) y un mínimo de genes por vía.
```r  
selected_pathways <- unique(results_df$ID[results_df$p.adjust < padj_threshold & results_df$Count > gene_count_threshold])
final_results <- results_df[results_df$ID %in% selected_pathways, ]
```
6. **Visualización de resultados:**
   - Se generan tablas para visualizar las vías enriquecidas.
```r 
filename = "GO.RDS"
write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)
```

## Implementación en R

Para realizar este análisis en R, se utilizó el paquete `clusterProfiler`, basado en el tutorial de BiostatSquid [(enlace al tutorial)](https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/).

El código completo utilizado para este análisis se encuentra en el siguiente script: 
- [clusterProfiler.R](./Recursos/clusterProfiler.R)

En este archivo se incluyen todos los pasos descritos anteriormente para realizar el análisis de sobre-representación en un conjunto de DEG.

Quiero agradecer a Laura, la persona detrás de biostatsquid, por su dedicación al crear estos tutoriales. Su esfuerzo por enseñar y divulgar conceptos complejos de forma sencilla es de gran ayuda para bioinformáticos recién graduados.

(Dentro de nada, actualizaré este readme para visualizar los [resultados obtenidos](./Recursos/GO.RDS_resclusterp.csv) en el script.) 
---


