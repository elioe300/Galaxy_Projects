# Análisis de Enriquecimiento de Conjuntos de Genes (GSEA)

## Introducción

El análisis de enriquecimiento de conjuntos de genes (Gene Set Enrichment Analysis, GSEA) es una técnica utilizada para determinar si un conjunto predefinido de genes muestra diferencias significativas en su expresión entre diferentes condiciones experimentales. A diferencia del análisis de sobre-representación (ORA), GSEA no requiere un umbral arbitrario para definir genes diferencialmente expresados, sino que utiliza toda la lista de genes ordenados según su nivel de expresión.

## ¿Por qué realizar un análisis GSEA?

El GSEA es útil porque:
- Considera la información de todos los genes en el experimento, evitando sesgos derivados de umbrales arbitrarios de significancia.
- Detecta cambios en vías biológicas incluso cuando los genes individuales no son altamente significativos.
- Es particularmente útil cuando los efectos biológicos son sutiles y afectan múltiples genes de manera coordinada.

Este análisis compara la distribución de un conjunto de genes dentro de una lista ordenada para determinar si están agrupados en la parte superior o inferior de la clasificación, lo que indicaría su enriquecimiento en una condición biológica específica.

## Pasos del Análisis GSEA

1. **Preprocesamiento de datos:**
   - Obtención de una lista de genes con su nivel de expresión y estadística de significancia.
   - Conversión de identificadores de genes a un formato compatible con bases de datos como KEGG, Reactome o Gene Ontology (GO).
```r
   library("fgsea")
   my_genes <- rownames(res_corrected)
   gmt_file <- "./m5.go.bp.v2024.1.Mm.entrez.gmt"
   bg_genes <- prepare_gmt(gmt_file, my_genes, savefile = FALSE)
```
2. **Construcción del ranking de genes:**
   - Se asigna un puntaje a cada gen en función de su expresión y significancia estadística.
   - Se ordenan los genes según este puntaje para evaluar la tendencia de enriquecimiento.
```r
   rankings <- sign(res_corrected$table$logFC) * (-log10(res_corrected$table$PValue))
   names(rankings) <- rownames(res_corrected)
   rankings <- sort(rankings, decreasing = TRUE)
   plot(rankings) # Visualizar distribución del ranking
```

3. **Definición del fondo de genes:**
   - Se utiliza un conjunto de referencia que incluye todos los genes detectados en el experimento.
```r
   max_ranking <- max(rankings[is.finite(rankings)])
   min_ranking <- min(rankings[is.finite(rankings)])
   rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
   rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
   rankings <- sort(rankings, decreasing = TRUE)
```

4. **Ejecución de GSEA:**
   - Se compara la distribución de los genes en una vía específica con el ranking global.
   - Se calcula un puntaje de enriquecimiento normalizado (NES) para cada conjunto de genes.
```r
   GSEA_results <- fgsea(
     pathways = bg_genes,  # Lista de conjuntos de genes a evaluar
     stats = rankings,
     scoreType = 'std',  # Tipo de ranking (std: valores positivos y negativos)
     minSize = 10,
     maxSize = 500,
     nproc = 1  # Número de núcleos para paralelización
   )
```

5. **Filtrado de resultados:**
   - Se seleccionan las vías significativamente enriquecidas según un umbral de p-ajustada (< 0.05).
   - Se eliminan conjuntos de genes redundantes para obtener resultados más claros.
```r
   collapsed_pathways <- collapsePathways(GSEA_results[order(pval)][pval < 0.05], bg_genes, rankings)
   main_pathways <- GSEA_results[pathway %in% collapsed_pathways$mainPathways][order(-NES), pathway]
```

6. **Visualización de resultados:**
   - Se generan gráficos de enriquecimiento y tablas de vías biológicas afectadas.
```r
   plotGseaTable(bg_genes[main_pathways], rankings, GSEA_results, gseaParam = 0.5)
   plotEnrichment(bg_genes[[head(GSEA_results[order(padj), ], 1)$pathway]], rankings) + 
     labs(title = head(GSEA_results[order(padj), ], 1)$pathway)
```
## Implementación en R

Para realizar este análisis en R, se utilizó el paquete `fgsea`. Para realizar este análisis en R, se utilizó el paquete clusterProfiler, basado en el tutorial de BiostatSquid [(enlace al tutorial)](https://biostatsquid.com/fgsea-tutorial-gsea/).
El código completo con todos los pasos necesarios se encuentra en el siguiente script:
- [GSEA.R](./Recursos/GSEA.R)

En este archivo se incluyen todos los pasos descritos anteriormente para realizar el análisis de enriquecimiento de conjuntos de genes utilizando datos de expresión diferencial y bases de datos funcionales.

Quiero agradecer a Laura, la persona detrás de biostatsquid, por su dedicación al crear estos tutoriales. Su esfuerzo por enseñar y divulgar conceptos complejos de forma sencilla es de gran ayuda para bioinformáticos recién graduados.

---


