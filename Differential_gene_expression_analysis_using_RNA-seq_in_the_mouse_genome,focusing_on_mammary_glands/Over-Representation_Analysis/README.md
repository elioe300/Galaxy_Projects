# Over-Representation Analysis in Differentially Expressed Genes (DEG)

## Introduction

In a differential expression analysis, we identify genes that show significant changes in their expression between experimental conditions. However, the list of DEGs alone does not provide information about the affected biological processes. To interpret these results in a functional context, an Over-Representation Analysis (ORA) is performed.

---

## Why Perform Over-Representation Analysis?

ORA helps answer key questions, such as:
- Which metabolic pathways or biological processes are enriched in the DEGs?
- What are the most represented functional categories in upregulated or downregulated genes?

This analysis compares the set of DEGs against reference databases such as Gene Ontology, KEGG, or Reactome to identify whether any gene set is significantly represented in the studied list.

---

## Steps for Over-Representation Analysis

1. **Data Preprocessing:**
   - Obtain DEGs.
```r
# Get the genes present in the dataset
genes_in_dataset <- rownames(res_corrected)
```

2. **Gene Filtering:**
   - Remove non-significant genes based on adjusted p-values, in this case **FDR > 0.05**.
   - Classify genes into upregulated and downregulated lists.
```r
# Filter significant genes
significant_genes <- data[data$DE != 'NO', ]

# Split genes into lists of upregulated and downregulated
deg_list <- split(significant_genes, significant_genes$DE)

# Define filtering thresholds
padj_threshold <- 0.05   # Significance threshold for adjusted p-values
gene_count_threshold <- 5 # Minimum number of genes in a pathway
```

3. **Definition of Gene Background:**
   - Use the total set of genes measured in the experiment as a reference to avoid biases.
```r  
# Read .gmt file with pathway annotations
gmt_file <- "./m5.go.bp.v2024.1.Mm.entrez.gmt"
pathways_data <- read.gmt(gmt_file)

# Filter pathways containing genes present in the dataset
filtered_pathways <- pathways_data[pathways_data$gene %in% genes_in_dataset, ]
```

4. **Enrichment Analysis:**
   - Use tools such as `clusterProfiler` in R to perform the over-representation test.
   - Compare DEGs with pathway and biological process databases, in this case with GO.
```r  
enrichment_results <- lapply(names(deg_list), function(category) {
  enricher(gene = rownames(deg_list[[category]]), TERM2GENE = filtered_pathways)
})
names(enrichment_results) <- names(deg_list)
```

5. **Filtering Results:**
   - Select only significantly enriched terms based on an adjusted p-value threshold (< 0.05) and a minimum number of genes per pathway.
```r  
selected_pathways <- unique(results_df$ID[results_df$p.adjust < padj_threshold & results_df$Count > gene_count_threshold])
final_results <- results_df[results_df$ID %in% selected_pathways, ]
```

6. **Result Visualization:**
   - Generate tables to visualize the enriched pathways.
```r 
filename = "GO.RDS"
write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)
```

---

## Implementation in R

This analysis was carried out in R using the `clusterProfiler` package, based on the tutorial by BiostatSquid [(link to tutorial)](https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/).

The complete code used for this analysis can be found in the following script: 
- [clusterProfiler.R](./Recursos/clusterProfiler.R)

This file includes all the steps described above to perform the over-representation analysis in a set of DEGs.

I would like to express my gratitude to Laura, the creator of BiostatSquid, for her dedication in creating these tutorials. Her efforts in teaching and simplifying complex concepts are a great help for newly graduated bioinformaticians.

---
