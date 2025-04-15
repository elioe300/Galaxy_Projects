# Gene Set Enrichment Analysis (GSEA)

## Introduction

Gene Set Enrichment Analysis (GSEA) is a technique used to determine whether a predefined set of genes shows significant differences in expression between different experimental conditions. Unlike over-representation analysis (ORA), GSEA does not require an arbitrary threshold to define differentially expressed genes; instead, it utilizes the entire list of genes ranked by their expression levels.

---

## Why Perform GSEA?

GSEA is valuable because:
- It considers information from all genes in the experiment, avoiding biases from arbitrary significance thresholds.
- It detects changes in biological pathways even when individual genes are not highly significant.
- It is particularly useful when biological effects are subtle and affect multiple genes in a coordinated manner.

This analysis compares the distribution of a gene set within a ranked list to determine whether they are grouped at the top or bottom of the ranking, indicating enrichment under a specific biological condition.

---

## Steps of GSEA

1. **Data Preprocessing:**
   - Obtain a gene list with expression levels and significance statistics.
   - Convert gene identifiers to a format compatible with databases such as KEGG, Reactome, or Gene Ontology (GO).
```r
   library("fgsea")
   my_genes <- rownames(res_corrected)
   gmt_file <- "./m5.go.bp.v2024.1.Mm.entrez.gmt"
   bg_genes <- prepare_gmt(gmt_file, my_genes, savefile = FALSE)
```

2. **Gene Ranking Construction:**
   - Assign a score to each gene based on its expression and statistical significance.
   - Sort the genes by this score to evaluate the enrichment trend.
```r
   rankings <- sign(res_corrected$table$logFC) * (-log10(res_corrected$table$PValue))
   names(rankings) <- rownames(res_corrected)
   rankings <- sort(rankings, decreasing = TRUE)
   plot(rankings) # Visualize the ranking distribution
```

3. **Definition of Gene Background:**
   - Use a reference set that includes all genes detected in the experiment.
```r
   max_ranking <- max(rankings[is.finite(rankings)])
   min_ranking <- min(rankings[is.finite(rankings)])
   rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
   rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
   rankings <- sort(rankings, decreasing = TRUE)
```

4. **Running GSEA:**
   - Compare the distribution of genes in a specific pathway to the global ranking.
   - Calculate a normalized enrichment score (NES) for each gene set.
```r
   GSEA_results <- fgsea(
     pathways = bg_genes,  # Gene sets to evaluate
     stats = rankings,
     scoreType = 'std',  # Ranking type (std: positive and negative values)
     minSize = 10,
     maxSize = 500,
     nproc = 1  # Number of parallel threads
   )
```

5. **Result Visualization:**
   - Generate a table of affected biological pathways.
```r
# Define the number of enriched pathways (positive and negative) to include.
number_of_top_pathways_up = 10
number_of_top_pathways_down = 10

# Select the most relevant enriched pathways, both positive and negative.
topPathwaysUp <- GSEA_results[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEA_results[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]

# Combine positive and negative enriched pathways.
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# Save a plot with the selected enriched pathways.
pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
dev.off()
```

---

## Implementation in R

This analysis was conducted in R using the `fgsea` package. The steps are based on a tutorial from BiostatSquid [(link to tutorial)](https://biostatsquid.com/fgsea-tutorial-gsea/). The complete code, detailing all the necessary steps for gene set enrichment analysis using differential expression data and functional databases, is included in the following script:
- [GSEA.R](./Recursos/GSEA.R)

In this file, all the previously described steps are included. A visualization of the top 30 significantly enriched pathways based on the project's data is provided at this [link](./Recursos/GO.RDS_gsea_top30pathways.pdf).

I want to express my gratitude to Laura, the creator of BiostatSquid, for her dedication to creating these tutorials. Her effort to teach and explain complex concepts in a simplified manner is immensely helpful for newly graduated bioinformaticians.

---
