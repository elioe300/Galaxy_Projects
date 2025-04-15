# Introduction

This project aims to delve into the fundamental aspects of **RNA-seq** and best practices for transcriptomic analysis using the Galaxy platform. Through this analysis, the objective is to **develop advanced knowledge on differential expression analysis**, addressing key topics such as:

- Creating graphs to evaluate the sample status prior to analysis.
- Visualizing the effects of filtering low-expression genes using specific graphs.
- Using graphs like **MDS** and **barplots** to explore the distribution of p-values and other relevant parameters.

Additionally, the project experiments with the **pseudomapping** method using the **Salmon** tool as part of the workflow. A basic differential gene expression (DGE) analysis employing **R** is also included.

This project is based on data from the study [*Long and short-read transcriptome profiling of human lung cancer cell lines*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545141/), providing a practical context to illustrate each addressed aspect.

---

## Sample Information

The study utilizes transcriptomic data generated in the framework of the study titled [*Long and short-read transcriptome profiling of human lung cancer cell lines*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545141/) by Dong et al. (2020). This study analyzed the transcriptomic profiles of **HCC827** and **H1975** cell lines, which exhibit different mutations in the **EGFR** gene.

The data includes paired-end **RNA-seq** reads generated using **Illumina™** technology. These datasets are available in the SRA repository:

- **H1975**:
  - [SRR14286057](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286057)
  - [SRR14286058](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286058)
  
- **HCC_827**:
  - [SRR14286065](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286065)
  - [SRR14286066](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286066)

---

# Methodology

## Quality Control Analysis

To keep this README concise, note that the reads exhibit good quality. GC values are within the expected range for RNA-seq data, duplicates are acceptable, adapter content is negligible, and the number of sequences is adequate for differential expression analysis.

Related quality control images can be viewed in the following link: [Quality Control](./Control_Calidad)

---

## Pseudomapping

The default settings of the **Salmon** tool were used, adjusting the **Kmer length** parameter to a value of 31. This decision was based on the recommendations in the [official Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode), which states: _We find that a k of 31 seems to work well for reads of 75bp or longer_. In this case, the reads are 80bp in length, fulfilling this condition.

The analysis was carried out using the transcriptome from **GENCODE v47**, which is available at the following link: [GENCODE v47](https://www.gencodegenes.org/human/).

Once results were generated with **Salmon**, further processing was done in **R** to perform additional analyses, including exploratory and differential gene expression analysis.

---

## Data Exploration

During data exploration, sample names and their experimental groups were defined. Salmon quantification files were located, and a transcript-to-gene map based on GENCODE v47 was generated. The data was imported using tximport to obtain a count matrix, which was cleaned and renamed. Finally, genes were annotated using org.Hs.eg.db, removing duplicates and genes without valid annotations. The data is now ready for DGE analysis.

---

### Preliminary Graph

To explore the initial characteristics of RNA-seq data before applying differential analyses or methodological corrections, a violin plot was used. On the **y-axis**, expression values are displayed in **log2 Counts per Million (logCPM)**, while the **x-axis** shows the different samples: **H1975_1**, **H1975_2**, **HCC827_1**, and **HCC827_2**.

Each violin illustrates the density of expression data for each sample, providing an idea of gene distribution. The black line in the center of each violin indicates the **median** expression level, serving as a reference for central levels in each sample.

<p align="center">
  <img src="./Recursos/LogCPMnofiltnonTMM.png">
</p>

A significant accumulation of genes with low expression below the median is observed. These genes often introduce noise into the analysis and do not provide biologically relevant information. The **edgeR** package offers a `filterbyExpr` option to eliminate low-expression genes and ensure that only biologically active and relevant ones are included in differential expression analysis.

---

### Graph After Applying `filterByExpr`

After applying the `filterByExpr` filter, significant changes in the RNA-seq data distribution are observed, as reflected in the graph:

<p align="center">
  <img src="./Recursos/LogCPMfilterednonTMM.png">
</p>

Genes with low expression, which generate noise and affect interpretation, have been removed. This results in cleaner and more homogeneous distributions. After gene filtering, the median is adjusted to logCPM expression values where genes are no longer noisy, providing a suitable reference for sample comparisons.
While the data has been filtered, normalization by TMM is necessary to correct technical differences between samples, ensuring that observed variations are of biological origin.

---

### Plot After Applying TMM

After applying normalization using TMM (Trimmed Mean of M-values), RNA-seq data exhibits an adjusted and comparable distribution across samples, reflected in the graph:

<p align="center">
  <img src="./Recursos/LogCPMfilteredTMM.png" alt="TMM Normalized Data">
</p>

The logCPM distributions are now more consistent between samples, indicating that technical differences, such as sequencing depth or library composition, have been corrected. The black lines representing the medians of the distributions are now aligned across samples, showing that normalization has balanced expression levels for reliable comparisons.

Thanks to TMM normalization, non-biological variations between samples have been mitigated, ensuring that observed differences are predominantly biological.

*(_Acknowledgment to Daniel Beiting for his comprehensive explanation of transcriptomic analysis on his website [DIYtranscriptomics](https://diytranscriptomics.com/), part of his scripts were used to visualize all previously mentioned steps_)*

---

### What is an MDS Plot?

The MDS plot (*Multidimensional Scaling Plot*) is a visualization tool used to graphically represent relationships among samples in a lower-dimensional space, typically two dimensions. It is often performed before proceeding with differential analyses, as it helps confirm whether observed clustering matches the expectations of the experimental design. It allows us to detect whether sample groupings reflect biological similarities or whether outliers exist, caused by technical errors or significant differences in underlying biology.

---

### MDS Plot Analysis for the 4 Samples

The MDS plot generated for **H1975_1**, **H1975_2**, **HCC827_1**, and **HCC827_2** shows a clear separation into two groups:

<p align="center">
  <img src="./Recursos/PlotMDS.png" alt="MDS Plot for 4 Samples">
</p>

- **H1975_1** and **H1975_2** samples are close in the MDS space, indicating that they have highly similar gene expression profiles.
- **HCC827_1** and **HCC827_2** samples form a separate group, reflecting the consistency between the biological replicates of the HCC827 cell line.

The MDS plot confirms that observed differences among samples are consistent with experimental expectations: biological replicates of each cell line are similar to each other, but the H1975 and HCC827 cell lines exhibit distinctly different expression profiles. This validates the experimental design and suggests that the data is reliable for differential expression analysis.

---

### Checking p-Value Distribution

Another recommended step is verifying the distribution of p-values. This is relevant, especially as we are conducting a statistical test with hundreds, thousands, or even millions of p-values. By reviewing the histogram, we can immediately identify whether the statistical test has captured significant differences through many very small p-values or, conversely, if the distribution is entirely uniform, which could indicate a lack of signal or issues in the analysis. This graph provides an immediate overview of the performance of tests across all hypotheses and facilitates the detection of potential problems.

Below are six approximate versions of what the histogram might look like:

<p align="center">
  <img src="./Recursos/Ejemplos_pvalores.png" alt="p-Value Distribution Examples">
</p>

*Image taken from [How to interpret a p-value histogram](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/) by David Robinson.*

In our data, the p-value histogram shows a distribution very similar to example **A**, in which there is a high frequency of small p-values that quickly decay toward larger values.

<p align="center">
  <img src="./Recursos/DistribucionFDR.png" alt="p-Value Histogram for FDR Distribution">
</p>

This suggests that:
- Many tests result in small p-values, indicating significant evidence against the null hypothesis in numerous cases.
- The shape of the histogram is as expected under a scenario of genuine discoveries, providing confidence in the robustness of analyses performed prior to multiple testing correction.

*(_Acknowledgment to Ming ‘Tommy’ Tang for the wonderful explanation regarding [p-values](https://divingintogeneticsandgenomics.com/post/understanding-p-value-multiple-comparisons-fdr-and-q-value/) and [the importance of reviewing histograms](https://divingintogeneticsandgenomics.com/post/downstream-of-bulk-rnaseq-read-in-salmon-output-using-tximport-and-then-deseq2/)_)*

---

### Statistical Test

After completing all checkpoints and verifying the implementation of procedures to minimize technical variability and retain only relevant genes, the statistical test is conducted. The Benjamini-Hochberg method is employed to control the false discovery rate, with a p-value threshold of 0.05 and a logFoldChange of 1, the default value.
```r
is.de1 <- decideTests(result, adjust.method = "BH", p.value = 0.05, lfc=1)
summary(is.de1)
```

<p align="center">
  <img src="./Recursos/Resultado_Test.png" alt="Result Summary of Test">
</p>

The complete script can be accessed [here](./Recursos/DEG_Salmon.R)


