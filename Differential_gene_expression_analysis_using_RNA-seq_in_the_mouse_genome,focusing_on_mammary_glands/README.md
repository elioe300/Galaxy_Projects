# Transcriptomic Analysis of Mammary Glands in Mice

## Introduction

This project adheres to the fundamental practices for transcriptomic analysis using the Galaxy platform. The analyzed samples originate from mammary tissue of lactating and pregnant mice. The main goal of this analysis is to understand molecular mechanisms associated with lactation and pregnancy through the identification of gene expression patterns.

---

## Sample Information

The analyzed samples are derived from mammary tissue of mice featured in the study titled [_EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival_](https://www.nature.com/articles/ncb3117) by Fu et al. (2015). This study examined expression profiles of basal and luminal cells in the mammary glands of virgin, pregnant, and lactating mice. Specifically, transcriptomic data from basal cells in mammary tissue of lactating and pregnant mice were used. These were generated as *single-end* reads using the Illumina HiSeq 2000 technology. The datasets are available in the SRA repository:

- **Basal pregnant**:
  - [SRR1552452](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552452)
  - [SRR1552453](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552453)

- **Basal lactate**:
  - [SRR1552454](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552454)
  - [SRR1552455](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552455)

---

# Methodology

## Quality Control Analysis

### Basic Statistics

To evaluate the quality of the generated reads, the widely recognized **FastQC** tool was used for transcriptomic analysis. Initially, individual analyses were conducted for each sample, followed by a consolidated analysis using **MultiQC**, which enabled comparative and general insights into the raw RNA-seq data.

<p align="center">
  <img src="./Recursos/FastQC_general_statistics.png" alt="FastQC General Statistics">
</p>

#### 1. **GC Content**  
All samples exhibit a GC content of **50.0%**, which falls within the typical and acceptable range for RNA-seq analysis.

#### 2. **Number of Sequences**  
The total reads per sample are as follows:
- **SRR1552452**: 31.7M
- **SRR1552453**: 29.6M
- **SRR1552454**: 27.2M
- **SRR1552455**: 25.4M  

These figures comply with RNA-seq recommendations. According to Illumina guidelines, in gene expression experiments, [_Gene expression profiling experiments that are looking for a quick snapshot of highly expressed genes may only need 5 million to 25 million reads per sample_](https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000001243).

#### 3. **Duplicates**  
The observed duplicate percentages are somewhat high, but this is common in RNA-seq data due to the high expression of certain genes. Although this does not invalidate the analysis, it is advisable to confirm that the duplicates originate from real biological transcripts rather than technical artifacts. To validate this, the relationship between duplicate percentages and GC content was analyzed.  
<p align="center"> <img src="./Recursos/FastQC_GC_content.png" alt="GC Content"> </p>

The GC content analysis shows that reads follow a modal distribution, without sharp peaks or double-modal peaks. This supports the hypothesis that duplicates are biological rather than due to contamination or technical errors. This is further validated by the following image, which shows minimal adapter content presence in the sequences:  
<p align="center"> <img src="./Recursos/FastQC_Adapter_content.png" alt="Adapter Content"> </p>

Additionally, it is generally recommended not to remove duplicates when working with RNA-seq data unless unique molecular identifiers are used, as indicated in the study [(Klepikova et al. 2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5357343/).

---

### Per-Base Sequence Quality

According to the study titled [_A survey of best practices for RNA-seq data analysis_](https://pmc.ncbi.nlm.nih.gov/articles/PMC4728800/#Sec3), *"as a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability."* However, the raw read quality is adequate and does not exhibit significant drops toward the 3’ end. Therefore, quality control on raw reads is not necessary.  
<p align="center"> <img src="./Recursos/FastQC_quality_scores.png" alt="Per-Base Sequence Quality Graph"> </p>

---

### Sequence Content per Base

In the **Per base sequence content** section, a **deregulation in nucleotide percentages** was detected in the initial positions.
<p align="center"> <img src="./Recursos/FastQC_GC_sequence_content.png" alt="GC Content by Base Graph"> </p>

As noted by the **Center for Cancer Research (CCR)** of the **National Cancer Institute (NCI)** in their guide [_Bioinformatics for Beginners 2022_](https://bioinformatics.ccr.cancer.gov/docs/b4b/Module2_RNA_Sequencing/Lesson10/), the variability observed in the first 10-15 bases of reads is a normal phenomenon in RNA-seq experiments. This fluctuation can be attributed to the use of *random primers* during library preparation, which generates a biased start in the reads.

---

## Read Mapping

The mapping process was conducted using **HISAT2**, a widely utilized tool for read alignment in transcriptomic analyses, using the mouse genome (**GRCh39/mm39**) as a reference. Default options were applied, with the `-k` parameter set to **5**. In Galaxy, the mapped reads are sorted by coordinates by default. Below are the general alignment statistics, generated by **MultiQC**:
<p align="center">
  <img src="./Recursos/Mapping_quality_scores.png" alt="Alignment Quality Statistics from MultiQC">
</p>

These results reflect high alignment quality, as all samples far exceed the threshold of **75% uniquely aligned reads**, which is considered an indicator of good quality in RNA-seq experiments, as established in the [Harvard Chan Repository](https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon/blob/master/lessons/05_multiQC.md#assessing-the-quality-control-metrics).

---

## Quality Control of Read Mapping

### 5'-3' Bias

According to best practices described in the same repository, the 5'-3' bias was analyzed. A value close to **1** was obtained, indicating a uniform distribution of reads along the transcripts.
<p align="center">
  <img src="./Recursos/Qualimap_5'3'_bias.png" alt="Qualimap 5'3' Bias">
</p>

---

### Junction Analysis

The junction analysis indicates that **81.7%** of the identified junctions correspond to known junctions. This high percentage demonstrates that most reads align with previously annotated regions in the reference genome, which is ideal for this study as it does not aim to identify novel transcripts.
<p align="center">
  <img src="./Recursos/Qualimap_Junction_Analysis.png" alt="Qualimap Junction Analysis">
</p>

---

### Genomic Origin of Reads

The [Harvard Chan Repository](https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon/blob/master/lessons/05_multiQC.md#assessing-the-quality-control-metrics) indicates that _"Generally in a good library, we expect over 60% of reads to map to exons for mouse and human organisms."_ In this analysis, **73%** of the reads aligned with exonic regions, exceeding the expected threshold. This confirms that the reads primarily originate from mRNA.

<p align="center">
  <img src="./Recursos/Qualimap_Reads_Genomic_Origin.png" alt="Qualimap Reads Genomic Origin">
</p>

---

## Read Counting

Read counting was conducted using Galaxy's featureCounts tool to process the BAM file generated by **HISAT2**. The tool's configuration options were as follows:
- **Specify strand information (option -s)**: **Unstranded**
- **Gene annotation file**: **Mus_musculus.GRCm39.113.gtf**
- **Feature type (option -t)**: **exon**
- **Gene ID attribute (option -g)**: **gene_id**

---

## Differential Expression Analysis (DEG)

The complete script can be found [here](./Recursos/DEG_analysis.R).

### Methodology

- *featureCounts* data was used, obtained from alignments to the mm39 (GRCm39) reference genome.
- The **biomaRt** package was used for gene annotation, mapping *Ensembl Gene IDs* to *Entrez IDs*, as the latter are required by the **goana** package for Gene Ontology analysis.
- The **edgeR** package was used for data normalization and identification of differentially expressed genes.
- A selection criterion of **logFC ≥ 2** or **logFC ≤ -2** with **FDR < 0.05** was established to define significantly regulated genes.

---

### Results

#### 1. Number of Differentially Expressed Genes

Below is a summary of the differentially expressed genes after adjusting *p-values* using the Benjamini-Hochberg (BH) method:

<p align="center">
  <img src="./Recursos/summary_BH_test.png" alt="Summary of BH Test">
</p>

---

#### 2. Volcano Plot

The following plot illustrates the distribution of differentially expressed genes based on their **log Fold Change (logFC)** and statistical significance (-log10(*FDR*)):

<p align="center">
  <img src="./Recursos/Rplot.png" alt="Volcano Plot">
</p>

Downregulated genes are represented in **red**, while upregulated genes appear in **blue**.

---

#### 3. Gene Ontology (GO)

The **Gene Ontology (GO)** analysis was conducted using the **goana** package.

Below are the five most representative categories:

<p align="center">
  <img src="./Recursos/top_GO.png" alt="Top GO Categories">
</p>

The analysis of differentially expressed genes identifies genes with significant expression changes between experimental conditions. To interpret these changes in a broader biological context, enrichment analyses such as Over-Representation Analysis and Gene Set Enrichment Analysis are employed. Using the links below, you can view the results obtained with these methodologies:
- [Over-Representation Analysis](./Over-Representation_Analysis)
- [Gene Set Enrichment Analysis](./Gene_Set_Enrichment_Analysis)

---
