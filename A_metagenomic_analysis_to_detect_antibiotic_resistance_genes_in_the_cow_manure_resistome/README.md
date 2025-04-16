# A metagenomic analysis to detect antibiotic resistance genes in the cow manure resistome

## Introduction

This project focuses on metagenomic analysis using the Galaxy platform to replicate the work described in the article [Aerobic composting as an effective cow manure management strategy for reducing the dissemination of antibiotic resistance genes: An integrated meta-omics study](https://www.sciencedirect.com/science/article/pii/S0304389419318497).

The main objective of this project is using the Galaxy platform to analyze the microbial diversity and detect **antibiotic resistance genes (ARGs)** utilizing metagenomic data, inspired by the previously mentioned paper. This approach will strengthen my skills in handling and analyzing this type of data.

Here’s the updated version of your sample information in Markdown:

---

## Sample Information

This project analyzes metagenomic samples from aerobic composting processes, divided into two distinct conditions: *Win* and *Sum*. Each sample corresponds to specific days of sampling, represented by numbers indicating the sampling day. The objective is to analyze microbial diversity and detect **antibiotic resistance genes (ARGs)**.

### Available Data in SRA
- **Win Condition:**
  - [Win 0 - SRR8719648](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719648)
  - [Win 3 - SRR8719647](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719647)
  - [Win 14 - SRR8719645](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719645)

- **Sum Condition:**
  - [Sum 0 - SRR8719644](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719644)
  - [Sum 3 - SRR8719643](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719643)
  - [Sum 14 - SRR8719642](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8719642)

---

# Methodology

## Quality Control Analysis

To keep this README concise, note that the reads exhibit good quality. GC values are within the expected range for RNA-seq data, duplicates are acceptable, adapter content is negligible, and the number of sequences is adequate for differential expression analysis.

Related quality control images can be viewed in the following link: [Quality Control](./Quality_Control)

---


## Assembly Process

For the assembly step, the procedures described in the mentioned study were followed. Each sample was assembled de novo using the tool **MEGAHIT** with the option `--min-contig-len` set to 300 bp, following the authors assembly process: _Multiple_Megahit was used to combine contigs with minimum lengths longer than 300 bp_. Additionally, by consulting the [MEGAHIT documentation titled "Assembly Tips"](https://github.com/voutcn/megahit/wiki/Assembly-Tips), the recommendations of the authors for these cases were followed, since the study did not specify the k-mers used in the samples.

Key recommendations adhered to include:
- **For ultra-complex metagenomics data (e.g., soil):** A larger `kmin`, such as 27, is recommended to reduce the complexity of the de Bruijn graph. **Quality trimming** is also recommended.
- **K-mer list:** The following k-mer sizes were used: `27, 37, 47, 57, 67, 77, 87, 97, 107, 117, 127`.
- **Smaller k-mer step:** A k-mer step of 10 was used, following the recommendation: “Smaller `--k-step`, say 10, is more friendly to low-coverage datasets.”
- **No Mercy k-mer option:** The `--no-mercy` option was used to recover low-coverage sequences, specially designed for metagenomics assembly. For generic datasets with coverage ≥ 30x, MEGAHIT may generate better results with this option.

