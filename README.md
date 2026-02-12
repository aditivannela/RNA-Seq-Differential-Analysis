# RNA-Seq-Differential-Analysis
This project performs of differential analysis of RNA-sequencing data from mouse bone marrow derived macrophages pretreated with a KAT8 inhibitor and stimulated with R848. The purpose of this analysis is to identify genes that are expressed or repressed when KAT8 is inhibited during TLR7/8 stimulation.  
## KAT8 Inhibition in Mouse Bone Marrow–Derived Macrophages

---

## Project Overview

This project performs differential gene expression analysis using RNA-sequencing (RNA-seq) data from mouse bone marrow–derived macrophages (BMDMs) pretreated with a KAT8 inhibitor and stimulated with R848.
The goal of this analysis is to identify genes that are activated or repressed following KAT8 inhibition during Toll-like receptor (TLR7/8) stimulation.
Differential expression analysis was conducted using the DESeq2 package in R.

---

## Biological Background

KAT8 is a key epigenetic regulator involved in chromatin remodeling and transcriptional regulation that promotes inflammation and oxidative stress. Inhibition of KAT8 may alter macrophage inflammatory responses by modifying gene expression programs.

---

## Experimental Design

| Condition | Treatment | Replicates |
|-----------|-----------|------------|
| Control   | DMSO + R848 | 3 |
| Treatment | KAT8 inhibitor + R848 | 3 |

Total samples: 6

---

## Data Source

The RNA-seq dataset analyzed in this repository was obtained from the Gene Expression Omnibus (GEO):
- **Accession Number:** GSE263150  
- **Organism:** *Mus musculus*  
- **Sequencing Platform:** Illumina RNA-seq  
This repository contains a secondary analysis of publicly available data. All raw sequencing data and experimental details are available through GEO.

---

## Bioinformatics Workflow

Analysis was performed in R using DESeq2:
1. Import gene-level count matrix
2. Create metadata for experimental conditions
3. Filter low-count genes
4. Perform normalization and differential expression analysis
5. Shrink log2 fold changes using apeglm
6. Identify significant genes (padj < 0.05, |log2FC| > 1)
7. Visualization:
   - PCA plot
   - Volcano plot
   - Heatmap of top 50 differentially expressed genes
8. Export full and significant results tables

---

## Output Files
- `DESeq2_all_results_KAT8i_vs_DMSO.csv`
- `DESeq2_significant_genes_padj0.05_LFC1.csv`
- `PCA_plot.pdf`
- `Volcano_plot.pdf`
- `Heatmap_top_genes.pdf`

---

## Key Findings
Differential expression analysis identified thousands of significantly regulated genes following KAT8 inhibition, indicating that KAT8 plays a substantial role in regulating macrophage transcriptional responses to R848 stimulation.
These results suggest that chromatin acetylation via KAT8 is an important regulatory mechanism in inflammatory gene expression.

---

## Requirements

Install required packages before running the script:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
install.packages(c("ggplot2", "pheatmap"))
