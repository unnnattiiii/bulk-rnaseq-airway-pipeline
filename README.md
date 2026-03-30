# Bulk RNA-seq Pipeline: Dexamethasone Treatment in Airway Smooth Muscle Cells

## Overview
End-to-end bulk RNA-seq differential expression analysis of dexamethasone-treated airway smooth muscle cells using dataset GSE52778.

**Biological Question:** How does glucocorticoid treatment transcriptionally reprogram airway smooth muscle cells? Dexamethasone is a frontline corticosteroid therapy for asthma, but the full transcriptional mechanism in airway tissue is still being characterized.

**Key Result:** 823 significant differentially expressed genes (padj < 0.05, |log2FC| > 1), with clear separation between treated and untreated samples on PCA, and top hits consistent with known glucocorticoid response programs.

---

## Pipeline
```
Raw FASTQ → FastQC → MultiQC → fastp → HISAT2 → featureCounts → DESeq2
```

### Tools and Versions
- FastQC v0.12.1
- MultiQC v1.33
- fastp v1.3.0
- HISAT2 v2.2.2
- SAMtools v1.23.1
- featureCounts (Subread v2.1.1)
- DESeq2 v1.50.2
- R v4.5.2

### Reference
- Genome: GRCh38 chromosome 1 (Ensembl release 109)
- Annotation: Homo_sapiens.GRCh38.109.gtf

---

## Repository Structure
```
bulk-rnaseq-airway-pipeline/
├── qc/                         # FastQC and fastp QC reports
├── counts/                     # featureCounts output
├── results/
│   └── deseq2/
│       ├── deseq2_analysis.R   # DESeq2 analysis script
│       ├── significant_DEGs.csv
│       ├── all_results.csv
│       ├── pca_plot.pdf
│       ├── volcano_plot.pdf
│       └── heatmap.pdf
├── logs/                       # Alignment logs
├── .gitignore
└── README.md
```

---

## Methods Summary

**Quality Control:** Raw reads assessed with FastQC and MultiQC. All samples showed high per-base quality scores (Phred >33) and consistent GC content (~50%).

**Trimming:** Adapter sequences and low-quality bases removed using fastp. Average 97% of reads retained after filtering.

**Alignment:** Trimmed reads aligned to GRCh38 chromosome 1 using HISAT2 with splice-aware alignment. Note: chromosome 1 subset was used due to local compute constraints; full genome analysis requires HPC infrastructure.

**Quantification:** Read pairs counted per gene using featureCounts with reverse strand specificity (-s 2), consistent with dUTP library preparation protocol.

**Differential Expression:** DESeq2 analysis performed on the full 8-sample airway dataset (Bioconductor airway package) with design formula accounting for cell line and treatment effects (~ cell + dex).

---

## Results

| Metric | Value |
|--------|-------|
| Total genes tested | 22,369 |
| Significant DEGs | 823 |
| Upregulated | 454 |
| Downregulated | 369 |
| padj threshold | < 0.05 |
| LFC threshold | > 1 |

---

## Limitations
- Alignment performed on chromosome 1 only due to local disk space constraints (50GB available vs ~150GB required for full genome)
- Full genome analysis with all 8 samples through the complete upstream pipeline requires HPC or cloud infrastructure
- DESeq2 analysis used the Bioconductor airway package count matrix for statistical validity (n=4 per group)

---

## Dataset
- GEO Accession: GSE52778
- Organism: Homo sapiens
- Tissue: Airway smooth muscle cells
- Condition: Dexamethasone treatment vs untreated
- Samples: 8 total (4 treated, 4 untreated)
- Platform: Illumina HiSeq 2000, paired-end

---

## Author
Unnati Moradiya
MS Bioinformatics, Northeastern University
