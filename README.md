**Analysis of Proteomic and Transcriptomic Data**

**A. Required Environment and Software:**
1. Operating System: Ubuntu Version 18.04
    - **Software**: fastqc 0.12.1, hisat2 2.2.1, samtools 1.17, stringtie 2.2.1, trimmomatic 0.39
2. R Version 4.3.2
    - **CRAN Packages**: BiocManager 1.30.22, dplyr 1.1.4, matrixStats 1.2.0, pROC 1.18.5
    - **GitHub Package**: tidy_estimate 1.1.1
    - **Bioconductor Packages**: Biobase 2.64.0, DEXSeq 1.50.0, IsoformSwitchAnalyzeR 2.4.0, MatrixGenerics 1.16.0, SummarizedExperiment 1.34.0, edgeR 4.2.0, limma 3.60.2, rtracklayer 1.64.0

**B. How to Use the Code:**
1. Ensure you have the specified environment and required packages/software.
2. Execute the following scripts:
    - **RNA-seq-data-manipulation_linux.txt**: Processes RNAseq data (fastq files) for quality checking, read trimming, alignment, transcript assembly, and gene expression estimation (Ubuntu environment).
    - **RNA-differential-expression-analysis.R**: Creates gene expression tables and performs differential gene expression analysis (R environment).
    - **Protein-differential-expression-analysis.R**: Creates protein expression tables and performs differential protein expression analysis (R environment).
    - **Proteo-transcriptomic-analysis.R**: Integrates proteomic and transcriptomic data, including protein and mRNA annotation, urine-tissue mRNA correlation, and identification of genes with significant expression at both RNA and protein levels (R environment).
    - **Tidyestimate-(Tumor purity).R**: Estimates tumor purity in urinary cells and tissue using the RNA-based inference tool, ESTIMATE (R environment).
    - **Validation-of-the-result.R**: Performs validation analysis (R environment).
