# cancer_transcriptomics

# Overview
This project performs a comprehensive bioinformatics analysis of Kidney Cancer (TCGA) gene expression data. It demonstrates a complete R workflow from raw data processing to machine learning classification.

Key Analyses:
- Data QC & Normalization: Using edgeR and limma-voom.
- Dimensionality Reduction: PCA to visualize cancer subtypes.
- Survival Analysis: Kaplan-Meier curves for patient prognosis.
- Classification: A Random Forest model to predict cancer subtypes (KIRC, KICH, KIRP) based on gene signatures.

# Prerequisites 
To run this analysis, you need R installed.

Required R Packages
Copy and run this code in your R console to install all dependencies:

    # 1. Install CRAN packages
    install.packages(c("here", "tidyverse", "survival", "survminer", "pROC", "randomForest", "caret", "cowplot", "reshape2"))

    # 2. Install Bioconductor packages
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("biomaRt", "limma", "edgeR", "DESeq2", "tximport", "rhdf5"))


# Data Setup
1. Download the Data: https://drive.google.com/drive/folders/1yaz2O69_CDS9gxNJe-JiJaUH_i86A4_Q?usp=drive_link

2. Organize Your Folders: After cloning this repository, create a data folder and place the downloaded files inside as shown below:

```
/omics-kidney-analysis
  ├── omics_analysis.R       # The main script
  └── data/                  # Create this folder
       ├── TCGA_KIDNEY_Gene_ReadCounts.txt
       ├── TCGA_KIDNEY_ClinicalAnnotation.txt
       ├── VEGF_A_UP.V1_UP.v2024.1.Hs.txt
       ├── VEGF_A_UP.V1_DN.v2024.1.Hs.txt
       └── output/           # Create this sub-folder
            └── abundance.h5
```

# How to Run
1. Clone this repository.
2. Set up the data folder as shown above.
3. Open omics_analysis.R in RStudio.
4. Run the script! (The script automatically handles file paths).
