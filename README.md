# Analysis Code for: HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications

---

## Description

This repository contains all analysis code accompanying the HaemaTONIC manuscript. HaemaTONIC is a novel multilineage human haematopoiesis in vitro culture platform that supports the differentiation and functional characterisation of diverse blood cell lineages.

The code covers:
- Multimodal single-cell RNA-seq analysis (10x Chromium v2 and 10x Flex with RNA + ADT + HTO)
- Flow cytometry UMAP and Phenograph clustering
- Drug screen phenotypic analysis
- CRISPR knockout phenotypic analysis
- Clonal haematopoiesis (base-edit) mutation analysis

---

## Data Availability

Raw sequencing data (FASTQ files) and processed count matrices will be deposited at NCBI GEO.

Flow cytometry input data (Excel files) required for Figures 2, 4, 5, and 6 are included in Zenodo (https://doi.org/10.5281/zenodo.20756620) 

---

## Repository Structure

```
HaemaTONIC_code/
├── README.md
├── sessionInfo.txt                     ← R session info for reproducibility
│
├── 01_10Xv2_Analysis.R                 ← 10x v2 multimodal pipeline (Fig 1 / Extended)
├── 02_Figure_2.R                       ← Flow cytometry UMAP and clustering
├── 03_Figure_3_10X_flex.R              ← 10x Flex multiomics + stimulation DE
├── 04_Figure_4.R                       ← Drug screen analysis
├── 05_Figure_5.R                       ← CRISPR KO phenotypic analysis
└── 06_Figure_6.R                       ← Clonal haematopoiesis mutation analysis
```

---

## Software Requirements

All analysis was performed in R. Key packages and versions:

| Package | Version | Source |
|---|---|---|
| R | 4.3.x | https://www.r-project.org |
| Seurat | 5.x | CRAN |
| SeuratObject | 5.x | CRAN |
| SeuratDisk | 0.0.0.9021 | GitHub (mojaveazure/seurat-disk) |
| SeuratData | 0.2.2 | GitHub (satijalab/seurat-data) |
| AzimuthAPI | dev | GitHub (satijalab/AzimuthAPI) |
| SingleR | 2.x | Bioconductor |
| celldex | 1.x | Bioconductor |
| SingleCellExperiment | 1.x | Bioconductor |
| dplyr | 1.x | CRAN |
| ggplot2 | 3.x | CRAN |
| patchwork | 1.x | CRAN |
| reshape2 | 1.x | CRAN |
| tidyr | 1.x | CRAN |
| readxl | 1.x | CRAN |
| uwot | 0.x | CRAN |
| Rphenograph | 0.x | GitHub (JinmiaoChenLab/Rphenograph) |
| viridis | 0.x | CRAN |
| paletteer | 1.x | CRAN |
| scales | 1.x | CRAN |
| ggbreak | 0.x | CRAN |
| ggrepel | 0.x | CRAN |
| DoMultiBarHeatmap | dev | GitHub |
| hdf5r | 1.x | CRAN |
| stringr | 1.x | CRAN |
| rlang | 1.x | CRAN |

Full session info (R version, package versions, platform) is recorded in `sessionInfo.txt`.

### Installation

```r
# CRAN packages
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "reshape2", "tidyr",
  "readxl", "uwot", "viridis", "paletteer", "scales",
  "ggbreak", "ggrepel", "hdf5r", "stringr", "rlang"
))

# GitHub packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("satijalab/AzimuthAPI")
remotes::install_github("JinmiaoChenLab/Rphenograph")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleR", "celldex", "SingleCellExperiment",
  "GenomeInfoDb", "BSgenome", "ensembldb"
))
```

---

## How to Run

### Step 1 — Download the data

Download the Flow cytometry data or seurat objects from Zenodo (https://doi.org/10.5281/zenodo.20756620) and place them in a local directory.

### Step 2 — Update input paths

Each script contains an **INPUT PATHS** section near the top. Update every `"path/to/..."` placeholder to point to your local files before running.

For example, in `01_10Xv2_Analysis.R`:
```r
PATH_10X_ST158 <- "path/to/ST158/outs/filtered_feature_bc_matrix"
PATH_HTO       <- "path/to/ST158_HTO_counts_reseq.txt"
# etc.
```

### Step 3 — Run scripts in order

Scripts numbers refer to the Figure where the graph are displayed. Each script saves its outputs to a numbered subfolder under `output/`:

```
output/
├── 01_10Xv2/
├── 02_Figure2/
├── 03_Figure3/
├── 04_Figure4/
├── 05_Figure5/
└── 06_Figure6/
```

Scripts are independent of each other and can be run in any order.

### Step 4 — Notes on large intermediate files

- Script 01 produces `pbmc.singlet_annotated.rds` — this is the annotated Seurat object for the 10x v2 dataset. This file is available for download from Zenodo (https://doi.org/10.5281/zenodo.20756620).
- Script 03 produces `combined.singlets_annotated.rds` — the annotated Seurat object for the 10x Flex dataset, also available from Zenodo.
- If you download these pre-computed objects from Zenodo, you can skip the relevant processing steps in Scripts 01 and 03 and go directly to the visualisation sections.

---

## Notes on HTO Demultiplexing (Script 03)

Script 03 exports a raw metadata table (`metadata_raw.txt`) after HTO demultiplexing. This file must be manually annotated to resolve which hashtag corresponds to which donor/condition combination (this information is specific to the experimental design described in the paper). The annotated file is then re-imported as `metadata_resolved.txt`.

A pre-annotated version of `metadata_resolved.txt` is included in the Zenodo data folder.

---

## Contact

Sara Tomei
Walter and Eliza Hall Institute of Medical Research (WEHI)
Melbourne, Australia
tomei.s@wehi.edu.au


