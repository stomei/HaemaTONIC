# =============================================================================
# Script: 03_Figure_3_10X_flex.R
# Description: 10x Flex multiomics (RNA + ADT + HTO) pipeline for Figure 3
#              Includes HTO demultiplexing, WNN dimensionality reduction,
#              cell type annotation (SingleR + Azimuth), differential expression
#              across immune stimulation conditions, per-cell-type DE analysis
# HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications
# Author: Sara Tomei, WEHI
# =============================================================================

# --- Libraries ----------------------------------------------------------------
# Note: SingleR, celldex, SingleCellExperiment, AzimuthAPI must be installed
# via Bioconductor (see installation block below) before first use
#
# BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment",
#                         "GenomeInfoDb", "BSgenome", "ensembldb"))
# remotes::install_github("satijalab/AzimuthAPI")

library(Seurat)
library(SeuratObject)
library(dplyr)
library(hdf5r)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(AzimuthAPI)
library(ggplot2)
library(tidyr)


# --- INPUT PATHS (update these to your local paths) --------------------------
# 10x Flex per-sample H5 files (CellRanger multi output)
POOL_PATHS <- list(
  Pool1 = "path/to/per_sample_outs/Pool1/sample_filtered_feature_bc_matrix.h5",
  Pool2 = "path/to/per_sample_outs/Pool2/sample_filtered_feature_bc_matrix.h5",
  Pool3 = "path/to/per_sample_outs/Pool3/sample_filtered_feature_bc_matrix.h5",
  Pool4 = "path/to/per_sample_outs/Pool4/sample_filtered_feature_bc_matrix.h5",
  Pool5 = "path/to/per_sample_outs/Pool5/sample_filtered_feature_bc_matrix.h5",
  Pool6 = "path/to/per_sample_outs/Pool6/sample_filtered_feature_bc_matrix.h5"
)

# Metadata file with resolved HTO assignments (Donor + Condition columns)
# See README for format: this file is produced by manually resolving the
# HTO demultiplexing output (metadata.txt) and re-imported below
PATH_METADATA_RESOLVED <- "path/to/metadata_resolved.txt"

# Output directory
OUT_DIR <- "output/03_Figure3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# --- 1. HELPER: LOAD A SINGLE POOL -------------------------------------------
# Reads an H5 file, separates RNA / ADT / HTO, normalises, and returns a
# Seurat object with all three assays
load_flex_sample <- function(path, sample_name) {
  data     <- Read10X_h5(path)
  rna      <- data[["Gene Expression"]]
  antibody <- data[["Antibody Capture"]]

  # HTO features contain "hash" in the row name; everything else is ADT
  adt <- antibody[!grepl("hash", rownames(antibody)), ]
  hto <- antibody[ grepl("hash", rownames(antibody)), ]

  obj <- CreateSeuratObject(counts = rna, project = sample_name)
  obj[["ADT"]] <- CreateAssayObject(counts = adt)
  obj[["HTO"]] <- CreateAssayObject(counts = hto)

  obj <- NormalizeData(obj)
  obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR")

  obj$sample <- sample_name
  return(obj)
}


# --- 2. LOAD ALL POOLS AND MERGE ---------------------------------------------
pool_objects <- mapply(
  load_flex_sample,
  path        = POOL_PATHS,
  sample_name = names(POOL_PATHS),
  SIMPLIFY    = FALSE
)

combined <- merge(
  pool_objects[[1]],
  y           = pool_objects[-1],
  add.cell.ids = names(POOL_PATHS)
)
combined <- JoinLayers(combined)


# --- 3. QUALITY CONTROL ------------------------------------------------------
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
combined <- subset(combined, subset = nFeature_RNA > 200)


# --- 4. HTO DEMULTIPLEXING ---------------------------------------------------
combined <- NormalizeData(combined, assay = "HTO", normalization.method = "CLR")
combined <- HTODemux(combined, assay = "HTO", positive.quantile = 0.99)

table(combined$HTO_classification.global)
HTOHeatmap(combined, assay = "HTO", ncells = 5000)
Idents(combined) <- "HTO_maxID"
RidgePlot(combined, assay = "HTO",
          features = rownames(combined[["HTO"]])[1:4], ncol = 2)

# Export raw metadata for manual HTO resolution
metadata <- combined@meta.data
write.table(metadata, file = file.path(OUT_DIR, "metadata_raw.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)

# NOTE: After exporting metadata_raw.txt, manually assign Donor and Condition
# labels to each HTO hashtag combination and save as metadata_resolved.txt
# (columns: all original metadata columns + "resolved_sample" column)
# Then re-import:
metadata_resolved <- read.table(PATH_METADATA_RESOLVED,
                                header = TRUE, row.names = 1, sep = "\t")
combined@meta.data <- metadata_resolved


# --- 5. FILTER: RETAIN SINGLETS, SPLIT SAMPLE METADATA ----------------------
Idents(combined) <- "HTO_classification.global"
combined <- subset(combined, idents = "Negative", invert = TRUE)
combined.singlets <- subset(combined, idents = "Singlet")

# Split "Condition Donor" into two separate columns
metadata <- combined.singlets@meta.data
metadata <- separate(metadata, resolved_sample, c("Condition", "Donor"), sep = " ")
combined.singlets@meta.data <- metadata

# Remove PBMC hashtag (used only as demultiplexing control)
Idents(combined) <- "resolved_sample"
combined <- subset(combined, idents = "PBMC", invert = TRUE)

combined.singlets[["RNA"]] <- JoinLayers(combined.singlets[["RNA"]])


# --- 6. WNN DIMENSIONALITY REDUCTION -----------------------------------------
DefaultAssay(combined.singlets) <- "RNA"
combined.singlets <- NormalizeData(combined.singlets) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

DefaultAssay(combined.singlets) <- "ADT"
VariableFeatures(combined.singlets) <- rownames(combined.singlets[["ADT"]])
combined.singlets <- NormalizeData(combined.singlets,
                                   normalization.method = "CLR", margin = 2) %>%
  ScaleData() %>%
  RunPCA(reduction.name = "apca")

ElbowPlot(combined.singlets, reduction = "pca")  # use to decide dims

combined.singlets <- FindMultiModalNeighbors(
  combined.singlets,
  reduction.list       = list("pca", "apca"),
  dims.list            = list(1:20, 1:16),
  modality.weight.name = "RNA.weight"
)

combined.singlets <- RunUMAP(combined.singlets, nn.name = "weighted.nn",
                             reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combined.singlets <- RunUMAP(combined.singlets, reduction = "pca",  dims = 1:20,
                             assay = "RNA", reduction.name = "rna.umap",
                             reduction.key = "rnaUMAP_")
combined.singlets <- RunUMAP(combined.singlets, reduction = "apca", dims = 1:16,
                             assay = "ADT", reduction.name = "adt.umap",
                             reduction.key = "adtUMAP_")


# --- 7. CLUSTERING -----------------------------------------------------------
combined.singlets <- FindClusters(combined.singlets, graph.name = "wsnn",
                                  algorithm = 3, resolution = 1, verbose = FALSE)


# --- 8. CELL TYPE ANNOTATION — SingleR ---------------------------------------
sce.fig <- as.SingleCellExperiment(combined.singlets, assay = "RNA")
ref_singler <- celldex::fetchReference("novershtern_hematopoietic", "2024-02-26")
pred.main <- SingleR(test = sce.fig, ref = ref_singler, labels = ref_singler$label.main)

combined.singlets[["SingleR.labels"]] <- pred.main$labels
plotScoreHeatmap(pred.main,
                 clusters     = combined.singlets$seurat_clusters,
                 order.by.clusters = TRUE,
                 filename     = file.path(OUT_DIR, "heatmap_SingleR.png"),
                 width = 7, height = 8)


# --- 9. CELL TYPE ANNOTATION — Azimuth (Cloud API) ---------------------------
combined.singlets <- CloudAzimuth(combined.singlets)
combined.singlets <- RunUMAP(combined.singlets, dims = 1:128,
                             reduction = "azimuth_embed",
                             reduction.name = "azimuth_umap")
combined.singlets <- PrepLabel(combined.singlets, "azimuth_fine",
                               "azimuth_fine_filtered", cutoff = 10)

p_az <- DimPlot(combined.singlets, group.by = "azimuth_fine_filtered",
                label = TRUE, label.size = 3, reduction = "azimuth_umap") + NoLegend()
p_az


# --- 10. MANUAL CLUSTER ANNOTATION ------------------------------------------
combined.singlets$annotation <- plyr::mapvalues(
  combined.singlets$seurat_clusters,
  from = c("1","22","18","34","36","25","29","6","24","8","9","28",
           "0","11","10","33","35","20","4","31","3","15","13","7",
           "32","26","17","2","16","19","27","12","14","21","5","23","30"),
  to   = c("cDC1","cDC1","cDC1","cDC1","DC progenitor","cDC1","B cells",
           "cDC2","DC progenitor","cDC2","mregDC","NK cell",
           "Mast cell","Mast cell","Mast cell","Mast cell","Mast cell",
           "Basophils","HSPC","NK cell","Mix-lymphoid","Monocytes",
           "Monocytes","Macrophages","Macrophages","Granulocyte progenitor",
           "Eosinophils","Baso/Eos progenitor","Monocytes","Macrophages",
           "Macrophages","Neutrophils","pDC","pDC",
           "Stromal contaminant","Stromal contaminant","Stromal contaminant")
)

combined.singlets$annotation_generic <- plyr::mapvalues(
  combined.singlets$seurat_clusters,
  from = c("1","22","18","34","36","25","29","6","24","8","9","28",
           "0","11","10","33","35","20","4","31","3","15","13","7",
           "32","26","17","2","16","19","27","12","14","21","5","23","30"),
  to   = c("Dendritic cell","Dendritic cell","Dendritic cell","Dendritic cell",
           "Progenitor","Dendritic cell","Lymphoid","Dendritic cell",
           "Progenitor","Dendritic cell","Dendritic cell","Lymphoid",
           "Granulocyte","Granulocyte","Granulocyte","Granulocyte","Granulocyte",
           "Granulocyte","Progenitor","Lymphoid","Lymphoid","Myeloid",
           "Myeloid","Myeloid","Myeloid","Progenitor","Granulocyte",
           "Progenitor","Myeloid","Myeloid","Myeloid","Granulocyte",
           "Dendritic cell","Dendritic cell",
           "Stromal contaminant","Stromal contaminant","Stromal contaminant")
)

saveRDS(combined.singlets, file = file.path(OUT_DIR, "combined.singlets_annotated.rds"))


# --- 11. SAVE UMAP COLOURED BY CONDITION -------------------------------------
Idents(combined.singlets) <- "Condition"
pdf(file.path(OUT_DIR, "umap_condition.pdf"), height = 4, width = 5)
DimPlot(combined.singlets, reduction = "wnn.umap", label = FALSE,
        repel = TRUE, label.size = 2.5)
dev.off()


# --- 12. MARKER DETECTION (RNA and ADT) --------------------------------------
DefaultAssay(combined.singlets) <- "RNA"
markers_RNA <- FindAllMarkers(combined.singlets, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
markers_RNA %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_RNA

DefaultAssay(combined.singlets) <- "ADT"
markers_ADT <- FindAllMarkers(combined.singlets)
markers_ADT %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_ADT

write.table(top10_RNA, file = file.path(OUT_DIR, "RNA_top10.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(top10_ADT, file = file.path(OUT_DIR, "ADT_top10.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)


# --- 13. SUBSET MATURE CELL TYPES FOR DE ANALYSIS ----------------------------
Idents(combined.singlets) <- "annotation_generic"
combined.singlets <- subset(combined.singlets,
                            idents = "Stromal contaminant", invert = TRUE)
mature_cell <- subset(combined.singlets, idents = "Progenitor", invert = TRUE)
Idents(mature_cell) <- "annotation"
mature_cell <- subset(mature_cell, idents = c("Mix-lymphoid", "B cells"), invert = TRUE)


# --- 14. HELPER: PER-CELL-TYPE DE ACROSS STIMULATION CONDITIONS --------------
run_de_by_celltype <- function(seurat_obj, assay_name, celltypes, conditions,
                               reference_condition = "Unstim") {
  DefaultAssay(seurat_obj) <- assay_name
  Idents(seurat_obj) <- "annotation_generic"
  all_de <- list()

  for (ct in celltypes) {
    obj <- subset(seurat_obj, subset = annotation_generic == ct)
    Idents(obj) <- "Condition"
    comparisons <- setdiff(unique(obj$Condition), reference_condition)

    for (cond in comparisons) {
      n1 <- sum(Idents(obj) == cond)
      n2 <- sum(Idents(obj) == reference_condition)

      if (n1 > 10 && n2 > 10) {
        de <- FindMarkers(obj, ident.1 = cond, ident.2 = reference_condition,
                          logfc.threshold = 0, min.pct = 0.1)
        de$gene       <- rownames(de)
        de$celltype   <- ct
        de$comparison <- paste0(cond, "_vs_", reference_condition)
        all_de[[paste(ct, cond, sep = "_")]] <- de
      }
    }
  }
  dplyr::bind_rows(all_de)
}


# --- 15. STIMULATION DE — RNA ------------------------------------------------
celltypes  <- c("Dendritic cell", "Myeloid", "Granulocyte", "Lymphoid")
conditions <- c("Unstim", "LPS", "CpG", "R848", "PolyI:C", "diAbizi", "Pan-stim")

de_RNA <- run_de_by_celltype(mature_cell, "RNA", celltypes, conditions)

top20_RNA <- de_RNA %>%
  filter(!grepl("^MT-", gene), avg_log2FC > 1) %>%
  group_by(celltype, comparison) %>%
  slice_head(n = 10) %>%
  ungroup()
top20_RNA_unique <- top20_RNA[!duplicated(top20_RNA$gene), ]

# Order cells for dotplot
mature_cell$cond_celltype <- paste(mature_cell$Condition,
                                   mature_cell$annotation_generic, sep = "_")
levels.use <- unlist(lapply(celltypes, function(ct) paste(conditions, ct, sep = "_")))
mature_cell$cond_celltype <- factor(mature_cell$cond_celltype, levels = levels.use)

pdf(file.path(OUT_DIR, "RNA_stim_dotplot.pdf"), height = 7, width = 17)
DotPlot(mature_cell, features = top20_RNA_unique$gene,
        group.by = "cond_celltype", cols = c("grey90", "darkred")) +
  RotatedAxis()
dev.off()


# --- 16. STIMULATION DE — ADT ------------------------------------------------
de_ADT <- run_de_by_celltype(mature_cell, "ADT", celltypes, conditions)

top20_ADT <- de_ADT %>%
  filter(!grepl("^MT-", gene), avg_log2FC > 1) %>%
  group_by(celltype, comparison) %>%
  slice_head(n = 10) %>%
  ungroup()
top20_ADT_unique <- top20_ADT[!duplicated(top20_ADT$gene), ]

# Clean up ADT feature names for display (strip "prot:" prefix and barcode suffix)
top20_ADT_unique$label <- gsub("^prot:", "", top20_ADT_unique$gene)
top20_ADT_unique$label <- sub("\\.[^.]+\\.[^.]+$", "", top20_ADT_unique$label)

pdf(file.path(OUT_DIR, "ADT_stim_dotplot.pdf"), height = 7, width = 16)
DotPlot(mature_cell, features = top20_ADT_unique$gene,
        group.by = "cond_celltype", cols = c("grey90", "darkgreen")) +
  RotatedAxis() +
  scale_x_discrete(labels = setNames(top20_ADT_unique$label,
                                     top20_ADT_unique$gene))
dev.off()

message("Script 03 complete. Output saved to: ", OUT_DIR)
