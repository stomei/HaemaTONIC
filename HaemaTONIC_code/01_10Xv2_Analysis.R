# =============================================================================
# Script: 01_10Xv2_Analysis.R
# Description: 10x Chromium v2 multimodal (RNA + ADT + HTO) analysis pipeline
#              Includes HTO demultiplexing, SNP-based donor deconvolution,
#              WNN dimensionality reduction, clustering, and cell type annotation
#              via reference mapping (Azimuth/CITE-seq reference)
# HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications
# Author: Sara Tomei, WEHI
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(hdf5r)
library(stringr)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(SeuratData)
library(readxl)
library(DoMultiBarHeatmap)
library(rlang)
library(ggplot2)


# --- INPUT PATHS (update these to your local paths) --------------------------
# Raw 10x count matrices
PATH_10X_ST158  <- "path/to/ST158/outs/filtered_feature_bc_matrix"
PATH_10X_MGI    <- "path/to/MGI_human/outs/filtered_feature_bc_matrix"

# HTO and ADT count tables (re-sequenced)
PATH_HTO        <- "path/to/ST158_HTO_counts_reseq.txt"
PATH_ADT        <- "path/to/ST158_ADT_counts_reseq.txt"

# TotalSeq-A antibody reference (for renaming ADT features)
PATH_ADT_REF    <- "path/to/TotalSeqA_reference_PDI.xls"

# Previously processed WNN Seurat object (contains original ADT counts)
PATH_WNN_RDS    <- "path/to/data_UMAP_WNN.rds"

# SNP-based donor deconvolution output (e.g. from Demuxlet/Vireo)
PATH_DONOR_IDS  <- "path/to/donor_ids.tsv"

# Azimuth CITE-seq reference (downloaded from Seurat website)
# 162,000 PBMCs measured with 228 antibodies
PATH_AZIMUTH_REF <- "path/to/multi.h5seurat"

# Output directory
OUT_DIR <- "output/01_10Xv2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# --- 1. LOAD RAW COUNT DATA --------------------------------------------------
pbmc.data <- Read10X(data.dir = PATH_10X_ST158)
mgi.data  <- Read10X(data.dir = PATH_10X_MGI)

pbmc     <- CreateSeuratObject(counts = pbmc.data, project = "10X_analysis",
                               min.cells = 3, min.features = 100)
pbmc.mgi <- CreateSeuratObject(counts = mgi.data,  project = "10X_analysis",
                               min.cells = 3, min.features = 100)


# --- 2. LOAD AND PROCESS HTO AND ADT COUNTS ----------------------------------

# HTO counts
hto <- read.table(PATH_HTO)
rownames(hto) <- hto[, 1]
hto <- hto[, -1]

# ADT counts: alternating columns contain barcodes (odd) and counts (even)
adt <- read.table(PATH_ADT, header = FALSE)
rownames(adt) <- adt[, 1]
adt <- adt[, -1]
colnames(adt) <- adt[1, ]

col_odd        <- seq_len(ncol(adt)) %% 2
adt_names      <- adt[, col_odd == 1]   # barcode sequences
adt.count      <- adt[, col_odd == 0]   # UMI counts
colnames(adt.count) <- colnames(adt_names)

# Load TotalSeq-A reference and clean up antibody names
TotalSeqA_ref <- read_excel(PATH_ADT_REF)
reference <- TotalSeqA_ref[, c(1, 5)]
reference <- reference %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed("anti-human"))) %>%
  mutate(across(where(is.character), str_trim)) %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" "))) %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed("/")))

# Map ADT barcodes to antibody names
adt.count <- t(adt.count)
adt.count <- as.data.frame(adt.count)
adt.count$sequence <- rownames(adt.count)
adt_named <- left_join(adt.count, reference, by = "sequence")
rownames(adt_named) <- adt_named[, ncol(adt_named)]
adt_named <- adt_named[, !colnames(adt_named) %in% c("sequence", colnames(reference)[1])]
colnames(adt_named) <- paste(colnames(adt_named), "1", sep = "-")


# --- 3. MERGE NEW AND OLD ADT SEQUENCING RUNS --------------------------------
data_UMAP_WNN <- readRDS(PATH_WNN_RDS)
pbmc.adt      <- data_UMAP_WNN@assays$ADT@counts

adt.old <- CreateSeuratObject(counts = pbmc.adt)
adt.new <- CreateSeuratObject(counts = adt_named)

adt.old <- NormalizeData(adt.old)
adt.new <- NormalizeData(adt.new)

adt.old.matrix <- as.data.frame(GetAssayData(adt.old, slot = "data"))
adt.new.matrix <- as.data.frame(GetAssayData(adt.new, slot = "data"))

# Clean antibody names in old matrix
adt.old.matrix$id <- rownames(adt.old.matrix)
adt.old.matrix <- adt.old.matrix %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed("anti-human"))) %>%
  mutate(across(where(is.character), str_trim))
rownames(adt.old.matrix) <- adt.old.matrix$id
adt.old.matrix <- adt.old.matrix[, -which(colnames(adt.old.matrix) == "id")]

# Remove features already present in new ADT matrix (avoid duplication)
adt.old.matrix <- adt.old.matrix[, !colnames(adt.old.matrix) %in% colnames(adt.new.matrix)]
adt.old.matrix$id <- rownames(adt.old.matrix)
adt.new.matrix$id <- rownames(adt.new.matrix)

adt <- left_join(adt.new.matrix, adt.old.matrix, by = "id") %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" "))) %>%
  mutate(across(where(is.character), str_remove_all, pattern = fixed("/")))
rownames(adt) <- adt$id
adt <- subset(adt, select = -id)


# --- 4. PROCESS HTO MATRIX ---------------------------------------------------
# Rename HTO columns to meaningful sample names
colnames(hto)[c(2,4,6,8,10,12,14,16,18,20)] <- c("J-frozen","E","F-frozen","E-frozen","J","PBMC","F","H-frozen","J-frozen-2","H")
hto <- hto[, -c(1,3,5,7,9,11,13,15,17,19)]  # remove barcode-only columns

hto <- t(hto)
hto <- as.data.frame(hto)
colnames(hto) <- paste(colnames(hto), "1", sep = "-")


# --- 5. SUBSET TO JOINT BARCODES AND BUILD MULTIMODAL SEURAT OBJECT ---------
joint.bcs <- Reduce(intersect, list(colnames(pbmc.mgi), colnames(hto), colnames(adt)))
pbmc.mgi  <- pbmc.mgi[, joint.bcs]
hto       <- as.matrix(hto[, joint.bcs])
adt       <- as.matrix(adt[, joint.bcs])

# Add HTO and ADT assays
pbmc.mgi[["HTO"]] <- CreateAssayObject(counts = hto)
pbmc.mgi[["ADT"]] <- CreateAssayObject(counts = adt)

# Normalise RNA
pbmc.mgi <- NormalizeData(pbmc.mgi)
pbmc.mgi <- FindVariableFeatures(pbmc.mgi, selection.method = "mean.var.plot")
pbmc.mgi <- ScaleData(pbmc.mgi, features = VariableFeatures(pbmc.mgi))

# Normalise HTO with centred log-ratio
pbmc.mgi <- NormalizeData(pbmc.mgi, assay = "HTO", normalization.method = "CLR")


# --- 6. HTO DEMULTIPLEXING ---------------------------------------------------
pbmc.mgi <- HTODemux(pbmc.mgi, assay = "HTO", positive.quantile = 0.99)

# QC: global classification
table(pbmc.mgi$HTO_classification.global)
Idents(pbmc.mgi) <- "HTO_classification.global"
VlnPlot(pbmc.mgi, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
HTOHeatmap(pbmc.mgi, assay = "HTO")

# Retain singlets only; exclude PBMC hashtag (used as a reference)
pbmc.singlet <- subset(pbmc.mgi, idents = "Singlet")
Idents(pbmc.singlet) <- "HTO_classification"
pbmc.singlet <- subset(pbmc.singlet, idents = "PBMC", invert = TRUE)
HTOHeatmap(pbmc.singlet, assay = "HTO")


# --- 7. QUALITY CONTROL ------------------------------------------------------
pbmc.singlet[["percent.mt"]] <- PercentageFeatureSet(pbmc.singlet, pattern = "^MT-")
VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# NOTE: QC filtering thresholds (nFeature_RNA > 500, percent.mt < 7) were assessed
# but not applied at this stage; filtering is performed downstream after SNP demultiplexing


# --- 8. SNP-BASED DONOR DECONVOLUTION ----------------------------------------
# Adds donor identity from external tool (e.g. Demuxlet or Vireo)
donors <- read.table(PATH_DONOR_IDS, header = TRUE, sep = "\t")
donors <- dplyr::select(donors, cell, donor_id)

pbmc.singlet@meta.data$cell <- rownames(pbmc.singlet@meta.data)
pbmc.singlet@meta.data <- left_join(pbmc.singlet@meta.data, donors, by = "cell")
rownames(pbmc.singlet@meta.data) <- pbmc.singlet@meta.data$cell

# Remove doublets and unassigned cells
Idents(pbmc.singlet) <- "donor_id"
pbmc.singlet <- subset(pbmc.singlet, idents = c("doublet", "donor0"), invert = TRUE)


# --- 9. DONOR COMPOSITION VISUALISATION --------------------------------------
Idents(pbmc.singlet) <- "HTO_classification"
pbmc.singlet_E <- subset(pbmc.singlet, idents = c("E-frozen", "E"))
pbmc.singlet_H <- subset(pbmc.singlet, idents = c("H-frozen", "H"))
pbmc.singlet_F <- subset(pbmc.singlet, idents = c("F-frozen", "F"))
pbmc.singlet_J <- subset(pbmc.singlet, idents = c("J-frozen", "J", "J-frozen-2"))

donor_cols <- c(donor1 = "#F8766D", donor2 = "#0CB702",
                donor3 = "#00A9FF", donor4 = "#E68613")

plot_donor_umap <- function(obj, label) {
  Idents(obj) <- "donor_id"
  DimPlot(obj, reduction = "wnn.umap", label = FALSE, repel = TRUE,
          label.size = 5, cols = donor_cols) + ggtitle(label)
}

p1 <- plot_donor_umap(pbmc.singlet_E, "E")
p2 <- plot_donor_umap(pbmc.singlet_H, "H")
p3 <- plot_donor_umap(pbmc.singlet_F, "F")
p4 <- plot_donor_umap(pbmc.singlet_J, "J")

pdf(file.path(OUT_DIR, "donor_id_umap.pdf"), height = 6, width = 9)
p1 + p2 + p3 + p4
dev.off()

# Donor composition bar chart
THEME2 <- theme(
  text              = element_text(size = 7, colour = "black"),
  plot.title        = element_text(size = 8, face = "bold"),
  axis.text.x       = element_text(angle = 45, hjust = 1, colour = "black"),
  axis.text.y       = element_text(colour = "black"),
  axis.line         = element_line(colour = "black"),
  panel.background  = element_blank(),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  strip.background  = element_blank(),
  legend.position   = "none"
)

make_demux_table <- function(obj, donor_label) {
  data.frame(table(obj@meta.data$donor_id), donor = donor_label)
}

demultiplexing <- rbind(
  make_demux_table(pbmc.singlet_E, "E"),
  make_demux_table(pbmc.singlet_H, "H"),
  make_demux_table(pbmc.singlet_F, "F"),
  make_demux_table(pbmc.singlet_J, "J")
) %>%
  group_by(donor) %>%
  mutate(percentage = Freq / sum(Freq) * 100)

demultiplexing$donor <- factor(demultiplexing$donor, levels = c("H", "J", "E", "F"))

pdf(file.path(OUT_DIR, "donor_barplot.pdf"), height = 4, width = 6)
ggplot(demultiplexing, aes(x = donor, y = percentage)) +
  geom_col(aes(fill = Var1), col = "black", width = 0.7) +
  scale_fill_manual(values = donor_cols) +
  scale_y_discrete(expand = c(0, 0),
                   limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title = "% Barcode", x = "", y = "Percentage") +
  THEME2
dev.off()


# --- 10. WNN DIMENSIONALITY REDUCTION ----------------------------------------

# RNA
DefaultAssay(pbmc.singlet) <- "RNA"
pbmc.singlet <- NormalizeData(pbmc.singlet) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# ADT
DefaultAssay(pbmc.singlet) <- "ADT"
VariableFeatures(pbmc.singlet) <- rownames(pbmc.singlet[["ADT"]])
pbmc.singlet <- ScaleData(pbmc.singlet) %>%
  RunPCA(reduction.name = "apca")

# HTO
DefaultAssay(pbmc.singlet) <- "HTO"
VariableFeatures(pbmc.singlet) <- rownames(pbmc.singlet[["HTO"]])
pbmc.singlet <- ScaleData(pbmc.singlet) %>%
  RunPCA(reduction.name = "hpca")

# Weighted nearest neighbours (RNA + ADT)
pbmc.singlet <- FindMultiModalNeighbors(
  pbmc.singlet,
  reduction.list      = list("pca", "apca"),
  dims.list           = list(1:30, 1:18),
  modality.weight.name = "RNA.weight"
)

# UMAPs
DefaultAssay(pbmc.singlet) <- "RNA"
pbmc.singlet <- RunUMAP(pbmc.singlet, nn.name = "weighted.nn",
                        reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = "pca",  dims = 1:20, assay = "RNA",
                        reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = "apca", dims = 1:18, assay = "ADT",
                        reduction.name = "adt.umap", reduction.key = "adtUMAP_")

pdf(file.path(OUT_DIR, "umaps.pdf"), height = 7, width = 15)
DimPlot(pbmc.singlet, reduction = "rna.umap", label = TRUE, repel = TRUE, label.size = 5) + NoLegend() +
DimPlot(pbmc.singlet, reduction = "adt.umap", label = TRUE, repel = TRUE, label.size = 5) + NoLegend() +
DimPlot(pbmc.singlet, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
dev.off()


# --- 11. CLUSTERING ----------------------------------------------------------
ElbowPlot(pbmc.singlet)
pbmc.singlet <- FindClusters(pbmc.singlet, graph.name = "wsnn",
                             algorithm = 3, resolution = 1, verbose = FALSE)


# --- 12. CELL TYPE ANNOTATION ------------------------------------------------
current.cluster.ids <- 0:25
new.cluster.ids <- c(
  "cDC2", "Eosinophils", "cDC2", "cDC1", "Mast cells", "Mast cells",
  "cDC1", "CD34+ progenitors", "cDC1", "DC3", "Neutrophils", "cDC1",
  "Eosinophils", "pDC", "Neutrophils", "activated DCs", "T cells",
  "Megakaryocytes", "cDC2", "Neutrophils", "cDC2", "cDC2",
  "cDC2", "cDC1", "Mast cells", "cDC2"
)
names(new.cluster.ids) <- levels(pbmc.singlet)
pbmc.singlet <- RenameIdents(pbmc.singlet, new.cluster.ids)


# --- 13. DIFFERENTIAL EXPRESSION MARKERS ------------------------------------
DefaultAssay(pbmc.singlet) <- "RNA"
rna_markers <- FindAllMarkers(pbmc.singlet, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

DefaultAssay(pbmc.singlet) <- "ADT"
adt_markers <- FindAllMarkers(pbmc.singlet, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25, assay = "ADT")

write.table(rna_markers, file = file.path(OUT_DIR, "wnn_rna_markers.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adt_markers, file = file.path(OUT_DIR, "wnn_adt_markers.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)

# Heatmaps: top 10 markers per cluster
top10_adt <- adt_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_rna <- rna_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

jpeg(file.path(OUT_DIR, "top10_adt_markers_heatmap.jpg"), height = 7000, width = 6000, res = 500)
DefaultAssay(pbmc.singlet) <- "ADT"
DoHeatmap(pbmc.singlet, features = top10_adt$gene) + NoLegend()
dev.off()

jpeg(file.path(OUT_DIR, "top10_rna_markers_heatmap.jpg"), height = 8000, width = 6000, res = 500)
DefaultAssay(pbmc.singlet) <- "RNA"
DoHeatmap(pbmc.singlet, features = top10_rna$gene) + NoLegend()
dev.off()


# --- 14. REFERENCE MAPPING WITH AZIMUTH -------------------------------------
# NOTE: Reference must be downloaded from: https://zenodo.org/record/6328210
# (Hao et al. 2021 CITE-seq PBMC reference, 162,000 cells, 228 antibodies)
ref <- LoadH5Seurat(PATH_AZIMUTH_REF)

pbmc.singlet <- SCTransform(pbmc.singlet, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference           = ref,
  query               = pbmc.singlet,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims                = 1:50
)

pbmc.singlet <- MapQuery(
  anchorset         = anchors,
  query             = pbmc.singlet,
  reference         = ref,
  refdata           = list(
    celltype.l1   = "celltype.l1",
    celltype.l2   = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model   = "wnn.umap"
)

p1 <- DimPlot(pbmc.singlet, reduction = "wnn.umap", group.by = "predicted.celltype.l1",
              label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
p2 <- DimPlot(pbmc.singlet, reduction = "wnn.umap", group.by = "predicted.celltype.l2",
              label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
p1 + p2


# --- 15. SAVE FINAL OBJECT ---------------------------------------------------
saveRDS(pbmc.singlet, file = file.path(OUT_DIR, "pbmc.singlet_annotated.rds"))

message("Script 01 complete. Output saved to: ", OUT_DIR)
