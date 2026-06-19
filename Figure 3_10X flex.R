#Use singleR and Immgen
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")
BiocManager::install("GenomeInfoDb")
BiocManager::install("BSgenome")
BiocManager::install("ensembldb")
install.packages("remotes")
devtools::install_github("satijalab/AzimuthAPI")
library(Seurat)
library(SeuratObject)
library(dplyr)
library(hdf5r)
library(celldex)
library(SingleR)
library(SeuratObject)
library(SingleCellExperiment)
library(AzimuthAPI)
library(ggplot2)


load_flex_sample <- function(path, sample_name) {
  data <- Read10X_h5(path)
  
  rna <- data[["Gene Expression"]]
  antibody <- data[["Antibody Capture"]]
  
  # Remove HTOs
  adt <- antibody[!grepl("hash", rownames(antibody)), ]
  hto <- antibody[grepl("hash", rownames(antibody)), ]
  
  obj <- CreateSeuratObject(counts = rna, project = sample_name)
  obj[["ADT"]] <- CreateAssayObject(counts = adt)
  obj[["HTO"]] <- CreateAssayObject(counts = hto)
  
  obj <- NormalizeData(obj)
  obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR")
  
  obj$sample <- sample_name
  
  return(obj)
}

s1 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool1/sample_filtered_feature_bc_matrix.h5", "Pool1")
s2 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool2/sample_filtered_feature_bc_matrix.h5", "Pool2")
s3 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool3/sample_filtered_feature_bc_matrix.h5", "Pool3")
s4 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool4/sample_filtered_feature_bc_matrix.h5", "Pool4")
s5 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool5/sample_filtered_feature_bc_matrix.h5", "Pool5")
s6 <- load_flex_sample("/vast/projects/ST227/multipro/multipro/Analysis/multipro_hash_new/outs/per_sample_outs/Pool6/sample_filtered_feature_bc_matrix.h5", "Pool6")

combined <- merge(s1, y = list(s2, s3, s4, s5, s6),
                  add.cell.ids = c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"))

combined <- JoinLayers(combined)

# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
combined <- subset(combined, subset = nFeature_RNA > 200)

#Demultiplexing using HTOs
combined <- NormalizeData(combined, assay = "HTO", normalization.method = "CLR")
combined <- HTODemux(combined, assay = "HTO", positive.quantile = 0.99)

table(combined$HTO_classification.global)
HTOHeatmap(combined, assay = "HTO", ncells = 5000)
Idents(combined) <- "HTO_maxID"
RidgePlot(combined, assay = "HTO", features = rownames(combined[["HTO"]])[1:4], ncol = 2)

#Add Donor and conditions to the metadata based on the Pool and the hashtag used
metadata <- combined@meta.data
write.table(metadata, file= "metadata.txt", sep="\t", row.names=TRUE, col.names=TRUE)

metadata_resolved=read.table("metadata_resolved.txt",header = TRUE, row.names = 1, sep = "\t")
combined@meta.data <- metadata_resolved

Idents(combined) <- "HTO_classification.global"
combined <- subset(combined, idents = "Negative", invert = TRUE)
combined.singlets <- subset(combined, idents = "Singlet")
metadata <- combined.singlets@meta.data
metadata <- separate(metadata, resolved_sample, c("Condition", "Donor"), sep=" ")
combined.singlets@meta.data <- metadata

Idents(combined) <- "resolved_sample"
pbmc <- subset(combined, idents = "PBMC")
combined <- subset(combined, idents = "PBMC", invert = TRUE)

combined.singlets[["RNA"]] <- JoinLayers(combined.singlets[["RNA"]])

#Wnn
DefaultAssay(combined.singlets) <- 'RNA'
combined.singlets <- NormalizeData(combined.singlets) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(combined.singlets) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(combined.singlets) <- rownames(combined.singlets[["ADT"]])
combined.singlets <- NormalizeData(combined.singlets, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#Check howm many PC to use
ElbowPlot(combined.singlets, reduction = "pca")

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
combined.singlets <- FindMultiModalNeighbors(
  combined.singlets, reduction.list = list("pca", "apca"), 
  dims.list = list(1:20, 1:16), modality.weight.name = "RNA.weight"
)

combined.singlets <- RunUMAP(combined.singlets, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

combined.singlets <- RunUMAP(combined.singlets, reduction = 'pca', dims = 1:20, assay = 'RNA', 
                             reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
combined.singlets <- RunUMAP(combined.singlets, reduction = 'apca', dims = 1:16, assay = 'ADT', 
                             reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

#Idents(combined.singlets) <- "resolved_sample"
Idents(combined.singlets) <- "seurat_clusters"
Idents(combined.singlets) <- "annotation"
Idents(combined.singlets) <- "Donor"
Idents(combined.singlets) <- "Condition"
Idents(combined.singlets) <- "azimuth_medium"
p1 <- DimPlot(combined.singlets, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 2.5)

p3 <- DimPlot(combined.singlets, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(combined.singlets, reduction = 'adt.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

pdf("umap_condition.pdf", height = 4, width = 5)
p1
dev.off()

DefaultAssay(combined.singlets) <- 'RNA'
markers_RNA <- FindAllMarkers(
  combined.singlets,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

DefaultAssay(combined.singlets) <- 'ADT'
markers_ADT <- FindAllMarkers(combined.singlets)

markers_RNA %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

markers_ADT %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> ADT_top10

write.table(top10, file= "RNA_top_10.txt", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(ADT_top10, file= "ADT_top_10.txt", sep="\t", row.names=TRUE, col.names=TRUE)

#Compare with BM dataset
sce.fig <- as.SingleCellExperiment(combined.singlets, assay = "RNA")
ref <- fetchReference("novershtern_hematopoietic", "2024-02-26")
pred.main <- SingleR(test = sce.fig, ref = ref, labels = ref$label.main)

combined.singlets[["SingleR.labels"]] <- pred.main$labels
Idents(combined.singlets) <- "SingleR.labels"

plotScoreHeatmap(pred.main,
                 clusters = combined.singlets$seurat_clusters, order.by.clusters = TRUE,
                 show.labels = TRUE, show.pruned = FALSE,
                 filename="heatmap_SingleR_main.png", width = 7, height = 8
)

# Annotate a Seurat object
combined.singlets <- CloudAzimuth(combined.singlets)
combined.singlets_qc <- RunUMAP(combined.singlets, dims = 1:128, reduction = "azimuth_embed", reduction.name = "azimuth_umap")
p3 <- DimPlot(combined.singlets_qc, group.by = "final_level_labels", label.size = 1.5, label = T, reduction = "azimuth_umap") + NoLegend()
p4 <- FeaturePlot(combined.singlets_qc, features = "final_level_confidence", reduction = "azimuth_umap")
combined.singlets_qc <- PrepLabel(combined.singlets_qc, "azimuth_fine", "azimuth_fine_filtered", cutoff = 10)
p7 <- DimPlot(combined.singlets_qc, group.by = "azimuth_fine_filtered", label.size = 3, label = T, reduction = "azimuth_umap") +
  NoLegend()
p7

plots <- make_azimuth_QC_heatmaps(combined.singlets_qc)
DefaultAssay(combined.singlets) <- 'RNA'
#HSPC
FeaturePlot(combined.singlets, c("CD34","SPINK2","GATA2","MEIS1","HLF"), reduction = 'wnn.umap')

#Myeloid
FeaturePlot(combined.singlets, c("MPO","AZU1","ELANE","PRTN3","CTSG"),reduction = 'wnn.umap')

#Dendritic cells
FeaturePlot(combined.singlets, c("FCER1A","CLEC10A", "CD1C","CLEC9A", "IRF8","BATF3"),reduction = 'wnn.umap')

#B cells
FeaturePlot(combined.singlets, c("MS4A1","CD79A", "CD79B","CD74"), reduction = 'wnn.umap')

#T cells
FeaturePlot(combined.singlets, c("CD3D","CD3E","TRBC1", "TRBC2","IL7R"), reduction = 'wnn.umap')

#NK Cells
FeaturePlot(combined.singlets, c("NKG7","GNLY","PRF1","KLRD1"), reduction = 'wnn.umap')

#Ery
FeaturePlot(combined.singlets, c("HBB","HBA1","HBA2","ALAS2"), reduction = 'wnn.umap')

#MK
FeaturePlot(combined.singlets, c("PPBP","PF4","ITGA2B","GP9"), reduction = 'wnn.umap')

#Monocytes
FeaturePlot(combined.singlets, c("LYZ","S100A8","S100A9","CTSS","FCN1","SAT1"), reduction = 'wnn.umap')

#Eosinophils
FeaturePlot(combined.singlets, c("KIT","CPA3","TPSAB1","TPSB2"), reduction = 'wnn.umap')

#Basophils
FeaturePlot(combined.singlets, c("IL3RA","GATA1","CLC","HDC"), reduction = 'wnn.umap')

#pDC
FeaturePlot(combined.singlets, features = c("IL3RA","CLEC4C", "GZMB", "JCHAIN", "PLD4", "TCF4"), reduction = 'wnn.umap')

#Misc
FeaturePlot(combined.singlets,c("KIT","TPSAB1","CPA3","MS4A2","FCER1A","CLC","PPBP","PF4","ITGA2B"),reduction = 'wnn.umap')

p5 <- FeaturePlot(combined.singlets, features = c("prot:CD14.60253.1","prot:CD14.65056.1","prot:CD15.65298.1", "prot:CD163.16646.1", "prot:CD163.65169.1", "prot:CD16.65090.1"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(combined.singlets, features = c("prot:CD7.65203.1","prot:CD56.65264.1", "prot:CD41a.65173.1", "prot:CD63.65255.1", "prot:CD42b.65163.1", "prot:CD36.65248.1", "prot:CD9.20597.1", "prot:CD45RA.65226.1"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 2)
p6

#Annotate
combined.singlets$annotation <- plyr::mapvalues(combined.singlets$seurat_clusters, from = c("1", "22", "18", "34", "36", "25", "29", "6", "24", "8", "9", "28", "0", "11", "10", "33", "35","20", "4", "31", "3", "15", "13", "7", "32", "26", "17", "2", "16", "19", "27", "12", "14", "21", "5", "23", "30"), 
                                                to   = c("cDC1","cDC1","cDC1","cDC1","DC progenitor","cDC1","B cells","cDC2","DC progenitor","cDC2","mregDC", "NK cell", "Mast cell","Mast cell","Mast cell","Mast cell","Mast cell", "Basophils", "HSPC", "NK cell", "Mix-lymphoid","Monocytes","Monocytes","Macrophages","Macrophages","Granulocyte progenitor", "Eosinophils", "Baso/Eos progenitor", "Monocytes", "Macrophages","Macrophages" ,"Neutrophils", "pDC", "pDC", "Stromal contaminant","Stromal contaminant","Stromal contaminant")
)

combined.singlets$annotation_generic <- plyr::mapvalues(combined.singlets$seurat_clusters, from = c("1", "22", "18", "34", "36", "25", "29", "6", "24", "8", "9", "28", "0", "11", "10", "33", "35","20", "4", "31", "3", "15", "13", "7", "32", "26", "17", "2", "16", "19", "27", "12", "14", "21", "5", "23", "30"), 
                                                to   = c("Dendritic cell","Dendritic cell","Dendritic cell","Dendritic cell","Progenitor","Dendritic cell","Lymphoid","Dendritic cell","Progenitor","Dendritic cell","Dendritic cell", "Lymphoid", "Granulocyte","Granulocyte","Granulocyte","Granulocyte","Granulocyte", "Granulocyte", "Progenitor", "Lymphoid", "Lymphoid","Myeloid","Myeloid","Myeloid","Myeloid","Progenitor", "Granulocyte", "Progenitor", "Myeloid", "Myeloid","Myeloid" ,"Granulocyte", "Dendritic cell", "Dendritic cell", "Stromal contaminant","Stromal contaminant","Stromal contaminant")
)

markers <- c(
  # Monocytes
  "LYZ","S100A8","S100A9", "SAT1",
  
  # Macrophages
  "C1QC","APOE","MSR1",
  
  # Neutrophils
  "MPO","ELANE","FCGR3B",
  
  # Eosinophils
  "CLC","EPX","PRG2",
  
  # Basophils
  "HDC","IL3RA","MS4A2"
)

DotPlot(
  combined.singlets.1,
  features = markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

dc.markers <- c(
  
  # cDC1
  "CLEC9A","XCR1","CADM1",
  
  # cDC2
  "CD1C","CLEC10A","FCER1A",
  
  # pDC
  "GZMB","JCHAIN","PLD4",
  
  # mregDC
  "LAMP3","CCR7","CD83","IDO1"
)

DotPlot(
  combined.singlets,
  features = dc.markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

saveRDS(combined.singlets, file = "/vast/projects/ST227/multipro/multipro/R_script/combined.singlets.rds")

#Remove progenitor and stromal cells
Idents(combined.singlets) <- "annotation_generic"
combined.singlets<- subset(combined.singlets, idents = "Stromal contaminant", invert = TRUE)
mature_cell<- subset(combined.singlets, idents = "Progenitor", invert = TRUE)
Idents(mature_cell) <- "annotation"
mature_cell<- subset(mature_cell, idents = "Mix-lymphoid", invert = TRUE)
mature_cell<- subset(mature_cell, idents = "B cells", invert = TRUE)

# Set identities to condition
Idents(mature_cell) <- "annotation_generic"
DefaultAssay(mature_cell) <- 'RNA'

celltypes <- unique(mature_cell$annotation_generic)

all.de <- list()

for(ct in celltypes){
  
  obj <- subset(
    mature_cell,
    subset = annotation_generic == ct
  )
  
  Idents(obj) <- "Condition"
  
  comparisons <- setdiff(
    unique(obj$Condition),
    "Unstim"
  )
  
  for(cond in comparisons){
    
    if(sum(Idents(obj) == cond) > 10 &
       sum(Idents(obj) == "Unstim") > 10){
      
      de <- FindMarkers(
        obj,
        ident.1 = cond,
        ident.2 = "Unstim",
        logfc.threshold = 0,
        min.pct = 0.1
      )
      
      de$gene <- rownames(de)
      de$celltype <- ct
      de$comparison <- paste0(cond, "_vs_Unstim")
      
      all.de[[paste(ct, cond, sep = "_")]] <- de
    }
  }
}

de.all <- dplyr::bind_rows(all.de)

markers_RNA <- de.all %>%
  filter(!grepl("^MT-", gene)) %>%
  filter(avg_log2FC > 1) %>%
  group_by(celltype, comparison) %>%
  slice_head(n = 10) %>%
  ungroup()

top20.unique <- markers_RNA[!duplicated(markers_RNA$gene), ]

# Define cell type order
celltypes <- c("cDC1", "cDC2", "pDC", "mregDC", "Neutrophils","Mast cell", "Eosinophils", "Basophils", "Monocytes", "Macrophages", "NK cell")
celltypes <- c("Dendritic cell", "Myeloid", "Granulocyte", "Lymphoid")
# Define desired condition order
conditions <- c("Unstim","LPS","CpG","R848", "PolyI:C", "diAbizi","Pan-stim")

# Create combined metadata column
mature_cell$cond_celltype <- paste(mature_cell$Condition, mature_cell$annotation_generic,sep = "_")

# Generate factor levels:
# Cell type first, then condition
levels.use <- unlist(
  lapply(
    celltypes,
    function(ct) paste(conditions, ct, sep = "_")
  )
)

# Apply ordering
mature_cell$cond_celltype <- factor(mature_cell$cond_celltype, levels = levels.use)
pdf("RNA_stim_dotplot_gen.pdf", height = 7, width= 17)
DotPlot(
  mature_cell,
  features = top20.unique$gene,
  group.by = "cond_celltype",
  cols = c("grey90", "darkred")
) +
  RotatedAxis()
dev.off()

#Do the same for ADT
DefaultAssay(mature_cell) <- 'ADT'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")
celltypes <- unique(mature_cell$annotation_generic)

all.de <- list()

for(ct in celltypes){
  
  obj <- subset(
    mature_cell,
    subset = annotation_generic == ct
  )
  
  Idents(obj) <- "Condition"
  
  comparisons <- setdiff(
    unique(obj$Condition),
    "Unstim"
  )
  
  for(cond in comparisons){
    
    if(sum(Idents(obj) == cond) > 10 &
       sum(Idents(obj) == "Unstim") > 10){
      
      de <- FindMarkers(
        obj,
        ident.1 = cond,
        ident.2 = "Unstim",
        logfc.threshold = 0,
        min.pct = 0.1
      )
      
      de$gene <- rownames(de)
      de$celltype <- ct
      de$comparison <- paste0(cond, "_vs_Unstim")
      
      all.de[[paste(ct, cond, sep = "_")]] <- de
    }
  }
}

de.all <- dplyr::bind_rows(all.de)

# Reorder columns

markers_ADT <- de.all %>%
  filter(!grepl("^MT-", gene)) %>%
  group_by(celltype, comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top20

top20.unique <- top20[!duplicated(top20$gene), ]

top20.unique$label <- top20.unique$gene
top20.unique$label <- gsub("^prot:", "", top20.unique$label)
top20.unique$label <- sub("\\.[^.]+\\.[^.]+$", "", top20.unique$label)

# Define cell type order
celltypes <- c("cDC1", "cDC2", "pDC", "mregDC", "Neutrophils","Mast cell", "Eosinophils", "Basophils", "Monocytes", "Macrophages", "NK cell", "B cells")
celltypes <- c("Dendritic cell", "Myeloid", "Granulocyte", "Lymphoid")
# Define desired condition order
conditions <- c("Unstim","LPS","CpG","R848", "PolyI:C", "diAbizi","Pan-stim")

# Create combined metadata column
mature_cell$cond_celltype <- paste(mature_cell$Condition, mature_cell$annotation_generic,sep = "_")

# Generate factor levels:
# Cell type first, then condition
levels.use <- unlist(
  lapply(
    celltypes,
    function(ct) paste(conditions, ct, sep = "_")
  )
)

# Apply ordering
mature_cell$cond_celltype <- factor(mature_cell$cond_celltype, levels = levels.use)

p2 <- DotPlot(
  mature_cell,
  features = top20.unique$gene,
  group.by = "cond_celltype",
  cols = c("grey90", "darkgreen")
) + RotatedAxis() 


pdf("ADT_stim_dotplot_gen.pdf", height = 7, width= 16)
p2 +
  scale_x_discrete(
    labels = setNames(
      top20.unique$label,
      top20.unique$feature
    )
  )
dev.off()
#Do the same for ADT
DefaultAssay(Monocytes) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Monocytes,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top20
DoHeatmap(Monocytes, features = top20$gene) + NoLegend()

markers_ADT_mono <- FindAllMarkers(Monocytes)
markers_ADT_mono %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

DoHeatmap(Monocytes, features = top5$gene) + NoLegend()


# Set identities to condition
Idents(Monocytes) <- "Condition"
DefaultAssay(Monocytes) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Monocytes,
    ident.1 = cond,
    ident.2 = "Unstim",
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Monocytes$Condition <- factor(Monocytes$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
Monocytes <- Monocytes[, order(Monocytes$Condition)]
DoHeatmap(Monocytes, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(Monocytes) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Monocytes,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top20
DoHeatmap(Monocytes, features = top20$gene) + NoLegend()

markers_ADT_mono <- FindAllMarkers(Monocytes)
markers_ADT_mono %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

DoHeatmap(Monocytes, features = top5$gene) + NoLegend()

#DC
# Set identities to condition
Idents(DC) <- "Condition"
DefaultAssay(DC) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(DC,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

DC$Condition <- factor(DC$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
DC <- DC[, order(DC$Condition)]
DoHeatmap(DC, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(DC) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(DC,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(DC, features = top20$gene) + NoLegend()

markers_ADT_mono <- FindAllMarkers(cDC1)
markers_ADT_mono %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

DoHeatmap(cDC1, features = top5$gene) + NoLegend()

#T cells
# Set identities to condition
Idents(T_cell) <- "Condition"
DefaultAssay(T_cell) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

T_cell$Condition <- factor(T_cell$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
T_cell <- T_cell[, order(T_cell$Condition)]
DoHeatmap(T_cell, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(T_cell) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(T_cell, features = top20$gene) + NoLegend()

#Mast cells
# Set identities to condition
Idents(Mast) <- "Condition"
DefaultAssay(Mast) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Mast,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Mast$Condition <- factor(Mast$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
Mast <- Mast[, order(Mast$Condition)]
DoHeatmap(Mast, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(T_cell) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(T_cell, features = top20$gene) + NoLegend()

#Eosinophils
# Set identities to condition
Idents(Eosinophils) <- "Condition"
DefaultAssay(Eosinophils) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Eosinophils,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Eosinophils$Condition <- factor(Eosinophils$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
Eosinophils <- Eosinophils[, order(Eosinophils$Condition)]
DoHeatmap(Eosinophils, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(T_cell) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(T_cell, features = top20$gene) + NoLegend()

#Eosinophils
# Set identities to condition
Idents(Eosinophils) <- "Condition"
DefaultAssay(Eosinophils) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Eosinophils,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Eosinophils$Condition <- factor(Eosinophils$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
Eosinophils <- Eosinophils[, order(Eosinophils$Condition)]
DoHeatmap(Eosinophils, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(T_cell) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(T_cell, features = top20$gene) + NoLegend()


#Basophils
# Set identities to condition
Idents(Basophils) <- "Condition"
DefaultAssay(Basophils) <- 'RNA'

# Conditions to compare against Unstim
conditions <- c("Pan-stim", "R848", "PolyI:C", "CpG", "diAbizi", "LPS")

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(Basophils,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_RNA <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Basophils$Condition <- factor(Basophils$Condition,levels = c("Unstim","LPS","CpG","R848","PolyI:C","diAbizi","Pan-stim"))
Basophils <- Basophils[, order(Basophils$Condition)]
DoHeatmap(Basophils, features = top20$gene, group.by = "Condition") + NoLegend()
#Do the same for ADT
DefaultAssay(T_cell) <- 'ADT'

# Run DE for each comparison
de.list <- lapply(conditions, function(cond) {
  
  message("Running: ", cond, " vs Unstim")
  
  res <- FindMarkers(T_cell,
                     ident.1 = cond,
                     ident.2 = "Unstim",
                     logfc.threshold = 0,
                     min.pct = 0.1,
                     test.use = "wilcox"
  )
  
  res$gene <- rownames(res)
  res$comparison <- paste0(cond, "_vs_Unstim")
  
  res
})

# Combine into one table
de.all <- dplyr::bind_rows(de.list)

# Reorder columns
de.all <- de.all %>% dplyr::select(comparison, gene, dplyr::everything())
markers_ADT <- de.all %>%
  group_by(comparison) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(T_cell, features = top20$gene) + NoLegend()




