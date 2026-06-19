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

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/wehisan/general/user_managed/grpu_naik.s_2/2021_Sequencing_Runs/ST158/Analysis/ST158/outs/filtered_feature_bc_matrix")
mgi.data <- Read10X(data.dir = "/stornext/General/data/user_managed/grpu_naik.s_2/2022_Sequencing_Runs/MGI_test_3/Test2_Test3_analsysis/10X_human_analysis/MGI_human/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the filtered data.
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "10X_analysis", min.cells = 3, min.features = 100)
pbmc.mgi <- CreateSeuratObject(counts = mgi.data, project = "10X_analysis", min.cells = 3, min.features = 100)

#merge


#Laod HTO and ADT re-seqeunced count
hto <- read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/HTO/ST158_HTO_counts_reseq.txt")
adt <- read.table(file="/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/HTO/ST158_ADT_counts_reseq.txt", header=FALSE)

rownames(hto) <- hto[,1]
hto <- hto[,-1]
rownames(adt) <- adt[,1]
adt <- adt[,-1]
colnames(adt) <- adt[1,]

#Select odd columns in adt data
col_odd <- seq_len(ncol(adt)) %% 2              # Create column indicator
adt_names <- adt[ , col_odd == 1]               # Subset odd columns
adt.count <- adt[ , col_odd == 0]               # Subset even columns
x <- colnames(adt_names)                        # Rename adt columns with Ab barcode
colnames(adt.count) <- x

TotalSeqA_reference_PDI <- read_excel("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/HTO/TotalSeqA_reference_PDI.xls")
reference <- TotalSeqA_reference_PDI[,c(1,5)]
reference <- reference %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed("anti-human")))
reference <- reference %>% 
  mutate(across(where(is.character), str_trim))
reference <- reference %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
reference <- reference %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed("/")))
adt.count <- t(adt.count)
adt.count <- as.data.frame(adt.count)
adt.count$sequence <- rownames(adt.count)
adt <- left_join(adt.count, reference, by="sequence")
rownames(adt) <- adt[,32536]
adt <- adt[,-32536]
adt = subset(adt, select = -sequence)
colnames(adt) <- paste(colnames(adt), "1", sep="-")


#Merge old ADT sequencing with new
data_UMAP_WNN <- readRDS("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/data_UMAP_WNN.rds")
pbmc.adt <- data_UMAP_WNN@assays$ADT@counts
adt.old <-  CreateSeuratObject(counts = pbmc.adt)
adt.new <-  CreateSeuratObject(counts = adt)
#adt.combined <- merge(adt.old, y = adt.new)
adt.old <- NormalizeData(adt.old)
adt.new <- NormalizeData(adt.new)
adt.old.matrix <- GetAssayData(object = adt.old, slot = "data")
adt.old.matrix <- as.data.frame(adt.old.matrix)
adt.new.matrix <- GetAssayData(object = adt.new, slot = "data")
adt.new.matrix <- as.data.frame(adt.new.matrix)
adt.old.matrix$id <- rownames(adt.old.matrix)
adt.old.matrix <- adt.old.matrix %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed("anti-human")))
adt.old.matrix <- adt.old.matrix %>% 
  mutate(across(where(is.character), str_trim))
rownames(adt.old.matrix) <-  adt.old.matrix[,37394]
adt.old.matrix <- adt.old.matrix[,-37394]
x <- colnames(adt.new.matrix)
adt.old.matrix = adt.old.matrix[,!(names(adt.old.matrix) %in% x)]
adt.old.matrix$id <- rownames(adt.old.matrix)
adt.new.matrix$id <- rownames(adt.new.matrix)
adt <- left_join(adt.new.matrix, adt.old.matrix, by= "id")
adt <- adt %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
adt <- adt %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed("/")))
rownames(adt) <- adt$id
adt = subset(adt, select = -c(id) )

#hto processing
colnames(hto)[2] <- "J-frozen"
colnames(hto)[4] <- "E"
colnames(hto)[6] <- "F-frozen"
colnames(hto)[8] <- "E-frozen"
colnames(hto)[10] <- "J"
colnames(hto)[12] <- "PBMC"
colnames(hto)[14] <- "F"
colnames(hto)[16] <- "H-frozen"
colnames(hto)[18] <- "J-frozen-2"
colnames(hto)[20] <- "H"
hto <- hto[,-c(1,3,5,7,9,11,13,15,17,19)]

hto <- t(hto)
hto <- as.data.frame(hto)
colnames(hto) <- paste(colnames(hto), "1", sep="-")
#colnames(hto1) <- glue::glue_collapse(x = colnames(hto1), "-1")
#joint.bcs <- intersect(colnames(pbmc.data@assays$RNA), colnames(hto))

# Subset RNA and HTO counts by joint cell barcodes
joint.bcs <- intersect(colnames(pbmc.mgi), colnames(hto))
joint.bcs <- intersect(joint.bcs, colnames(adt))
pbmc.mgi <- pbmc.mgi[, joint.bcs]
hto <- as.matrix(hto[, joint.bcs])
adt <- as.matrix(adt[, joint.bcs])
rownames(hto)

# Normalize RNA data with log normalization
pbmc.mgi <- NormalizeData(pbmc.mgi)
# Find and scale variable features
pbmc.mgi <- FindVariableFeatures(pbmc.mgi, selection.method = "mean.var.plot")
pbmc.mgi <- ScaleData(pbmc.mgi, features = VariableFeatures(pbmc.mgi))

#Add HTO to seurat object
pbmc.mgi[["HTO"]] <- CreateAssayObject(counts = hto)
pbmc.mgi[["ADT"]] <- CreateAssayObject(counts = adt)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.mgi <- NormalizeData(pbmc.mgi, assay = "HTO", normalization.method = "CLR")

#Demultiplex
pbmc.mgi <- HTODemux(pbmc.mgi, assay = "HTO", positive.quantile = 0.99)

#Visualise results of demultiplexing
# Global classification results
table(pbmc.mgi$HTO_classification.global)
Idents(pbmc.mgi) <- "HTO_classification.global"
VlnPlot(pbmc.mgi, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(pbmc.mgi, assay = "HTO")

# Extract the singlets
pbmc.singlet <- subset(pbmc.mgi, idents = "Singlet")
Idents(pbmc.singlet) <- "HTO_classification"
pbmc.singlet <- subset(pbmc.singlet, idents = "PBMC", invert = TRUE)
HTOHeatmap(pbmc.singlet, assay = "HTO")

#Filter low quality cells
pbmc.singlet[["percent.mt"]] <- PercentageFeatureSet(pbmc.singlet, pattern = "^MT-")
VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#pbmc.singlet <- subset(pbmc.singlet, subset = nFeature_RNA > 500 & nCount_RNA > 500 & percent.mt < 7) # Not done here

#Add demultiplex with SNP 
donors <- read.table("/wehisan/general/user_managed/grpu_naik.s_2/2022_Sequencing_Runs/MGI_test_3/Test2_Test3_analsysis/10X_human_analysis/GEX_demux/donor_ids.tsv", header = TRUE, sep ="\t")
donors <- select (donors, cell, donor_id)
pbmc.singlet@meta.data$cell=rownames(pbmc.singlet@meta.data)
pbmc.singlet@meta.data <- left_join(pbmc.singlet@meta.data, donors, by="cell")
rownames(pbmc.singlet@meta.data) <- pbmc.singlet@meta.data$cell

Idents(pbmc.singlet) <- "donor_id"
pbmc.singlet <- subset(pbmc.singlet, idents = "doublet", invert = TRUE)
pbmc.singlet <- subset(pbmc.singlet, idents = "donor0", invert = TRUE)

Idents(pbmc.singlet) <- "HTO_classification"
pbmc.singlet_E <- subset(pbmc.singlet, idents = c("E-frozen", "E"))
pbmc.singlet_H <- subset(pbmc.singlet, idents = c("H-frozen", "H"))
pbmc.singlet_F <- subset(pbmc.singlet, idents = c("F-frozen", "F"))
pbmc.singlet_J <- subset(pbmc.singlet, idents = c("J-frozen", "J", "J-frozen-2"))

Idents(pbmc.singlet_E) <- "donor_id"
p1 <- DimPlot(pbmc.singlet_E, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 5, cols= c(donor1= "#F8766D", donor2= "#0CB702", donor3="#00A9FF", donor4= "#E68613"))
Idents(pbmc.singlet_H) <- "donor_id"
p2 <- DimPlot(pbmc.singlet_H, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 5, cols= c(donor1= "#F8766D", donor2= "#0CB702", donor3="#00A9FF", donor4= "#E68613"))
Idents(pbmc.singlet_F) <- "donor_id"
p3 <- DimPlot(pbmc.singlet_F, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 5, cols= c(donor1= "#F8766D", donor2= "#0CB702", donor3="#00A9FF", donor4= "#E68613"))
Idents(pbmc.singlet_J) <- "donor_id"
p4 <- DimPlot(pbmc.singlet_J, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 5, cols= c(donor1= "#F8766D", donor2= "#0CB702", donor3="#00A9FF", donor4= "#E68613"))

pdf("donor_id_umap.pdf", height = 6, width = 9)
p1+p2+p3+p4
dev.off()

demultipexing_1 <- data.frame(table(pbmc.singlet_E@meta.data$donor_id), donor="E")
demultipexing_2 <- data.frame(table(pbmc.singlet_H@meta.data$donor_id), donor="H")
demultipexing_3 <- data.frame(table(pbmc.singlet_F@meta.data$donor_id), donor="F")
demultipexing_4 <- data.frame(table(pbmc.singlet_J@meta.data$donor_id), donor="J")
demultiplexing <- rbind(demultipexing_1, demultipexing_2, demultipexing_3, demultipexing_4)

demultiplexing <- demultiplexing%>%
  group_by(donor)%>%
  mutate(percentage=Freq/sum(Freq)*100)

THEME2=theme(text = element_text(size = 7, colour = "black"), 
             plot.title = element_text(size = 8, face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
             axis.text.y = element_text(colour = "black"), 
             axis.line = element_line(colour = "black"),
             panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             legend.position = "none"
)

donor_cols= c(donor1= "#F8766D", donor2= "#0CB702", donor3="#00A9FF", donor4= "#E68613")
demultiplexing$donor = factor(demultiplexing$donor, levels =c("H","J","E","F"))
pdf("donor_barplot.pdf", height = 4, width = 6)
ggplot(demultiplexing, aes(x=donor, y=percentage))+
  geom_col(aes(fill = Var1), col="black", width = 0.7)+
  scale_fill_manual(values = donor_cols)+
  scale_y_discrete(expand = c(0,0), limits = c(0,10,20,30,40,50,60,70,80,90,100))+
  labs(title="%Barcode",x="",y="Percentage")+
  THEME2
dev.off()

#Normalise data and make umap
#Now compute WNN analysis
DefaultAssay(pbmc.singlet) <- 'RNA'
pbmc.singlet <- NormalizeData(pbmc.singlet) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(pbmc.singlet) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(pbmc.singlet) <- rownames(pbmc.singlet[["ADT"]])
pbmc.singlet <- ScaleData(pbmc.singlet) %>% RunPCA(reduction.name = 'apca')

DefaultAssay(pbmc.singlet) <- 'HTO'
VariableFeatures(pbmc.singlet) <- rownames(pbmc.singlet[["HTO"]])
pbmc.singlet <- ScaleData(pbmc.singlet) %>% RunPCA(reduction.name = 'hpca')

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using pbmc.singlet[['weighted.nn']]
# The WNN graph can be accessed at pbmc.singlet[["wknn"]], 
# and the SNN graph used for clustering at pbmc.singlet[["wsnn"]]
# Cell-specific modality weights can be accessed at pbmc.singlet$RNA.weight
pbmc.singlet <- FindMultiModalNeighbors(
  pbmc.singlet, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
DefaultAssay(pbmc.singlet) <- 'RNA'
pbmc.singlet <- RunUMAP(pbmc.singlet, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#Clustering
#Find out how many PCA to use
DefaultAssay(pbmc.singlet) <- 'RNA'
ElbowPlot(pbmc.singlet)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc.singlet <- FindClusters(pbmc.singlet.2, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = 'pca', dims = 1:20, assay = 'RNA', 
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = 'hpca', dims = 1:18, assay = 'HTO', 
                        reduction.name = 'hto.umap', reduction.key = 'htoUMAP_')

p1 <- DimPlot(pbmc.singlet, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
p2 <- DimPlot(pbmc.singlet, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
p3 <- DimPlot(pbmc.singlet, reduction = 'adt.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()

pdf("umaps.pdf", height= 7, width = 15)
p2+p3+p1
dev.off()

p1

VlnPlot(pbmc.singlet, features = "RNA.weight", sort = TRUE, pt.size = 0.1) +
  NoLegend()

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator

# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(pbmc.singlet[["RNA"]])
Key(pbmc.singlet[["ADT"]])

DefaultAssay(pbmc.singlet) <- 'ADT'
FeatureScatter(pbmc.singlet, feature1 = "CD1c", feature2 = "CD370")

#We can identify cluster as expressing a aprticular surface marker
VlnPlot(pbmc.singlet, "CD14")

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
new.cluster.ids <- c("cDC2", "Eosinophils", "cDC2", "cDC1", "Mast cells", "Mast cells", "cDC1", "CD34+ progenitors", "cDC1", "DC3", "Neutrophils", "cDC1", "Eosinophils", "pDC", "Neutrophils", "activated DCs", "T cells", "Megakaryokytes", "cDC2", "Neutrophils", "cDC2", "cDC2","cDC2", "cDC1", "Mast cells", "cDC2")
names(new.cluster.ids) <- levels(pbmc.singlet)
pbmc.singlet <- RenameIdents(pbmc.singlet, new.cluster.ids)


#Heatmap
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
DefaultAssay(pbmc.singlet) <- 'RNA'
rna_markers <- FindAllMarkers(pbmc.singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
DefaultAssay(pbmc.singlet) <- 'ADT'
adt_markers <- FindAllMarkers(pbmc.singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "ADT")
#Save files
write.table(rna_markers, file= "wnn_rna_markers_pop.txt", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(adt_markers, file= "wnn_adt_markers_pop.txt", sep="\t", row.names=TRUE, col.names=TRUE)

jpeg("top_10_adt_markers_heatmap_pop.jpg", height= 7000, width = 6000, res = 500)
DefaultAssay(pbmc.singlet) <- 'ADT'
adt_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc.singlet, features = top10$gene) + NoLegend()
dev.off()

#write.table(top10, file= "wnn_top10_adt.txt", sep="\t", row.names=TRUE, col.names=TRUE)

jpeg("top_10_rna_markers_heatmap_pop.jpg", height= 8000, width = 6000, res = 500)
DefaultAssay(pbmc.singlet) <- 'RNA'
rna_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top5
DoHeatmap(pbmc.singlet, features = top5$gene) + NoLegend()
dev.off()
#write.table(top5, file= "wnn_top5_rna.txt", sep="\t", row.names=TRUE, col.names=TRUE)

DefaultAssay(pbmc.singlet) <- 'ADT'
cluster_9_markers <- FindMarkers(pbmc.singlet, ident.1 = 9 , ident.2 = NULL, only.pos = TRUE)

DefaultAssay(pbmc.singlet) <- 'RNA'
cluster_17_markers <- FindMarkers(pbmc.singlet, ident.1 = 17 , ident.2 = NULL, only.pos = TRUE)

DefaultAssay(pbmc.singlet) <- 'HTO'
p5 <- FeaturePlot(pbmc.singlet, features = c("F", "F-frozen"),
                  reduction = 'wnn.umap', max.cutoff = 2, ncol = 2)
p5
DefaultAssay(pbmc.singlet) <- 'ADT'
p6 <- FeaturePlot(pbmc.singlet, features = c("CD14","CD16", "CD1c", "CD370(CLEC9ADNGR1)", "CD123", "CD45RA", "CD56(NCAM)", "CD19", "CD7" ), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p7 <- FeaturePlot(pbmc.singlet, features = c("CD15(SSEA-1)", "CD11b", "CD163", "CD3", "TCRαβ" ,"CD34", "CD86", "CD40", "CD117(c-kit)"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p8 <- FeaturePlot(pbmc.singlet, features = c("CD19", "C5A", "CD89"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
CD64, XCR1 "CD117(c-kit)" 
pdf("new.2.pdf", height = 8, width = 10)
p7
dev.off()
pdf("new.pdf", height = 6, width = 10)
p7
dev.off()

#Add reference
#First load reference (downloaded from Seurat website): CITE-seq reference on 162,000 PBMC measured with 228 antibodies.
ref <- LoadH5Seurat("/stornext/Home/data/allstaff/t/tomei.s/10X Analysis/ST158/multi.h5seurat")

#Then normalize dataset with SCTransform since that is what was used for the reference
pbmc.singlet <- SCTransform(pbmc.singlet, verbose = FALSE)

#Find anchors between reference and query
anchors <- FindTransferAnchors(
  reference = ref,
  query = pbmc.singlet,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

#Transfer cell type labels from reference to query: l1 and l2 are just 2 different level of granularity
pbmc.singlet <- MapQuery(
  anchorset = anchors,
  query = pbmc.singlet,
  reference = ref,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(pbmc.singlet, reduction = "wnn.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc.singlet, reduction = "wnn.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 4 ,repel = TRUE) + NoLegend()
p1 + p2

VlnPlot(pbmc.singlet, features = "RNA.weight", group.by = 'predicted.celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend()

DoMultiBarHeatmap(object = m.bm.singlets.mgi, features = top5$gene, group.by="cluster", size = 2 ) + ggtitle("data2_mygenelist")
DoMultiBarHeatmap(object = pbmc.singlet, features = top10$gene, group.by="seurat_clusters", additional.group.by = "HTO_maxID", size = 2 )

saveRDS(pbmc.singlet, file = "/stornext/Home/data/allstaff/t/tomei.s/10X Analysis/ST158/ST158_R_files/pbmc.singlet.rds")
