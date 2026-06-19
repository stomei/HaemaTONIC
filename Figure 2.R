library(ggpubr)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridisLite)
library(viridis)
library(paletteer)
library(scales)
library(ggbreak)
library(rstatix)
library(uwot)
library(Rphenograph)
library(scales)
library(patchwork)
library(readr)
library(ggrepel)


THEME=theme(text = element_text(size = 15, colour = "black"), 
            plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
            axis.text.y = element_text(colour = "black"), 
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "right"
)

Sample_ID <- read_excel("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/Figure2/Sample ID.xlsx")
concat_1 <- read_csv("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/Figure2/concat_1.csv")

cells <- left_join(concat_1, Sample_ID, by= "SampleID")

new_names <- ifelse(
  grepl("::", names(cells)),
  sub("^.*::\\s*([^ ]+).*$", "\\1", names(cells)),
  names(cells)
)

names(cells) <- new_names


# Remove unwanted columns
exclude <- c(
  "FSC-A", "FSC-H",
  "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-H",
  "IgM", "IgD", "CD200", "Comp-BFP-A", "Comp-GFP-A", "Comp-mCherry-A",
  "Live", "Time",
  "SampleID", "Donor", "Day"
)

markers <- cells %>%
  select(-all_of(exclude))


# UMAP
set.seed(123)

umap <- uwot::umap(markers, n_neighbors = 30, min_dist = 0.3, metric = "cosine", n_threads = 8)

# Phenograph clustering
Rphenograph_out <- Rphenograph(
  markers,
  k = 30
)

clusters <- membership(Rphenograph_out[[2]])

cells$UMAP1 <- umap[,1]
cells$UMAP2 <- umap[,2]
cells$cluster <- factor(clusters)


library(ggplot2)
library(dplyr)
library(ggrepel)

# Calculate cluster centers
cluster_centers <- cells %>%
  group_by(cluster) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

pdf("umap_facs.pdf", height = 4, width = 5)
ggplot(cells, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.2) +
  geom_text_repel(
    data = cluster_centers,
    aes(x = UMAP1, y = UMAP2, label = cluster),
    color = "black",
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()

ggplot(cells,
       aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.2) +
  theme_classic()

marker_names <- names(markers)
plots <- lapply(marker_names, function(marker){
  
  upper <- quantile(cells[[marker]], 0.99, na.rm = TRUE)
  
  ggplot(
    cells,
    aes(
      UMAP1,
      UMAP2,
      colour = pmin(.data[[marker]], upper)
    )
  ) +
    geom_point(size = 0.1) +
    scale_colour_viridis_c() +
    ggtitle(marker) +
    theme_classic()
})

pdf("umap_facs_marker.pdf", height = 10, width = 15)
wrap_plots(plots, ncol = 5)
dev.off()

#Seurat
library(Seurat)

# Seurat expects features (markers) as rows and cells as columns
mat <- t(as.matrix(markers))

seu <- CreateSeuratObject(
  counts = mat,
  assay = "Protein",
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu) <- "Protein"

# Store the scaled values in the data slot
seu <- SetAssayData(
  seu,
  assay = "Protein",
  layer = "data",
  new.data = mat
)


seu <- ScaleData(
  seu,
  features = rownames(seu),
  verbose = FALSE
)


seu <- RunPCA(seu, features = rownames(seu))
seu <- RunUMAP(seu, dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu)

pdf("umap_seu.pdf", height = 4, width = 5)
DimPlot(seu, label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
dev.off()

pdf("umap_facs_marker_seu.pdf", height = 13, width = 18)
FeaturePlot(
  seu,
  features = rownames(seu),
  reduction = "umap",
  ncol = 4,
  cols = c("lightgrey","darkgreen")
)
dev.off()

markers_de <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)
