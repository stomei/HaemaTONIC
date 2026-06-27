# =============================================================================
# Script: 02_Figure_2.R
# Description: Flow cytometry data analysis for Figure 2
#              Loads concatenated flow cytometry data, computes UMAP and
#              Phenograph clustering, and generates marker feature plots
# Paper: HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications
# Author: Sara Tomei, WEHI
# =============================================================================

# --- Libraries ----------------------------------------------------------------
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
library(patchwork)
library(readr)
library(ggrepel)
library(Seurat)


# --- INPUT PATHS (update these to your local paths) --------------------------
PATH_SAMPLE_ID  <- "path/to/Figure2/Sample_ID.xlsx"
PATH_CONCAT_CSV <- "path/to/Figure2/Spectral_unmixed_values.csv"

# Output directory
OUT_DIR <- "output/02_Figure2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# --- THEME -------------------------------------------------------------------
THEME <- theme(
  text             = element_text(size = 15, colour = "black"),
  plot.title       = element_text(size = 20, face = "bold"),
  axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
  axis.text.y      = element_text(colour = "black"),
  axis.line        = element_line(colour = "black"),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  legend.position  = "right"
)


# --- 1. LOAD AND CLEAN DATA --------------------------------------------------
Sample_ID <- read_excel(PATH_SAMPLE_ID)
concat_1  <- read_csv(PATH_CONCAT_CSV)

cells <- left_join(concat_1, Sample_ID, by = "SampleID")

# Simplify column names: strip "Marker :: Description" format to just "Marker"
new_names <- ifelse(
  grepl("::", names(cells)),
  sub("^.*::\\s*([^ ]+).*$", "\\1", names(cells)),
  names(cells)
)
names(cells) <- new_names

# Remove non-marker columns (scatter, live/dead, sample metadata)
exclude <- c(
  "FSC-A", "FSC-H",
  "SSC-A", "SSC-B-A", "SSC-B-H", "SSC-H",
  "IgM", "IgD", "CD200",
  "Comp-BFP-A", "Comp-GFP-A", "Comp-mCherry-A",
  "Live", "Time",
  "SampleID", "Donor", "Day"
)
markers <- cells %>% select(-all_of(exclude))


# --- 2. UMAP -----------------------------------------------------------------
set.seed(123)
umap <- uwot::umap(markers, n_neighbors = 30, min_dist = 0.05,
                   metric = "cosine", n_threads = 8)


# --- 3. PHENOGRAPH CLUSTERING ------------------------------------------------
Rphenograph_out <- Rphenograph(markers, k = 70)
clusters        <- membership(Rphenograph_out[[2]])

cells$UMAP1    <- umap[, 1]
cells$UMAP2    <- umap[, 2]
cells$cluster  <- factor(clusters)


# --- 4. CLUSTER UMAP PLOT ---------------------------------------------------
cluster_centers <- cells %>%
  group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

pdf(file.path(OUT_DIR, "umap_facs.pdf"), height = 4, width = 5)
ggplot(cells, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.2) +
  geom_text_repel(
    data     = cluster_centers,
    aes(x = UMAP1, y = UMAP2, label = cluster),
    color    = "black",
    size     = 4,
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


# --- 5. MARKER FEATURE PLOTS -------------------------------------------------
marker_names <- names(markers)

plots <- lapply(marker_names, function(marker) {
  upper <- quantile(cells[[marker]], 0.99, na.rm = TRUE)
  ggplot(cells, aes(UMAP1, UMAP2, colour = pmin(.data[[marker]], upper))) +
    geom_point(size = 0.1) +
    scale_colour_viridis_c() +
    ggtitle(marker) +
    theme_classic()
})

pdf(file.path(OUT_DIR, "umap_facs_marker.pdf"), height = 10, width = 15)
wrap_plots(plots, ncol = 5)
dev.off()


