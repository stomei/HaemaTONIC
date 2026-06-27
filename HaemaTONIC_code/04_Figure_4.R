# =============================================================================
# Script: 04_Figure_4.R
# Description: Drug screen analysis for Figure 4
#              Loads flow cytometry counts/frequencies from drug treatment
#              experiments, normalises to DMSO control, and generates
#              heatmaps and line plots across drug concentrations and timepoints
# HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications
# Author: Sara Tomei, WEHI
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridisLite)
library(viridis)
library(paletteer)
library(scales)


# --- INPUT PATHS (update these to your local paths) --------------------------
PATH_DRUG_DATA <- "path/to/Drug_Screen/Drug_screen.xlsx"

# Output directory
OUT_DIR <- "output/04_Figure4"
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

# Colour scheme: each drug gets a gradient from dark (high dose) to light (low dose)
drug_col <- c(
  "DMSO_10uM"         = "grey60",  "DMSO_1nM"          = "grey90",
  "S4_10uM"           = "#0F6B99", "S4_1uM"            = "#2C85B2",
  "S4_100nM"          = "#51A3CC", "S4_10nM"           = "#7EC3E5",
  "S4_1nM"            = "#B2E5FF",
  "Panabinostat_10uM" = "#99540F", "Panabinostat_1uM"  = "#B26F2C",
  "Panabinostat_100nM"= "#CC8E51", "Panabinostat_10nM" = "#E5B17E",
  "Panabinostat_1nM"  = "#FFD8B2",
  "Venetoclax_10uM"   = "#6B990F", "Venetoclax_1uM"    = "#85B22C",
  "Venetoclax_100nM"  = "#A3CC51", "Venetoclax_10nM"   = "#C3E57E",
  "Venetoclax_1nM"    = "#E5FFB2",
  "Rux_10uM"          = "#990F0F", "Rux_1uM"           = "#B22C2C",
  "Rux_100nM"         = "#CC5151", "Rux_10nM"          = "#E57E7E",
  "Rux_1nM"           = "#FFB2B2"
)


# --- 1. LOAD DATA ------------------------------------------------------------
count       <- read_excel(PATH_DRUG_DATA, sheet = "Count")
frequencies <- read_excel(PATH_DRUG_DATA, sheet = "Frequencies")

# Reshape to long format and average technical replicates
count <- melt(count, id.vars = c("Donor", "Drug", "Concentration", "Day")) %>%
  group_by(Donor, Drug, Concentration, Day, variable) %>%
  summarise(value = mean(value), .groups = "drop")

frequencies <- melt(frequencies, id.vars = c("Donor", "Drug", "Concentration", "Day")) %>%
  group_by(Donor, Drug, Concentration, Day, variable) %>%
  summarise(value = mean(value), .groups = "drop")

# Use frequencies for downstream analysis
count <- frequencies


# --- 2. NORMALISE TO DMSO ----------------------------------------------------
count <- filter(count, variable != "Clec9A")

count_normalised <- filter(count, Drug != "No DMSO")
DMSO <- filter(count, Drug == "DMSO")

count_normalised <- unite(count_normalised, Label,
                          Donor, Concentration, Day, variable, sep = "_", remove = FALSE)
DMSO <- unite(DMSO, Label, Donor, Concentration, Day, variable, sep = "_")
DMSO <- dplyr::select(DMSO, -Drug)
names(DMSO) <- c("Label", "DMSO_value")

count_normalised <- left_join(count_normalised, DMSO, by = "Label") %>%
  mutate(normalised_value = value / DMSO_value)


# --- 3. HEATMAP: RUXOLITINIB (log2 scale) ------------------------------------
cell_order <- c("B_cells","T_cells","NK_cells","cDC1","cDC2","pDC",
                "Mast_cells","Monocytes","Neutrophils","CD34","Ery")
conc_order  <- c("10uM_3","10uM_20","10uM_22","1uM_3","1uM_20","1uM_22",
                 "100nM_3","100nM_20","100nM_22","10nM_3","10nM_20","10nM_22",
                 "1nM_3","1nM_20","1nM_22","0_3","0_20","0_22")

count_heatmap <- unite(count_normalised, Conc, Concentration, Donor, sep = "_")
count_heatmap$Conc     <- factor(count_heatmap$Conc,     levels = conc_order)
count_heatmap$variable <- factor(count_heatmap$variable, levels = cell_order)

count_heatmap_r <- filter(count_heatmap, Drug == "Rux")

pdf(file.path(OUT_DIR, "Heatmap_Ruxolitinib.pdf"), height = 8, width = 8)
ggplot(count_heatmap_r, aes(Conc, variable, fill = log2(normalised_value + 1))) +
  geom_tile() +
  scale_fill_viridis(option = "inferno", discrete = FALSE) +
  facet_wrap(~Day, ncol = 2) +
  THEME
dev.off()


# --- 4. HEATMAP: PANABINOSTAT (diverging colour scale centred on 1) ----------
max_dist <- max(abs(log2(count_heatmap$normalised_value + 1) - 1), na.rm = TRUE)

count_heatmap <- count_heatmap %>%
  mutate(
    log_values  = log2(normalised_value + 1),
    value_scaled = rescale(log_values,
                           from = c(1 - max_dist, 1 + max_dist),
                           to   = c(0, 1))
  )

count_heatmap_p <- filter(count_heatmap, Drug == "Panabinostat")
my_palette <- rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30))
max_dist_p <- max(abs(count_heatmap_p$value_scaled - 0.5), na.rm = TRUE)

pdf(file.path(OUT_DIR, "Heatmap_Panabinostat_scaled.pdf"), height = 8, width = 8)
ggplot(count_heatmap_p, aes(Conc, variable, fill = value_scaled)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = my_palette,
    limits  = c(0.5 - max_dist_p, 0.5 + max_dist_p),
    oob     = scales::squish,
    name    = "Fold change\nvs DMSO"
  ) +
  facet_wrap(~Day, ncol = 2) +
  THEME
dev.off()


# --- 5. LINE PLOTS: NORMALISED VALUE BY CONCENTRATION ------------------------
count_normalised$Donor         <- as.character(count_normalised$Donor)
count_normalised$Concentration <- factor(count_normalised$Concentration,
                                         levels = c("10uM","1uM","100nM","10nM","1nM"))

pdf(file.path(OUT_DIR, "line_plots_normalised.pdf"), height = 10, width = 20)
for (d in unique(count$Day)) {
  p <- count_normalised %>%
    filter(Day == d) %>%
    ggplot(aes(x = Concentration, y = log2(normalised_value + 1))) +
    geom_point(aes(colour = variable, shape = Donor), size = 3, alpha = 0.5) +
    stat_summary(aes(y = log2(normalised_value + 1), group = variable, colour = variable),
                 fun = mean, geom = "point", shape = 20, size = 5) +
    stat_summary(aes(y = log2(normalised_value + 1), group = variable, colour = variable),
                 fun = mean, geom = "line", linewidth = 0.7) +
    stat_summary(aes(y = log2(normalised_value + 1), group = variable),
                 fun.data = mean_se, geom = "errorbar",
                 width = 0.25, linewidth = 0.5, alpha = 0.7) +
    facet_wrap(~Drug, scales = "free") +
    ggtitle(paste("Day", d)) +
    THEME
  print(p)
}
dev.off()


# --- 6. LINE PLOTS: VENETOCLAX (raw frequency over time) ---------------------
count$Donor <- as.factor(count$Donor)
count$Day   <- as.factor(count$Day)

high_low <- unite(count, Condition, Drug, Concentration, sep = "_", remove = FALSE) %>%
  filter(!Condition %in% c("DMSO_1uM", "DMSO_100nM", "DMSO_10nM", "No DMSO_0"))

Ven <- filter(high_low, grepl("Venetoclax|DMSO", Drug))

pdf(file.path(OUT_DIR, "Venetoclax_perc.pdf"), height = 6, width = 9)
ggplot(Ven, aes(x = Day, y = value)) +
  geom_point(aes(colour = Condition, shape = Donor), size = 2, alpha = 0.5) +
  scale_colour_manual(values = drug_col) +
  stat_summary(aes(y = value, group = Condition, colour = Condition),
               fun = mean, geom = "point", shape = 20, size = 4) +
  stat_summary(aes(y = value, group = Condition, colour = Condition),
               fun = mean, geom = "line", linewidth = 0.7) +
  stat_summary(aes(y = value, group = Condition),
               fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  THEME
dev.off()

message("Script 04 complete. Output saved to: ", OUT_DIR)
