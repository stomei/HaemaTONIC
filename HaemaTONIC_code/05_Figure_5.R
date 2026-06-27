# =============================================================================
# Script: 05_Figure_5.R
# Description: CRISPR knockout phenotypic analysis for Figure 5
#              Loads flow cytometry counts from CRISPR KO experiments,
#              normalises to scramble control, and generates bar graphs
#              Also includes SCID experiment T cell analysis
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
library(ggbreak)


# --- INPUT PATHS (update these to your local paths) --------------------------
PATH_CRISPR_KO  <- "path/to/Figure5/CRISPR_KO.xlsx"
PATH_SCID_PHENO <- "path/to/Figure5/SCID_phenotyping.xlsx"
PATH_SCID_CORR  <- "path/to/Figure5/SCID_corrected.xlsx"

# Output directory
OUT_DIR <- "output/05_Figure5"
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

# Colour palette: all KO conditions in red, SC control in grey
condition_palette <- c(
  "BCL11B"     = "#A50021",
  "GATA2"      = "#A50021",
  "IRF8"       = "#A50021",
  "KLF1"       = "#A50021",
  "RUNX1"      = "#A50021",
  "SC control" = "grey70",
  "SPI1"       = "#A50021"
)

# Cell type display order
cell_order <- c("RBC","MK_platelets","Monocytes","Neutrophils",
                "cDC1","cDC2","pDC","Mast_cells",
                "NK","T_cells","B_cells","Progenitors")


# --- 1. LOAD AND RESHAPE DATA ------------------------------------------------
count       <- read_excel(PATH_CRISPR_KO, sheet = "Count")
frequencies <- read_excel(PATH_CRISPR_KO, sheet = "Frequencies")

count <- melt(count, id.vars = c("Donor", "Condition", "Day"))


# --- 2. NORMALISE TO SCRAMBLE CONTROL ----------------------------------------
sc_control_means <- count %>%
  filter(Condition == "SC control") %>%
  group_by(variable, Donor, Day) %>%
  summarise(sc_mean = mean(value, na.rm = TRUE), .groups = "drop")

data_normalised <- count %>%
  left_join(sc_control_means, by = c("variable", "Donor", "Day")) %>%
  mutate(value_normalised = value / sc_mean)

# Clean up edge cases
data_normalised$value_normalised[is.infinite(data_normalised$value_normalised)] <- 0
data_normalised[is.na(data_normalised)] <- 0

# Filter and factor
data_normalised <- data_normalised %>%
  filter(variable != "WBC", !is.na(variable))

data_normalised$Day      <- as.factor(data_normalised$Day)
data_normalised$variable <- factor(data_normalised$variable, levels = cell_order)


# --- 3. SPLIT BY EXPERIMENT --------------------------------------------------
# IRF8 KO was performed in a separate experiment (different cord blood donors)
donors_all  <- c("CB023", "CB027", "CB028")
donors_irf8 <- c("CB020", "CB021", "CB022")

data_normalised_all  <- filter(data_normalised,
                               Donor %in% donors_all)
data_normalised_irf8 <- filter(data_normalised,
                               Donor %in% donors_irf8)

# Subset by KO gene
IRF8   <- filter(data_normalised_irf8, grepl("IRF8|SC",   Condition))
BCL11B <- filter(data_normalised_all,  grepl("BCL11B|SC", Condition))
GATA2  <- filter(data_normalised_all,  grepl("GATA2|SC",  Condition))
KLF1   <- filter(data_normalised_all,  grepl("KLF1|SC",   Condition))
RUNX1  <- filter(data_normalised_all,  grepl("RUNX1|SC",  Condition))
SPI1   <- filter(data_normalised_all,  grepl("SPI1|SC",   Condition))

# NOTE: CB027 excluded from BCL11B and RUNX1 due to low KO efficiency
BCL11B <- filter(BCL11B, Donor != "CB027")
RUNX1  <- filter(RUNX1,  Donor != "CB027")


# --- 4. HELPER: BAR GRAPH WITH INDIVIDUAL DATA POINTS -----------------------
plot_ko_bar <- function(data, ko_data, day_filter, gene_name,
                        palette, out_file, height = 6, width = 8) {
  data_summary <- data %>%
    filter(Day == day_filter) %>%
    group_by(Day, Condition, variable) %>%
    summarise(
      Mean = mean(value_normalised),
      SEM  = sd(value_normalised) / sqrt(n()),
      .groups = "drop"
    )

  pdf(out_file, height = height, width = width)
  p <- ggplot(data_summary, aes(x = variable, y = Mean, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             color = "black", width = 0.8) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                  position = position_dodge(width = 0.9), width = 0.25) +
    geom_point(data = filter(ko_data, Day == day_filter),
               aes(x = variable, y = value_normalised, shape = as.factor(Donor)),
               position = position_dodge(width = 0.9),
               size = 1.5, color = "black") +
    scale_fill_manual(values = palette) +
    scale_y_break(breaks = c(4.5, 10), scales = c(0.8, 0.2)) +
    ggtitle(paste(gene_name, "KO — Day", day_filter)) +
    THEME
  print(p)
  dev.off()
}

plot_ko_bar(SPI1, SPI1, day_filter = "14", gene_name = "SPI1",
            palette  = condition_palette,
            out_file = file.path(OUT_DIR, "SPI1_bar_day14.pdf"))


# --- 5. SCID EXPERIMENT — T CELL FREQUENCY OVER TIME ------------------------
pheno <- read_excel(PATH_SCID_PHENO, sheet = "Count")
pheno <- melt(pheno, id.vars = c("Donor", "Day", "Mesure"))

pheno_freq <- filter(pheno, Mesure == "freq")
pheno_freq$variable <- factor(pheno_freq$variable,
                               levels = c("Mono","Neutrophils","CD34",
                                          "cDC1","cDC2","pDC",
                                          "B_cells","T_cells","NK cells"))

corrected <- read_excel(PATH_SCID_CORR, sheet = "Count")
corrected <- dplyr::select(corrected, T_cells, Cell, Day, Mesure) %>%
  filter(!Cell %in% c("BE", "WT"))
corrected <- melt(corrected, id.vars = c("Cell", "Day", "Mesure"))
corrected_freq <- filter(corrected, Mesure == "freq")

condition_palette_scid <- c(
  "SCID" = "#990F0F",
  "WT"   = "#99700F",
  "Nuc"  = "#1F990F",
  "IL2"  = "#710F99"
)

data_summary_scid <- corrected_freq %>%
  group_by(Day, Cell) %>%
  summarise(
    Mean = mean(value),
    SEM  = sd(value) / sqrt(n()),
    .groups = "drop"
  )

pdf(file.path(OUT_DIR, "SCID_Tcell_freq.pdf"), height = 3, width = 5)
ggplot(data_summary_scid, aes(x = as.factor(Day), y = Mean, fill = Cell)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_point(data = corrected_freq,
             aes(x = as.factor(Day), y = value),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black") +
  scale_fill_manual(values = condition_palette_scid) +
  THEME
dev.off()

message("Script 05 complete. Output saved to: ", OUT_DIR)
