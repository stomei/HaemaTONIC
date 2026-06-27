# =============================================================================
# Script: 06_Figure_6.R
# Description: Clonal haematopoiesis (CH) mutation analysis for Figure 6
#              Loads flow cytometry counts from base-edited CH mutation
#              experiments (DNMT3A, TET2, UBA1 variants), normalises to
#              PE control, generates line plots, heatmaps, bar graphs,
#              and performs t-tests vs PE control
#              Also analyses cytokine secretion data from stimulated CH cells
# HaemaTONIC: a multi-lineage model of human haematopoiesis with diverse applications
# Author: Sara Tomei, WEHI
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)


# --- INPUT PATHS (update these to your local paths) --------------------------
PATH_CH_DATA      <- "path/to/PE.xlsx"
PATH_STIMULI_DATA <- "path/to/stimuli.xlsx"

# Output directory
OUT_DIR <- "output/06_Figure6"
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

# Colour palette for CH mutation groups
group_palette <- c(
  "DNMT3A R736H" = "#CC5800",
  "DNMT3A R882H" = "#FFAD65",
  "TET2 H1904R"  = "#1E8E99",
  "TET2 Q888X"   = "#99F9FF",
  "UBA1 M41L"    = "#860086",
  "UBA1 M41V"    = "#F100F1",
  "PE Control"   = "grey80"
)

# Short-name palette for stimuli figure
group_palette_short <- c(
  "R736H"   = "#CC5800",
  "R882H"   = "#FFAD65",
  "H1904R"  = "#1E8E99",
  "Q888X"   = "#99F9FF",
  "M41L"    = "#860086",
  "M41V"    = "#F100F1",
  "Control" = "grey80"
)

# Factor level order for mutation groups
group_levels <- c("DNMT3A R736H","DNMT3A R882H","TET2 H1904R",
                  "TET2 Q888X","UBA1 M41L","UBA1 M41V","PE Control")


# --- 1. LOAD AND RESHAPE DATA ------------------------------------------------
CL014 <- read_excel(PATH_CH_DATA, sheet = "Count")
count <- melt(CL014, id.vars = c("Sample", "Group", "Well", "Day", "Donor"))


# --- 2. NORMALISE TO PE CONTROL ----------------------------------------------
pe_control_means <- count %>%
  filter(Group == "PE Control") %>%
  group_by(variable, Donor, Day) %>%
  summarise(pe_mean = mean(value, na.rm = TRUE), .groups = "drop")

data_normalised <- count %>%
  left_join(pe_control_means, by = c("variable", "Donor", "Day")) %>%
  mutate(value_normalised = value / pe_mean)

# Clean up edge cases
data_normalised$value_normalised[is.infinite(data_normalised$value_normalised)] <- 0
data_normalised[is.na(data_normalised)] <- 0

data_normalised$Day   <- as.factor(data_normalised$Day)
data_normalised$Group <- factor(data_normalised$Group, levels = group_levels)


# --- 3. SUBSET BY MUTATION GROUP ---------------------------------------------
DNMT3A <- filter(data_normalised, grepl("DNMT3A|PE", Group))
TET2   <- filter(data_normalised, grepl("TET2|PE",   Group))
UBA1   <- filter(data_normalised, grepl("UBA1|PE",   Group))


# --- 4. LINE PLOTS -----------------------------------------------------------
pdf(file.path(OUT_DIR, "UBA1_fold_change.pdf"), height = 6, width = 8)
ggplot(UBA1, aes(x = Day, y = value_normalised)) +
  geom_point(aes(colour = Group, shape = as.factor(Donor)),
             alpha = 0.5, size = 2) +
  scale_colour_manual(values = group_palette) +
  stat_summary(aes(y = value_normalised, group = Group, colour = Group),
               fun = mean, geom = "point", shape = 20, size = 3) +
  stat_summary(aes(y = value_normalised, group = Group, colour = Group),
               fun = mean, geom = "line", linewidth = 0.7) +
  stat_summary(aes(y = value_normalised, group = Group),
               fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  THEME
dev.off()


# --- 5. HEATMAP: FOLD CHANGE OVER TIME PER DONOR ----------------------------
data_avg <- data_normalised %>%
  group_by(variable, Donor, Group, Day) %>%
  summarise(mean_norm = mean(value_normalised, na.rm = TRUE), .groups = "drop") %>%
  mutate(col_label = paste(Group, "Donor", Donor))

# Column order: group-then-donor
col_levels <- c(
  paste("DNMT3A R736H Donor", 1:3), paste("DNMT3A R882H Donor", 1:3),
  paste("TET2 H1904R Donor",  1:3), paste("TET2 Q888X Donor",   1:3),
  paste("UBA1 M41L Donor",    1:3), paste("UBA1 M41V Donor",    1:3),
  paste("Wild Type Donor",    1:3), paste("PE Control Donor",   1:3)
)
data_avg$col_label <- factor(data_avg$col_label, levels = col_levels)
data_avg$mean_norm[is.infinite(data_avg$mean_norm)] <- 0

pdf(file.path(OUT_DIR, "heatmap_fold_change.pdf"), height = 6, width = 12)
ggplot(data_avg, aes(x = col_label, y = variable,
                     fill = log2(mean_norm + 1))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", name = "log2(FC+1)") +
  facet_wrap(~Day) +
  THEME +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()


# --- 6. BAR GRAPHS: TET2 MUTATIONS -------------------------------------------
TET2_summary <- TET2 %>%
  group_by(Day, Group, variable) %>%
  summarise(
    Mean = mean(value_normalised),
    SEM  = sd(value_normalised) / sqrt(n()),
    .groups = "drop"
  )
TET2_summary$Group <- factor(TET2_summary$Group,
                             levels = c("TET2 H1904R","TET2 Q888X","PE Control"))

pdf(file.path(OUT_DIR, "TET2_bar_graph.pdf"), height = 6, width = 8)
ggplot(TET2_summary, aes(x = Day, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_point(data = TET2,
             aes(x = Day, y = value_normalised, group = Group,
                 shape = as.factor(Donor)),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black") +
  scale_fill_manual(values = group_palette) +
  facet_wrap(~variable, scales = "free") +
  THEME
dev.off()


# --- 7. T-TESTS: EACH MUTATION GROUP VS PE CONTROL --------------------------
ttest_results <- data_normalised %>%
  group_by(Day, variable) %>%
  summarise(
    p_DNMT3A_R736H_vs_PE = tryCatch(
      t.test(value_normalised[Group == "DNMT3A R736H"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    p_DNMT3A_R882H_vs_PE = tryCatch(
      t.test(value_normalised[Group == "DNMT3A R882H"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    p_TET2_H1904R_vs_PE = tryCatch(
      t.test(value_normalised[Group == "TET2 H1904R"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    p_TET2_Q888X_vs_PE = tryCatch(
      t.test(value_normalised[Group == "TET2 Q888X"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    p_UBA1_M41L_vs_PE = tryCatch(
      t.test(value_normalised[Group == "UBA1 M41L"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    p_UBA1_M41V_vs_PE = tryCatch(
      t.test(value_normalised[Group == "UBA1 M41V"],
             value_normalised[Group == "PE Control"])$p.value,
      error = function(e) NA_real_),

    .groups = "drop"
  )

print(ttest_results)
write.table(ttest_results, file.path(OUT_DIR, "ttest_vs_PE_control.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE)


# --- 8. CYTOKINE SECRETION FROM STIMULATED CH CELLS --------------------------
stimuli <- read_excel(PATH_STIMULI_DATA)
stimuli <- melt(stimuli, id.vars = c("sample", "donor", "stimuli"))

stimuli$sample  <- factor(stimuli$sample,
                           levels = c("Control","R736H","R882H",
                                      "H1904R","Q888X","M41L","M41V"))
stimuli$stimuli <- factor(stimuli$stimuli,
                           levels = c("unstimulated", "stimulated"))

stimuli_summary <- stimuli %>%
  group_by(sample, donor, stimuli, variable) %>%
  summarise(
    Mean = mean(value),
    SEM  = sd(value) / sqrt(n()),
    .groups = "drop"
  )

pdf(file.path(OUT_DIR, "cytokines.pdf"), height = 8, width = 9)
ggplot(stimuli_summary, aes(x = stimuli, y = Mean, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_point(data = stimuli,
             aes(x = stimuli, y = value, group = sample,
                 shape = as.factor(donor)),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black") +
  scale_fill_manual(values = group_palette_short) +
  facet_wrap(~variable, scales = "free") +
  THEME
dev.off()


# --- 9. T-TESTS: STIMULATED VS UNSTIMULATED ----------------------------------
ttest_stim <- stimuli %>%
  group_by(sample, variable, stimuli) %>%
  summarise(
    p_stim_vs_unstim = tryCatch(
      t.test(value[stimuli == "stimulated"],
             value[stimuli == "unstimulated"])$p.value,
      error = function(e) NA_real_),
    .groups = "drop"
  )

write.table(ttest_stim, file.path(OUT_DIR, "ttest_stim_vs_unstim.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE)

message("Script 06 complete. Output saved to: ", OUT_DIR)
