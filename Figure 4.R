setwd("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/Drug_Screen")
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(viridisLite)
library(viridis)
library(paletteer)
library(scales)


THEME=theme(text = element_text(size = 15, colour = "black"), 
            plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
            axis.text.y = element_text(colour = "black"), 
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "right"
)


count <- read_excel("ST232_Numbers.xlsx", sheet = "Count")
frequencies <- read_excel("ST232_Numbers.xlsx", sheet = "Frequencies")

count <- melt(count, id.vars = c("Donor", "Drug", "Concentration", "Day"))
count <- count %>% group_by(Donor, Drug, Concentration, Day, variable) %>% summarise(value= mean(value)) %>% ungroup()
frequencies <- melt(frequencies, id.vars = c("Donor", "Drug", "Concentration", "Day"))
frequencies <- frequencies %>% group_by(Donor, Drug, Concentration, Day, variable) %>% summarise(value= mean(value)) %>% ungroup()

count <- frequencies
#Heatmap normalised to DMSO
count <- filter(count, !variable=="Clec9A")
count_normalised <- filter(count, !Drug=="No DMSO")
DMSO <- filter(count, Drug=="DMSO")
count_normalised <- unite(count_normalised, Label, Donor, Concentration, Day, variable, sep = "_", remove = FALSE)
DMSO <- unite(DMSO, Label, Donor, Concentration, Day, variable, sep = "_")
DMSO <- dplyr::select(DMSO, -Drug)
names(DMSO) <- c("Label", "DMSO_value")
count_normalised <- left_join(count_normalised, DMSO, by= "Label")
count_normalised <- mutate(count_normalised, normalised_value= value/DMSO_value)


#Heatmap
count_heatmap <- unite(count_normalised, Conc, Concentration, Donor, sep="_")
count_heatmap_norm <- count_heatmap %>%
  group_by(Day, Drug, variable) %>%
  mutate(colnorm = normalised_value / sum(normalised_value, na.rm = TRUE)*100) %>%
  ungroup()
count_heatmap_norm$Conc <- factor(count_heatmap_norm$Conc, levels = c("10uM_3", "10uM_20", "10uM_22", "1uM_3", "1uM_20", "1uM_22", "100nM_3", "100nM_20", "100nM_22", "10nM_3", "10nM_20", "10nM_22","1nM_3", "1nM_20", "1nM_22", "0_3", "0_20", "0_22"))
count_heatmap_norm$variable <- factor(count_heatmap_norm$variable, levels = c("B_cells", "T_cells", "NK_cells",  "cDC1", "cDC2", "pDC", "Mast_cells", "Monocytes", "Neutrophils", "CD34", "Ery"))
count_heatmap_r <- filter(count_heatmap_norm, Drug=="Rux")
pdf("ST257_Haetmap_count_r.pdf", height = 8, width = 8)
ggplot(count_heatmap_r, aes(Conc, variable, fill= log2(normalised_value+1))) + 
  geom_tile()+
  scale_fill_viridis(option = "inferno", discrete=FALSE)+
  facet_wrap(~ Day, ncol=2)+
  THEME
dev.off()


# Find the maximum distance from 1 in your data
count_heatmap_norm$log_values= log2(count_heatmap_norm$normalised_value+1)
max_dist <- max(abs(count_heatmap_norm$log_values - 1), na.rm = TRUE)

# Rescale so that 1 = 0.5 (middle of the colour scale)
count_heatmap_norm <- count_heatmap_norm %>%
  mutate(value_scaled = rescale(log_values, 
                                from = c(1 - max_dist, 1 + max_dist), 
                                to = c(0, 1)))

count_heatmap_p <- filter(count_heatmap_norm, Drug=="Panabinostat")
my_palette <- rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30))
max_dist <- max(abs(count_heatmap_p$value_scaled - 0.5), na.rm = TRUE)

pdf("ST257_Haetmap_count_scaled_p.pdf", height = 8, width = 8)
ggplot(count_heatmap_p, aes(Conc, variable, fill= value_scaled)) + 
  geom_tile()+
  scale_fill_gradientn(
    colours = my_palette,
    limits = c(0.5 - max_dist, 0.5 + max_dist),  # symmetric around 0.5
    oob = scales::squish,                          # squish any out of range values
    name = "Fold change\nvs DMSO"
  ) +
  facet_wrap(~ Day, ncol=2)+
  THEME
dev.off()

#line plot
count$Donor <- as.character(count$Donor)
count$Concentration <- factor(count$Concentration, levels = c("10uM", "1uM", "100nM", "10nM", "1nM", "0"))

pdf("line_plots_frequencies.pdf", height = 10, width = 20) 
for(d in unique(count$Day)) {
  p <- count %>%
    filter(Day == d) %>%
    ggplot(aes(x=Concentration, y=normalised_value)) +
    geom_point(aes(colour= Drug, shape= Donor, size= Donor), alpha=0.5) +
    #scale_colour_manual(values = CD200_col)+
    #scale_shape_manual(values=c(15, 17, 18))+
    #scale_size_manual(values=c(3,3,4))+
    stat_summary(aes(y = value, group = Drug, colour= Drug), fun.y=mean, geom="point", shape=20, size=5) +
    stat_summary(aes(y = value, group = Drug, colour= Drug), fun.y=mean, geom="line", linewidth= 0.7)+
    stat_summary(aes(y = value, group = Drug), fun.data = mean_se, geom='errorbar', width = 0.25, size = 0.5, alpha = 0.7, linetype = "solid")+
    #stat_compare_means(method = "t.test", label = "p.signif", label.y = c(22, 22.5, 23, 23.5, 24, 24.5), comparisons = y, hide.ns = FALSE)+
    #ylim(0,25)+
    facet_wrap(~variable, scales = "free")+
    THEME
  print(p)  # print to PDF (each print = new page)
}
dev.off()

count_normalised$Donor <- as.character(count_normalised$Donor)
pdf("line_plots_frequencies.pdf", height = 10, width = 20) 
count_normalised$Concentration <- factor(count_normalised$Concentration, levels = c("10uM", "1uM", "100nM", "10nM", "1nM"))
for(d in unique(count$Day)) {
  p <- count_normalised %>%
    filter(Day == d) %>%
    ggplot(aes(x=Concentration, y=log2(normalised_value+1))) +
    geom_point(aes(colour= variable, shape= Donor),  size= 3, alpha=0.5) +
    #scale_colour_manual(values = CD200_col)+
    #scale_shape_manual(values=c(15, 17, 18))+
    #scale_size_manual(values=c(3,3,4))+
    stat_summary(aes(y = log2(normalised_value+1), group = variable, colour= variable), fun.y=mean, geom="point", shape=20, size=5) +
    stat_summary(aes(y = log2(normalised_value+1), group = variable, colour= variable), fun.y=mean, geom="line", linewidth= 0.7)+
    stat_summary(aes(y = log2(normalised_value+1), group = variable), fun.data = mean_se, geom='errorbar', width = 0.25, size = 0.5, alpha = 0.7, linetype = "solid")+
    #stat_compare_means(method = "t.test", label = "p.signif", label.y = c(22, 22.5, 23, 23.5, 24, 24.5), comparisons = y, hide.ns = FALSE)+
    #ylim(0,25)+
    facet_wrap(~Drug, scales = "free")+
    THEME
  print(p)  # print to PDF (each print = new page)
}
dev.off()

#line plot with 2 concentrations
count$Donor= as.factor(count$Donor)
count$Day= as.factor(count$Day)
high_low <- unite(count, Condition, Drug, Concentration, sep = "_", remove = FALSE)
high_low <- filter(high_low, !Condition == "DMSO_1uM" & !Condition =="DMSO_100nM" & !Condition =="DMSO_10nM" & !Condition =="No DMSO_0")
Pana= high_low[grepl("Panabinostat", high_low$Drug) |grepl("DMSO", high_low$Drug),]
Rux= high_low[grepl("Rux", high_low$Drug) |grepl("DMSO", high_low$Drug),]
S4= high_low[grepl("S4", high_low$Drug) |grepl("DMSO", high_low$Drug),]
Ven= high_low[grepl("Venetoclax", high_low$Drug) |grepl("DMSO", high_low$Drug),]

drug_col= c("DMSO_10uM"= "grey60", "DMSO_1nM"= "grey90", 
            "S4_10uM"= "#0F6B99", "S4_1uM"= "#2C85B2", "S4_100nM"= "#51A3CC", "S4_10nM"= "#7EC3E5", "S4_1nM"= "#B2E5FF", 
            "Panabinostat_10uM"= "#99540F", "Panabinostat_1uM"= "#B26F2C", "Panabinostat_100nM"= "#CC8E51", "Panabinostat_10nM"= "#E5B17E", "Panabinostat_1nM"="#FFD8B2", 
            "Venetoclax_10uM"= "#6B990F", "Venetoclax_1uM"= "#85B22C", "Venetoclax_100nM"= "#A3CC51", "Venetoclax_10nM"= "#C3E57E", "Venetoclax_1nM"= "#E5FFB2", 
            "Rux_10uM"= "#990F0F", "Rux_1uM"= "#B22C2C", "Rux_100nM"= "#CC5151", "Rux_10nM"= "#E57E7E", "Rux_1nM"= "#FFB2B2")


pdf("Venetoclax_perc.pdf", height = 6, width = 9)
ggplot(Ven, aes(x=Day, y=value)) +
  geom_point(aes(colour= Condition, shape= Donor),  size= 2, alpha=0.5) +
  scale_colour_manual(values = drug_col)+
  #scale_shape_manual(values=c(15, 17, 18))+
  #scale_size_manual(values=c(3,3,4))+
  stat_summary(aes(y = value, group = Condition, colour= Condition), fun.y=mean, geom="point", shape=20, size=4) +
  stat_summary(aes(y = value, group = Condition, colour= Condition), fun.y=mean, geom="line", linewidth= 0.7)+
  stat_summary(aes(y = value, group = Condition), fun.data = mean_se, geom='errorbar', width = 0.25, size = 0.5, alpha = 0.7, linetype = "solid")+
  #stat_compare_means(method = "t.test", label = "p.signif", label.y = c(22, 22.5, 23, 23.5, 24, 24.5), comparisons = y, hide.ns = FALSE)+
  #ylim(0,25)+
  facet_wrap(~variable, scales = "free")+
  THEME
dev.off()