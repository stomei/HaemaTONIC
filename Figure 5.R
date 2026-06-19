install.packages("ggbreak")
setwd("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/CRISPR_KO")
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


THEME=theme(text = element_text(size = 15, colour = "black"), 
            plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
            axis.text.y = element_text(colour = "black"), 
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "right"
)

condition_palette= c("BCL11B"= "#A50021", "GATA2"= "#A50021", "IRF8"= "#A50021", "KLF1"= "#A50021", "RUNX1"= "#A50021", "SC control"= "grey70", "SPI1"= "#A50021")
count <- read_excel("ML099_CRISPR_KO.xlsx", sheet = "Count")
frequencies <- read_excel("ML099_CRISPR_KO.xlsx", sheet = "Frequencies")

count <- melt(count, id.vars = c("Donor","Condition", "Day"))

#Average technical replicates for PE control
sc_control_means <- count %>%
  filter(Condition == "SC control") %>%
  group_by(variable, Donor, Day) %>%
  summarise(sc_mean = mean(value, na.rm = TRUE), .groups = "drop")

#Join PE control means back to the full dataset and normalise
data_normalised <- count %>%
  left_join(sc_control_means, by = c("variable", "Donor", "Day")) %>%
  mutate(value_normalised = value / sc_mean)

data_normalised$value_normalised[is.infinite(data_normalised$value_normalised)]=0
data_normalised[is.na(data_normalised)]=0
data_normalised= filter(data_normalised, !variable =="WBC")
data_normalised <- filter(data_normalised, !is.na(variable))
data_normalised$Day= as.factor(data_normalised$Day)
data_normalised$variable = factor(data_normalised$variable , levels = c("RBC","MK_platelets", "Monocytes","Neutrophils", "cDC1", "cDC2", "pDC", "Mast_cells", "NK", "T_cells", "B_cells", "Progenitors"))

#Divide SC and WT control by experiments
data_normalised_all= data_normalised[grepl("CB023", data_normalised$Donor) | grepl("CB027", data_normalised$Donor) | grepl("CB028", data_normalised$Donor),]
data_normalised_irf8= data_normalised[grepl("CB020", data_normalised$Donor) | grepl("CB021", data_normalised$Donor) | grepl("CB022", data_normalised$Donor),]

IRF8= data_normalised_irf8[grepl("IRF8", data_normalised_irf8$Condition) |grepl("SC", data_normalised_irf8$Condition),]
BCL11B= data_normalised_all[grepl("BCL11B", data_normalised_all$Condition) |grepl("SC", data_normalised_all$Condition),]
GATA2= data_normalised_all[grepl("GATA2", data_normalised_all$Condition) |grepl("SC", data_normalised_all$Condition),]
KLF1= data_normalised_all[grepl("KLF1", data_normalised_all$Condition) |grepl("SC", data_normalised_all$Condition),]
RUNX1= data_normalised_all[grepl("RUNX1", data_normalised_all$Condition) |grepl("SC", data_normalised_all$Condition),]
SPI1= data_normalised_all[grepl("SPI1", data_normalised_all$Condition) |grepl("SC", data_normalised_all$Condition),]

#Remove donor CB027 from BCL11B and RUNX1 because KO effieciency is too low
BCL11B= filter(BCL11B, !Donor == "CB027")
RUNX1= filter(RUNX1, !Donor == "CB027")

#Bar graph
data_summary <- SPI1 %>%
  group_by(Day, Condition, variable) %>%
  summarise(
    Mean = mean(value_normalised),
    SEM  = sd(value_normalised) / sqrt(n()),
    .groups = "drop"
  )

#data_summary$Group = factor(data_summary$Group , levels = c("TET2 H1904R", "TET2 Q888X", "PE Control"))

pdf("TET2_bar_graph_numb.pdf", height = 6, width = 8)
ggplot(filter(data_summary, Day== "14"), aes(x = variable, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_point(data = filter(SPI1, Day== "14"), aes(x = variable, y = value_normalised, Condition = Condition, shape= as.factor(Donor)),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black")  +
  scale_fill_manual(values = condition_palette)+
  scale_y_break(breaks = c(4.5, 10),   # where to cut
                scales = c(0.8, 0.2)) + # 70% lower, 30% upper
  #facet_wrap(~variable)+
  THEME
dev.off()  


#SCID
pheno <- read_excel("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/CRISPR_KO/ST235_Numbers.xlsx", sheet = "Count")
corrected <- read_excel("/vast/projects/ST227/Sara/HaemaTONIC_protocol_manuscript/CRISPR_KO/ST258_Numbers.xlsx", sheet = "Count")

pheno <- melt(pheno, id.vars = c("Donor", "Day", "Mesure"))
pheno_numb <- filter(pheno, Mesure=="count")
pheno_freq <- filter(pheno, Mesure=="freq")
pheno_numb$variable = factor(pheno_numb$variable , levels = c("Mono", "Neutrophils", "CD34", "cDC1", "cDC2", "pDC", "B_cells", "T_cells", "NK cells" ))
pheno_freq$variable = factor(pheno_freq$variable , levels = c("Mono", "Neutrophils", "CD34", "cDC1", "cDC2", "pDC", "B_cells", "T_cells", "NK cells" ))

corrected <- select(corrected, T_cells, Cell, Day, Mesure)
corrected <- filter(corrected, !Cell =="BE")
corrected <- filter(corrected, !Cell =="WT")
corrected <- melt(corrected, id.vars = c("Cell", "Day", "Mesure"))
corrected_numb <- filter(corrected, Mesure=="count")
corrected_freq <- filter(corrected, Mesure=="freq")

#Bar graph
data_summary <- corrected_freq %>%
  group_by(Day, Cell) %>%
  summarise(
    Mean = mean(value),
    SEM  = sd(value) / sqrt(n()),
    .groups = "drop"
  )

condition_palette= c("SCID"= "#990F0F", "WT"= "#99700F", "Nuc"= "#1F990F", "IL2"= "#710F99")
pheno_numb$variable = factor(pheno_numb$variable , levels = c("Mono", "Neutrophils", "CD34", "cDC1", "cDC2", "pDC", "B_cells", "T_cells", "NK cells" ))


pdf("corrected_freq.pdf", height = 3, width = 5)
ggplot(data_summary, aes(x = as.factor(Day), y = Mean, fill = Cell)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_point(data= corrected_freq, aes(x = as.factor(Day), y = value),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black")  +
  scale_fill_manual(values = condition_palette)+
  #scale_y_break(breaks = c(7000, 7000),   # where to cut
               # scales = c(0.8, 0.2)) + # 70% lower, 30% upper
  #facet_wrap(~variable)+
  THEME
dev.off()




