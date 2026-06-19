setwd("/vast/projects/ST227/Sara/Christina")
library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)

THEME=theme(text = element_text(size = 15, colour = "black"), 
            plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
            axis.text.y = element_text(colour = "black"), 
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "right"
)

group_palette= c("DNMT3A R736H"= "#CC5800","DNMT3A R882H"= "#FFAD65","TET2 H1904R"= "#1E8E99", "TET2 Q888X"= "#99F9FF", "UBA1 M41L"= "#860086", "UBA1 M41V"= "#F100F1", "PE Control"= "grey80")
CL014 <- read_excel("CL014.xlsx", sheet = "Count")

count <- melt(CL014, id.vars = c("Sample","Group", "Well", "Day", "Donor"))


#Average technical replicates for PE control
pe_control_means <- count %>%
  filter(Group == "PE Control") %>%
  group_by(variable, Donor, Day) %>%
  summarise(pe_mean = mean(value, na.rm = TRUE), .groups = "drop")

#Join PE control means back to the full dataset and normalise
data_normalised <- count %>%
  left_join(pe_control_means, by = c("variable", "Donor", "Day")) %>%
  mutate(value_normalised = value / pe_mean)

data_normalised$value_normalised[is.infinite(data_normalised$value_normalised)]=0
data_normalised[is.na(data_normalised)]=0
data_normalised$Day= as.factor(data_normalised$Day)
data_normalised$Group = factor(data_normalised$Group , levels = c("DNMT3A R736H","DNMT3A R882H","TET2 H1904R", "TET2 Q888X", "UBA1 M41L", "UBA1 M41V", "PE Control"))
#write.table(data_normalised, "CL014_normalised.txt", sep="\t", row.names=FALSE, col.names=TRUE)
#Plot
ggplot(data_normalised[!grepl("UBA1", data_normalised$Group),], aes(x=Day, y=log2(value_normalised+1))) +
  geom_point(aes(colour= Group, shape= as.factor(Donor)), alpha=0.5, size= 3) +
  #scale_colour_manual(values = pop_palette)+
  stat_summary(aes(y = log2(value_normalised), group = Group, colour= Group), fun.y=mean, geom="point", shape=20, size=4) +
  stat_summary(aes(y = log2(value_normalised), group = Group, colour= Group), fun.y=mean, geom="line", linewidth= 0.7)+
  stat_summary(aes(y = log2(value_normalised), group = Group), fun.data = mean_se, geom='errorbar', width = 0.25, size = 0.5, alpha = 0.7, linetype = "solid")+
  facet_wrap(~variable, scales = "free")+
  THEME

DNMT3A= data_normalised[grepl("DNMT3A", data_normalised$Group) |grepl("PE", data_normalised$Group),]
TET2= data_normalised[grepl("TET2", data_normalised$Group) |grepl("PE", data_normalised$Group),]
UBA1= data_normalised[grepl("UBA1", data_normalised$Group) |grepl("PE", data_normalised$Group),]

pdf("UBA1_fold_change_numb.pdf", height = 6, width = 8)
ggplot(UBA1, aes(x=Day, y=value_normalised)) +
  geom_point(aes(colour= Group, shape= as.factor(Donor)), alpha=0.5, size= 2) +
  scale_colour_manual(values = group_palette)+
  stat_summary(aes(y = value_normalised, group = Group, colour= Group), fun.y=mean, geom="point", shape=20, size=3) +
  stat_summary(aes(y = value_normalised, group = Group, colour= Group), fun.y=mean, geom="line", linewidth= 0.7)+
  stat_summary(aes(y = value_normalised, group = Group), fun.data = mean_se, geom='errorbar', width = 0.25, size = 0.5, alpha = 0.7, linetype = "solid")+
  facet_wrap(~variable, scales = "free")+
  THEME
dev.off()

#Heatmap
#Average across donors, Group and Day
data_avg <- data_normalised %>%
  group_by(variable, Donor, Group, Day) %>%
  summarise(mean_norm = mean(value_normalised, na.rm = TRUE), .groups = "drop")

data_avg <- data_avg %>%
  mutate(col_label = paste(Group, "Donor", Donor))

data_avg <- data_avg %>%
  group_by(Donor, Group, variable) %>%
  mutate(colnorm = mean_norm / sum(mean_norm, na.rm = TRUE)*100) %>%
  ungroup()

data_avg$mean_norm[is.infinite(data_avg$mean_norm)]=0
data_avg$col_label = factor(data_avg$col_label, levels = c(
    "DNMT3A R736H Donor 1", "DNMT3A R736H Donor 2", "DNMT3A R736H Donor 3",
    "DNMT3A R882H Donor 1", "DNMT3A R882H Donor 2", "DNMT3A R882H Donor 3",
    "TET2 H1904R Donor 1",  "TET2 H1904R Donor 2",  "TET2 H1904R Donor 3",
    "TET2 Q888X Donor 1",   "TET2 Q888X Donor 2",   "TET2 Q888X Donor 3",
    "UBA1 M41L Donor 1",    "UBA1 M41L Donor 2",     "UBA1 M41L Donor 3",
    "UBA1 M41V Donor 1",    "UBA1 M41V Donor 2",     "UBA1 M41V Donor 3",
    "Wild Type Donor 1",    "Wild Type Donor 2",      "Wild Type Donor 3",
    "PE Control Donor 1",   "PE Control Donor 2",     "PE Control Donor 3"
  ))
ggplot(data_avg, aes(x = col_label, y = variable, fill = log2(mean_norm+1))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", name = "% of total") +
  THEME +
  facet_wrap(~Day)

#Bar graph
data_summary <- TET2 %>%
  group_by(Day, Group, variable) %>%
  summarise(
    Mean = mean(value_normalised),
    SEM  = sd(value_normalised) / sqrt(n()),
    .groups = "drop"
  )

data_summary$Group = factor(data_summary$Group , levels = c("TET2 H1904R", "TET2 Q888X", "PE Control"))

pdf("TET2_bar_graph_numb.pdf", height = 6, width = 8)
ggplot(data_summary, aes(x = Day, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_point(data = TET2, aes(x = Day, y = value_normalised, group = Group, shape= as.factor(Donor)),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black")  +
  scale_fill_manual(values = group_palette)+
  facet_wrap(~variable, scale="free")+
  THEME
dev.off()  

#t-test
ttest_results <- data_normalised %>%
  group_by(Day, variable) %>%
  summarise(
    # DNMT3A R736H vs PE control
    p_R763H_vs_PE  = t.test(value_normalised[Group == "DNMT3A R736H"],
                          value_normalised[Group == "PE Control"])$p.value,
    
    # DNMT3A R882H vs PE control
    p_R882H_vs_PE  = t.test(value_normalised[Group == "DNMT3A R882H"],
                          value_normalised[Group == "PE Control"])$p.value,
    
    # TET2 H1904R vs PE control
    p_H1904R_vs_PE  = t.test(value_normalised[Group == "TET2 H1904R"],
                            value_normalised[Group == "PE Control"])$p.value,

    # TET2 Q888X vs PE control
    p_Q888X_vs_PE  = t.test(value_normalised[Group == "TET2 Q888X"],
                            value_normalised[Group == "PE Control"])$p.value,
   
    # UBA1 M41L vs PE control
    p_M41L_vs_PE  = t.test(value_normalised[Group == "UBA1 M41L"],
                            value_normalised[Group == "PE Control"])$p.value,
    
    # UBA1 M41V vs PE control
    p_M41V_vs_PE  = t.test(value_normalised[Group == "UBA1 M41V"],
                            value_normalised[Group == "PE Control"])$p.value,
  )

print(ttest_results)
write.table(ttest_results, "t-test_PE.txt", sep="\t", row.names=FALSE, col.names=TRUE)
#Stimuli
stimuli <- read_excel("/vast/projects/ST227/Sara/Christina/stimuli.xlsx")
stimuli <- melt(stimuli, id.vars = c("sample", "donor", "stimuli"))
stimuli$sample = factor(stimuli$sample , levels = c("Control","R736H","R882H", "H1904R", "Q888X", "M41L", "M41V"))
stimuli$stimuli = factor(stimuli$stimuli , levels = c("unstimulated","stimulated"))


group_palette= c("R736H"= "#CC5800","R882H"= "#FFAD65","H1904R"= "#1E8E99", "Q888X"= "#99F9FF", "M41L"= "#860086", "M41V"= "#F100F1", "Control"= "grey80")

#Bar graph
data_summary <- stimuli %>%
  group_by(sample, donor, stimuli, variable) %>%
  summarise(
    Mean = mean(value),
    SEM  = sd(value) / sqrt(n()),
    .groups = "drop"
  )

pdf("cytokines.pdf", height = 8, width = 9)
ggplot(data_summary, aes(x = stimuli, y = Mean, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_point(data = stimuli, aes(x = stimuli, y = value, group = sample, shape= as.factor(donor)),
             position = position_dodge(width = 0.9),
             size = 1.5, color = "black")  +
  scale_fill_manual(values = group_palette)+
  facet_wrap(~variable, scale="free")+
  THEME
dev.off()  
x <- stimuli
#t-test
ttest_results <- x %>%
  group_by(sample, variable, stimuli) %>%
  summarise(
    # DNMT3A R736H vs PE control
    p_R763H_vs_PE  = t.test(value[sample == "R736H"],
                            value[sample == "Control"])$p.value,
    
    # DNMT3A R882H vs PE control
    p_R882H_vs_PE  = t.test(value[sample == "R882H"],
                            value[sample == "Control"])$p.value,
    
    # TET2 H1904R vs PE control
    p_H1904R_vs_PE  = t.test(value[sample == "H1904R"],
                             value[sample == "Control"])$p.value,
    
    # TET2 Q888X vs PE control
    p_Q888X_vs_PE  = t.test(value[sample == "Q888X"],
                            value[sample == "Control"])$p.value,
    
    # UBA1 M41L vs PE control
    p_M41L_vs_PE  = t.test(value[sample == "M41L"],
                           value[sample == "Control"])$p.value,
    
    # UBA1 M41V vs PE control
    p_M41V_vs_PE  = t.test(value[sample == "M41V"],
                           value[sample == "Control"])$p.value,
  )

print(ttest_results)
write.table(ttest_results, "t-test_PE.txt", sep="\t", row.names=FALSE, col.names=TRUE)

ttest_results <- stimuli %>%
  group_by(sample, variable, stimuli) %>%
  summarise(
    p_stim_vs_unstim = t.test(
      value[stimuli == "stimulated"],
      value[stimuli == "unstimulated"]
    )$p.value,
    .groups = "drop"
  )



