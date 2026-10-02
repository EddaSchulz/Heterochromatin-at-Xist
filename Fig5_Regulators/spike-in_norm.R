## C&T analysis: reads at spike-ins ##
## Script by Nelly Kanata, OWL Schulz ##
## Created on: 11.08.2023
## Modified on: 06.05.2026 ##

library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(egg)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank(),
                  plot.title = element_text(size=8)))



# Read all the .tab files
wdir="./"

setwd(wdir)

### Activators ####
qc_samples_activators <- read.csv("./fig_data/spike-in_activators.txt", sep="\t")
qc_samples_activators %<>%
  separate(Sample, c("cell_line", "timepoint", "replicate", "antibody", "sample_number"))

qc_samples_activators$KD <- if_else(qc_samples_activators$cell_line == "G0", "Rnf12 KD",
                                    if_else(qc_samples_activators$cell_line == "G500", "Rnf12 ctrl", 
                                            if_else(qc_samples_activators$cell_line == "T0", "Zic3 KD","Zic3 ctrl")))

my_comparisons <- list(c("Rnf12 KD", "Rnf12 ctrl"), c("Zic3 KD", "Zic3 ctrl"))

qc_samples_activators$KD <- factor(qc_samples_activators$KD, levels=c("Rnf12 ctrl", "Rnf12 KD", "Zic3 ctrl", "Zic3 KD"))

stat.test <- compare_means(
  AmpR_normalized ~ KD,
  data = qc_samples_activators,
  method = "t.test", 
  paired = TRUE)

# format pvalue
format_p <- function(p) {
  if (p < 0.001) {
    "P < 0.001"
  } else if (p < 0.01) {
    sprintf("P = %.3f", p)
  } else {
    sprintf("P = %.2f", p)
  }
}

stat.test$pval <- vapply(stat.test$p, format_p, character(1))

stat.test %<>%
  filter((group1 == "Rnf12 ctrl" & group2 == "Rnf12 KD") |
           (group1 == "Zic3 ctrl" & group2 == "Zic3 KD"))

stat.test$y.position <- c(29, 29)


plot <- qc_samples_activators %>%
  ggplot( aes(x = KD, y = AmpR_normalized)) +
  geom_point(size=1.5,  aes(shape = replicate, color=KD)) +
  labs(title = "Spike-in") +
  ylab("Norm. reads")+
  xlab("")+
  scale_color_manual(values = c("Rnf12 ctrl" ="#3f414b",`Rnf12 KD`= "#f37748",
                                "Zic3 ctrl" ="#3f414b",`Zic3 KD`= "#f37748"))+
  ylim(c(0,32))+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)

ggsave(paste0("../Fig4/Figures/Spike-in_H3K9me3_activators.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

### repressors ####

castuner <- read.csv("fig_data/spike-in_castuner_repressors.txt", sep="\t")
krab <- read.csv("fig_data/spike-in_krab_repressors.txt", sep="\t")

qc_samples <- rbind(castuner, krab)

qc_samples %<>%
  separate(Sample, c("cell_line", "replicate", "antibody", "sample_number"), sep = "_")

qc_samples$KD <- if_else(qc_samples$cell_line == "NK06", "Rex1 KD",
                                    if_else(qc_samples$cell_line == "NK11", "Zfp281 KD", 
                                            if_else(qc_samples$cell_line == "LR15", "STC ctrl",
                                             if_else(qc_samples$cell_line == "TS35+", "Zfp280c ctrl",
                                             if_else(qc_samples$cell_line == "TS35-", "Zfp280c KD",
                                              if_else(qc_samples$cell_line == "TS38+", "Zfp36l1 ctrl", "Zfp36l1 KD"))))))


qc_samples$KD <- factor(qc_samples$KD, levels=c("Zfp280c ctrl", "Zfp280c KD", 
                                                "Zfp36l1 ctrl", "Zfp36l1 KD",
                                                "STC ctrl", "Rex1 KD", "Zfp281 KD"))

stat.test <- compare_means(
  AmpR_normalized ~ KD,
  data = qc_samples,
  method = "t.test", 
  paired = TRUE)

# format pvalue
stat.test$pval <- vapply(stat.test$p, format_p, character(1))

stat.test %<>%
  filter((group1 == "Zfp280c ctrl" & group2 == "Zfp280c KD") |
           (group1 == "Zfp36l1 ctrl" & group2 == "Zfp36l1 KD") |
           (group1 == "STC ctrl" & group2 == "Rex1 KD") |
           (group1 == "STC ctrl" & group2 == "Zfp281 KD"))

stat.test$y.position <- c(11, 12, 8, 10)


plot <- qc_samples %>%
  ggplot( aes(x = KD, y = log2(AmpR_normalized))) +
  geom_point(size=1.5,  aes(shape = replicate, color=KD)) +
  labs(title = "Spike-in") +
  ylab("Norm. reads")+
  xlab("")+
  scale_color_manual(values = c("Zfp280c ctrl" ="#3f414b",`Zfp280c KD`= "#3b9ad9",
                                "Zfp36l1 ctrl" ="#3f414b",`Zfp36l1 KD`= "#3b9ad9",
                                "NTC ctrl" ="#3f414b",`Rex1 KD`= "#3b9ad9",`Zfp281 KD`= "#3b9ad9"))+
  ylim(c(0,14))+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/Spike-in_H3K9me3_repressors.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

### RNF12-KO DKO ####

qc_samples <- read.csv("fig_data/spike_in_DKO.csv")

my_comparisons <- list(c("WT", "Rnf12 KO"), c("WT", "DKO"))

qc_samples$Sample <- factor(qc_samples$Sample, levels=c("WT", "Rnf12 KO", "DKO"))

qc_samples %<>%
  filter(Replicate!="Rep2") #to be able to do paired test

stat.test <- compare_means(
  AmpR_normalized ~ Sample,
  data = qc_samples,
  method = "t.test", 
  paired = TRUE)

# format pvalue
format_p <- function(p) {
  if (p < 0.001) {
    "P < 0.001"
  } else if (p < 0.01) {
    sprintf("P = %.3f", p)
  } else {
    sprintf("P = %.2f", p)
  }
}

stat.test$pval <- vapply(stat.test$p, format_p, character(1))

stat.test %<>%
  filter(group1 == "WT")

stat.test$y.position <- c(13, 15.5)



plot <- qc_samples %>%
  ggplot( aes(x = Sample, y = AmpR_normalized)) +
  geom_point(size=1.5,  aes(shape = Replicate, color=Sample)) +
  labs(title = "Spike-in") +
  ylab("Norm. reads")+
  xlab("")+
  scale_color_manual(values = c("WT" ="#3f414b",`Rnf12 KO`= "#f37748", DKO="#751f58"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  ylim(c(0,18))+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/Spike-in_H3K9me3_DKO.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
