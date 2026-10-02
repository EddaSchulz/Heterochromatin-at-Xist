## RNA-FISH
## Nelly Kanata ~ OWL Schulz
## Created: 07.07.2025
## Modified: 11.05.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)


setwd("./")


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

# DKO ####
# Load data

raw_data <- read.csv("./fig_data/RNF12-REX1-KOs_RNA-FISH.csv", sep=";")

# Calculate percentages
raw_data %<>%
  mutate(Xist_cloud_perc = Cloud/ Total*100,
         Xist_dots_perc = Dots/Total*100,
         Xist_positive = (Cloud+Dots+Biallelic)/Total *100,
         Xist_Biallelic_perc = Biallelic/Total*100,
         Xist_neg = 100-Xist_cloud_perc - Xist_dots_perc - Xist_Biallelic_perc)

raw_data$Cell_line <- ifelse(raw_data$Cell_line == "RNF12", "Rnf12 KO", raw_data$Cell_line)
raw_data$Cell_line <- factor(raw_data$Cell_line, levels = c("WT", "Rnf12 KO", "DKO"))


# stats
stat.test <- compare_means(
  Xist_positive ~ Cell_line,
  data = raw_data,
  method = "t.test"
)

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


stat.test %<>% filter(group1 == "WT")

stat.test$y.position <- c(62, 85)




xist_plot <- raw_data %>%
  ggplot(aes(x = Cell_line, y = Xist_positive)) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(width = 0.8), fill="#f47748") +
  stat_summary(aes(group = Cell_line, y = Xist_Biallelic_perc), 
               geom = "bar", fun = "mean", position = position_dodge(width = 0.8), color = NA, fill = "#e2360e") +
  geom_point(aes(x = Cell_line, y = Xist_positive),size = 1, position = position_dodge(width = 1), shape = 21) +
  scale_y_continuous(limits = c(0, 100), name = "Xist positive (%cells)") +
  scale_x_discrete(name = NULL)+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)
  

fix <- set_panel_size(xist_plot, height = unit(2.4, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/DKO_FISH.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

# Repressors ####
#### Day 4 ####
# Load data

aba_samples <- read.csv("./fig_data/D4_TFKDs_ABA_RNA-FISH.csv", sep=";")
dtag_samples <- read.csv("./fig_data/D4_TFKDs_dTAG_RNA-FISH.csv", sep=";")

raw_data_d4 <- rbind(aba_samples, dtag_samples)

#### Day 2 ####
# Load data

raw_data_d2 <- read.csv("./fig_data/D2_TFKDs_RNA-FISH.csv", sep=";")

raw_data <- rbind(raw_data_d2, raw_data_d4)

# Calculate percentages
raw_data %<>%
  mutate(Xist_cloud_perc = Cloud/ Total*100,
         Xist_dots_perc = Dots/Total*100,
         Xist_positive = (Cloud+Dots+Biallelic)/Total *100,
         Xist_Biallelic_perc = Biallelic/Total*100,
         Xist_neg = 100-Xist_cloud_perc - Xist_dots_perc - Xist_Biallelic_perc)

raw_data$KD <- if_else(raw_data$Cell_line == "NK06", "Rex1 KD", 
                       if_else(raw_data$Cell_line == "NK11", "Zfp281 KD", 
                               if_else(raw_data$Cell_line == "LR15", "NTC ctrl", 
                                       if_else(raw_data$Cell_line == "TS35+", "Zfp280c ctrl", 
                                               if_else(raw_data$Cell_line == "TS35-", "Zfp280c KD", 
                                                       if_else(raw_data$Cell_line == "TS38+", "Zfp36l1 ctrl", 
                                                               if_else(raw_data$Cell_line == "TS38-", "Zfp36l1 KD","" 
                                                               )))))))

raw_data$KD <- factor(raw_data$KD, levels = c("Zfp280c ctrl", "Zfp280c KD",
                                              "Zfp36l1 ctrl", "Zfp36l1 KD", 
                                              "NTC ctrl", "Rex1 KD", "Zfp281 KD"))

raw_data$Gene <- gsub(" .*","", raw_data$KD) 

# duplicate NTC to have as separate entry for each KD 
NTC <- raw_data[raw_data$Gene == "NTC" & raw_data$Day == "D4",]
NTC$Gene <- "Rex1" 

raw_data$Gene[raw_data$Gene == "NTC"] <- "Zfp281"

raw_data <- rbind(raw_data, NTC)
raw_data$Gene <- paste0(raw_data$Gene, " KD")
raw_data$Gene <- factor(raw_data$Gene, levels = c("Zfp280c KD", 
                                                  "Zfp36l1 KD", 
                                                  "Rex1 KD", "Zfp281 KD"))

raw_data$KD_vs_ctrl <- gsub(".* ", "", raw_data$KD)


# stats
stat.test <- compare_means(
  Xist_positive ~ KD_vs_ctrl,
  data = raw_data,
  method = "t.test",
  paired = TRUE,
  group.by = c("Day", "Gene")
)

stat.test$label <- vapply(stat.test$p, format_p, character(1))

stat.test$y.position <- c(95, 95, 90, 80, 85, 85)



xist_plot <- raw_data %>%
  ggplot(aes(x = KD_vs_ctrl, y = Xist_positive)) +
  facet_grid(rows=vars(Gene), cols=vars(Day), scales = "free")+
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(width = 0.6), 
               aes(fill=KD_vs_ctrl), width=0.7) +
  scale_fill_manual(values = c("ctrl" = "gray70", "KD" = "#3b9ad9")) +
  ggnewscale::new_scale_fill() +
  stat_summary(aes(group = KD, y = Xist_Biallelic_perc, fill = KD_vs_ctrl), 
               geom = "bar", fun = "mean", position = position_dodge(width = 0.6), 
               color = NA , width=0.7
  ) +
  scale_fill_manual(values = c("ctrl" ="gray20", "KD" = "#0e6fa5")) +
  geom_point(aes(x = KD_vs_ctrl, y = Xist_positive),size = 1, position = position_dodge(width = 1), 
             shape = 21, fill="white") +
  scale_y_continuous(limits = c(0, 100), name = "Cells (%)",
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(name = NULL)+
  stat_pvalue_manual(stat.test, label = "label", size = 6/2.8, vjust = -0.2)


fix <- set_panel_size(xist_plot, height = unit(1.8, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/TFKDs_FISH.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
