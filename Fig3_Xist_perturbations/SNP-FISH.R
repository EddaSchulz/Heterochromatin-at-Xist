# SNP FISH
# Nelly Kanata, OWL Schulz
# Created: 28.05.2026

library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(magrittr)
library(ggpubr)


# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

setwd("./")


# Read data
allelic_data <- read.csv("fig_data/transient_Dox_SNP-FISH.csv", sep = ";")

allelic_data$Day <- if_else(allelic_data$Day == "24h", "Day 1",
                                            if_else(allelic_data$Day == "96h", "Day 4", NA))

allelic_data$Condition <- gsub("\"", "", allelic_data$Condition)
allelic_data$Dox  <- if_else(allelic_data$Condition == "-Dox", "no Dox", 
                             "2d Dox")
allelic_data$Dox <- factor(allelic_data$Dox, levels = c("no Dox", "2d Dox"))

allelic_data$Cast <- allelic_data$cy5 / allelic_data$Total
allelic_data$B6 <- allelic_data$white / allelic_data$Total
allelic_data$Biallelic_perc <- allelic_data$Biallelic / allelic_data$Total


allelic_data %<>%
  pivot_longer(cols = c(Cast, B6, Biallelic_perc),
               names_to = "allele",
               values_to = "fraction")

allelic_data$allele <- factor(allelic_data$allele, levels = c("Cast", "B6", "Biallelic_perc"))

allelic_data$fraction <- allelic_data$fraction * 100

# remove NA rows
allelic_data %<>%
  filter(!is.na(Total))

# for errorbars
summary_df <- allelic_data %>%
  group_by(Day, Dox, allele) %>%
  summarise(
    mean = mean(fraction),
    se = sd(fraction) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(Day, Dox, desc(allele)) %>%
  group_by(Day, Dox) %>%
  mutate(
    ymax = cumsum(mean),
    ymin = ymax - mean,
    error_center = ymax
  )


# Plot
my_comparisons <- list(
  c("no Dox", "48h Dox")
)

plot <- ggplot(summary_df,
       aes(x = Dox, y = mean,
           fill = interaction(allele, Dox))) +
  
  facet_grid(cols = vars(Day)) +
  
  geom_col(
    position = "stack",
    width = 0.8,
    alpha = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = error_center - se,
      ymax = error_center + se
    ),
    width = 0.05
  )+
  scale_y_continuous(limits =c(0,100),     breaks = seq(0, 100, 25)) +
  
  scale_fill_manual(values = c(
    "B6.no Dox" = "gray20",      
    "Cast.no Dox" = "gray",      
    "Biallelic_perc.no Dox" = "black",
    
    
    "B6.2d Dox" = "#F47748",   
    "Cast.2d Dox" = "gray",   
    "Biallelic_perc.2d Dox" = "red"
  )) +
  labs(x = "",
       y ="Cells (%)") +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.border = element_blank()
  )

fix <- set_panel_size(plot, height = unit(2.2, "cm"), width = unit(1.2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/SNP-FISH_byDay.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


