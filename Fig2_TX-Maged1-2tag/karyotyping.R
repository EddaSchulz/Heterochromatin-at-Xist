## qPCR Analysis
## Nelly Kanata ~ OWL Schulz
## Script from Melissa Bothe
## Created: 30.04.2026
## Modified: 30.04.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

setwd("./")

theme_set(theme_minimal() + 
            theme(legend.text = element_text(size = 6), 
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_text(size=6)))


# Define basic functions
div_fun <- function(a,b) {
  log2(a / b)
}

sum_fun <- function(a) {
  a / sum(a)
}


# Read-in counts
counts <- read.delim("./fig_data/karyotyping_reads.txt")
control <- "SC02_A3_con"

# Normalize to control
karyo_df <- counts %>%
  select(-start, -end) %>% 
  mutate_at(vars(-chr), ~sum_fun(.)) %>%
  mutate_at(vars(-chr, -!!sym(control)), ~div_fun(., !!sym(control)))  %>%
  select(-!!sym(control)) %>%
  pivot_longer(-chr, names_to = "line", values_to = "normed")

# Order chromosomes
karyo_df$chr <- factor(karyo_df$chr, levels = unique(karyo_df$chr))
karyo_df$line <- gsub("SC106_", "", karyo_df$line)

# Plot the heatmap
heatmap_plot <- ggplot(karyo_df, aes(x = chr, y = line, fill = normed)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  labs(fill = "log2(Clone / XX control)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  theme(legend.key.width=unit(0.5,"cm"),
        legend.key.height =unit(0.3,"cm"))



fix <- set_panel_size(heatmap_plot, height = unit(1, "cm"), width = unit(8, "cm"))
grid.arrange(fix)


ggsave("./Figures/karyotyping_heatmap_plot.pdf", fix, dpi = 300)
