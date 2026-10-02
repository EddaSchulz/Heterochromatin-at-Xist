## New peaks
## Nelly Kanata ~ OWL Schulz
## Created: 07.02.2026
## Modified: 30.06.2026

library(GenomicRanges)
library(ggplot2)
library(egg)
library(magrittr)
library(rtracklayer)
library(ggpubr)
library(dplyr)

setwd("./")


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_text(size=6),
                  title = element_text(size=8)))



DKO <- read.csv(file = 
                  "./fig_data/deseq2/DKOvsWT_deseq2.csv")

##### How many of the called peaks are new upon KO? (not called in WT) #####

DKO$new_peak <- DKO$`Called2` == 0
DKO$new_peak <- factor(DKO$new_peak, levels = c("TRUE", "FALSE"))
DKO$new_peak <- if_else(DKO$new_peak == "TRUE", "New", "Old")
plot <- DKO %>%
  filter(`direction`=="Up") %>%
  ggstatsplot::ggbarstats(
    x = new_peak, y = `direction`,
    results.subtitle = FALSE,
    label = "counts",
    title = "",
    legend.title = "",
    label.args = list(size = 2, fill = "white"),
    sample.size.label.args =  c("", "", "", "")) +
  scale_fill_manual(values = c("New"= "#f37748", "Old" = "white"))+
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.11))) +
  ylab("% of peaks") +
  xlab("Peaks gaining\nH3K9me3 upon DKO")+
  theme(
    legend.key.size = unit(0.4, 'cm'),
    text = element_text(size = 6),
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )

fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/new_peaks_DKO.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

