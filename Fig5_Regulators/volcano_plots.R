# diffbind volcano plots
# Nelly Kanata, OWL Schulz
# Created: 02.02.2026
# Modified: 16.02.2026


library(magrittr)
library(dplyr)
library(EnhancedVolcano)
library(egg)
library(GenomicRanges)

setwd("./")


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


#### RNF12 ####

# Read data
Rnf12KO <- read.csv(file = 
                      "./fig_data/deseq2/RNF12vsWT_deseq2.csv")

#sort on FC
Rnf12KO %<>%
  arrange(-Fold)

Rnf12KO_GR <- GRanges(
  seqnames = Rnf12KO[[1]],
  ranges   = IRanges(start = Rnf12KO[[2]], end = Rnf12KO[[3]]),
  strand   = "*",
  mcols    = Rnf12KO[, -(1:3)]
)

# where is RE57 chrX:103481047-103482684
query <- GRanges("chrX", IRanges(103481047, 103487599))
hits <- subsetByOverlaps(Rnf12KO_GR, query)
hits

best_hit <- hits[which.min(mcols(hits)$mcols.FDR)]

RE57 <-   if_else( length(best_hit) == 0, NA, 
                   paste0(seqnames(best_hit), ":", start(best_hit), "-\n", end(best_hit)))

regions_to_label <- c( RE57)

if (length(best_hit) == 0) {
  Rnf12KO$Label <- NA_character_
} else {
  Rnf12KO$Label <- if_else(
    !is.na(RE57) &
      as.vector(Rnf12KO$seqnames == seqnames(best_hit)) &
      Rnf12KO$start == start(best_hit),
    "XistP",
    NA_character_
  )
}
Rnf12KO$Label <- if_else(as.vector(Rnf12KO$seqnames == "chr7") &
                           Rnf12KO$start == 142576400, "H19", Rnf12KO$Label )

Rnf12KO$Label <- if_else(as.vector(Rnf12KO$seqnames == "chr6") &
                           Rnf12KO$start == 30735200, "Mest", Rnf12KO$Label )

n_peaks <- nrow(Rnf12KO)

plot <- print(EnhancedVolcano(Rnf12KO,
                              lab = Rnf12KO$Label,
                              x = 'Fold',
                              y = 'FDR',
                              #subtitle=saveprefix,
                              pointSize = 0.5,
                              colAlpha=0.8,
                              
                              pCutoff = 0.05,
                              FCcutoff = 0.5,
                              #ylim = c(-0, 2.5),
                              xlim = c(-4, 4),
                              col=c('#3f414b', '#3f414b', '#3f414b', '#f37748'),
                              selectLab = c("XistP", "H19", "Mest"),
                              boxedLabels = FALSE,
                              labSize = 6.0/ 2.845, # because it's not the same scale as element_text
                              drawConnectors = TRUE,
                              max.overlaps = Inf,
                              xlab = bquote(~Log[2]~ " H3K9me3 fold change KD/WT"),
                              ylab = bquote(~-Log[10]~italic("FDR")),
                              legendPosition = 'right',
                              title="Rnf12 KO vs WT",
                              titleLabSize=8/2.8,
                              captionLabSize = 8/2.8,
                              subtitle = "",
                              caption = paste0(n_peaks, " peaks")))+
  annotate("text", x = -2.5, y = 32, 
           label = paste0("Down: ", sum(Rnf12KO$direction == "Down") ), size=6/ 2.845)+
  annotate("text", x = 2, y = 32, 
           label = paste0("Up: ", sum(Rnf12KO$direction == "Up") ), size=6/ 2.845)+
  theme_classic() + 
  theme(legend.text = element_text(size = 6), axis.text = element_text(size = 6),
        axis.title = element_text(size = 6), strip.text = element_text(size = 6),
        strip.background = element_blank(), legend.title = element_blank(), plot.title = element_text(size=8),
        plot.caption = element_text(size = 6))

fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)


ggsave(paste0("Figures/RNF12KOvsWT_VolcanoPlot.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



#### DKO ####

# Read data
DKO <- read.csv(file = 
                      "./fig_data/deseq2/DKOvsWT_deseq2.csv")

#sort on FC
DKO %<>%
  arrange(-Fold)

DKO_GR <- GRanges(
  seqnames = DKO[[1]],
  ranges   = IRanges(start = DKO[[2]], end = DKO[[3]]),
  strand   = "*",
  mcols    = DKO[, -(1:3)]
)

# where is RE57 chrX:103481047-103482684
query <- GRanges("chrX", IRanges(103481047, 103487599))
hits <- subsetByOverlaps(DKO_GR, query)
hits

best_hit <- hits[which.min(mcols(hits)$mcols.FDR)]

RE57 <-   if_else( length(best_hit) == 0, NA, 
                   paste0(seqnames(best_hit), ":", start(best_hit), "-\n", end(best_hit)))

regions_to_label <- c( RE57)

if (length(best_hit) == 0) {
  DKO$Label <- NA_character_
} else {
  DKO$Label <- if_else(
    !is.na(RE57) &
      as.vector(DKO$seqnames == seqnames(best_hit)) &
      DKO$start == start(best_hit),
    "XistP",
    NA_character_
  )
}

DKO$Label <- if_else(as.vector(DKO$seqnames == "chr7") &
                       DKO$start == 143259200, "Kcnq1ot1", DKO$Label )

n_peaks <- nrow(DKO)


plot <- print(EnhancedVolcano(DKO,
                              lab = DKO$Label,
                              x = 'Fold',
                              y = 'FDR',
                              pointSize = 0.5,
                              colAlpha=0.8,
                              
                              pCutoff = 0.05,
                              FCcutoff = 0.5,
                              xlim = c(-4, 4),
                              col=c('#3f414b', '#3f414b', '#3f414b', '#f37748'),
                              selectLab = c("XistP", "Kcnq1ot1"),
                              boxedLabels = FALSE,
                              labSize = 6.0/ 2.845, # because it's not the same scale as element_text
                              drawConnectors = TRUE,
                              max.overlaps = Inf,
                              xlab = bquote(~Log[2]~ " H3K9me3 fold change KD/WT"),
                              ylab = bquote(~-Log[10]~italic("FDR")),
                              legendPosition = 'right',
                              title="DKO vs WT",
                              titleLabSize=8/2.8,
                              captionLabSize = 8/2.8,
                              subtitle = "",
                              caption = paste0(n_peaks, " peaks")))+
  annotate("text", x = -2.5, y = 28, 
           label = paste0("Down: ", sum(DKO$direction == "Down") ), size=6/ 2.845)+
  annotate("text", x = 2, y = 28, 
           label = paste0("Up: ", sum(DKO$direction == "Up") ), size=6/ 2.845)+
  theme_classic() + 
  theme(legend.text = element_text(size = 6), axis.text = element_text(size = 6),
        axis.title = element_text(size = 6), strip.text = element_text(size = 6),
        strip.background = element_blank(), legend.title = element_blank(), plot.title = element_text(size=8),
        plot.caption = element_text(size = 6))

fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(2.5, "cm"))

grid.arrange(fix)

ggsave(paste0("Figures/DKOvsWT_VolcanoPlot.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
