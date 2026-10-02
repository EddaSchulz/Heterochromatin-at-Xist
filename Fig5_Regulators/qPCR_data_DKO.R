## qPCR Analysis
## Nelly Kanata ~ OWL Schulz
## Created: 26.10.2025
## Modified: 26.10.2025

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


# load data

rel.data_long <- read.csv("./fig_data/qPCR_Rnf12-DKO.csv")

# pre-processing
rel.data_long %<>%
  select(-Sample)


rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix","Dnmt3b", "Fgf5",  "Otx2", "Setdb1",
                                                             "Zfp280c", "Zfp36l1", "Zfp281",  "Rex1",   "Myc",  "Nanog" ,
                                                             "Esrrb",   "Prdm14" , "Oct4" ))

rel.data_long$Cell_line <- factor(rel.data_long$Cell_line, levels=c("WT", "RNF12", "DKO"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix"))
  
filtered_rel.data_long$KD <- if_else(filtered_rel.data_long$Cell_line == "WT", "WT", "KD")
  
  # calculate significance
  test_df <- expand.grid(Cell_line = unique(filtered_rel.data_long$Cell_line), 
                         gene = unique(filtered_rel.data_long$gene),
                         Day = unique(filtered_rel.data_long$Day))

  
test_df$pval <- mapply(function(cl, g, d) {
  t.test(
    filtered_rel.data_long$log2_re[
      filtered_rel.data_long$Cell_line == cl &
        filtered_rel.data_long$gene == g &
        filtered_rel.data_long$Day == d
    ],
    filtered_rel.data_long$log2_re[
      filtered_rel.data_long$Cell_line == "WT" &
        filtered_rel.data_long$gene == g &
        filtered_rel.data_long$Day == d
    ]
  )$p.value
},
test_df$Cell_line, test_df$gene, test_df$Day
)

test_df %<>%
  filter(Cell_line !="WT")



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

test_df$pval <- vapply(test_df$pval, format_p, character(1))



plot <- filtered_rel.data_long %>%
  ggplot(aes(x=Day, y=log2_re, group=Cell_line)) + 
  facet_grid(rows=vars(gene), scales = "free") +
  geom_point(aes(color = Cell_line), alpha=1, position = position_dodge(width=0.75), shape=16) +
  xlab("Differentiation timepoint (days)") + ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("WT"= "lightgray", "RNF12"= "#f37748", "DKO" = "#800080"), name = '') +
  geom_text(data = test_df, aes(label = pval, 
                                y = c(-7, -7, -4.3, -4.3, -6,-1, -4.4, -4.4, -6, 0, -4.5, -4.5)), 
            size=6/2.8, angle=90, hjust=-0.1, vjust=1.2,
            position = position_dodge(width = 0.5), show.legend = FALSE)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 1))) +
    theme(legend.position = "bottom")


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_Tsix_Rnf12-DKO.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



#### effect on differentiation #####

filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Nanog", "Esrrb", "Fgf5", "Otx2")) %>%
  mutate(log2FC_NTC = log2(re/re_NTC)) 

filtered_rel.data_long$Day <- factor(filtered_rel.data_long$Day, levels=c("D0", "D2", "D4"))

# calculate significance
test_df <- expand.grid(Cell_line = unique(filtered_rel.data_long$Cell_line), 
                       gene = unique(filtered_rel.data_long$gene),
                       Day = unique(filtered_rel.data_long$Day))


test_df$pval <- mapply(function(cl, g, d) {
  t.test(
    filtered_rel.data_long$log2_re[
      filtered_rel.data_long$Cell_line == cl &
        filtered_rel.data_long$gene == g &
        filtered_rel.data_long$Day == d
    ],
    filtered_rel.data_long$log2_re[
      filtered_rel.data_long$Cell_line == "WT" &
        filtered_rel.data_long$gene == g &
        filtered_rel.data_long$Day == d
    ]
  )$p.value
},
test_df$Cell_line, test_df$gene, test_df$Day
)
test_df %<>%
  filter(Cell_line !="WT")

test_df$pval <- vapply(test_df$pval, format_p, character(1))



plot <- filtered_rel.data_long %>%
  ggplot(aes(x=Day, y=log2_re, group =Cell_line)) + 
  facet_grid(rows=vars(gene), scales = "free_y") +
  geom_point(aes(color = Cell_line), alpha=1, position = position_dodge(width=0.75), shape=16) +
  xlab("Differentiation timepoint (days)") + ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("WT"= "lightgray", "RNF12"= "#f37748", "DKO" = "#800080"), name = '',) +
  geom_text(data = test_df, aes(label = pval, y = c(-5, -5, 1, 1, -1,-1, -2, -2, -1, -1, 0, 0)), size=6/2.8, angle=90, hjust=-0.1, vjust=1.2,
            position = position_dodge(width = 0.5), show.legend = FALSE)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.8))) +
  theme(legend.position = "bottom")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_KO_differentiation.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

