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
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

# format pvalue function
format_p <- function(p) {
  if (p < 0.001) {
    "P < 0.001"
  } else if (p < 0.01) {
    sprintf("P = %.3f", p)
  } else {
    sprintf("P = %.2f", p)
  }
}


##### Xist KD XX #####
Xist_XX <- read.csv("./fig_data/qPCR_XistKD-XX.csv")
Xist_XX %<>%
  select(-Sample) 

rel.data_long <- Xist_XX

rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix_spliced", "Nanog" ,
                                                          "Esrrb",   "Prdm14","Rnf12"  ))
rel.data_long$Cell_line <- paste0(rel.data_long$Gene, " KD")
rel.data_long$Cell_line <- factor(rel.data_long$Cell_line, 
                                  levels=c("Tsix KD", "Xist KD"))
rel.data_long$KD <- if_else(rel.data_long$KD == "WT", "ctrl", "KD")
rel.data_long$KD <- factor(rel.data_long$KD, 
                           levels=c("ctrl", "KD"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix_spliced")) 


# calculate significance
test_df <- expand.grid(Cell_line = unique(filtered_rel.data_long$Cell_line), 
                       gene = unique(filtered_rel.data_long$gene))

test_df$pval <- mapply(function(a,b)
{t.test(filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                 filtered_rel.data_long$gene==b & 
                                 filtered_rel.data_long$KD=="ctrl",]$log2_re,
        filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                 filtered_rel.data_long$gene==b & 
                                 filtered_rel.data_long$KD=="KD",]$log2_re)$p.value},
test_df$Cell_line, test_df$gene)

test_df$pval <- vapply(test_df$pval, format_p, character(1))


test_df$KD <- "KD"

filtered_rel.data_long$gene <- factor(filtered_rel.data_long$gene, levels=c("Nanog", "Esrrb", "Prdm14", "Rnf12", "Xist", "Tsix_spliced"))

plot <- filtered_rel.data_long %>%
  filter(Cell_line=="Xist KD") %>%
  ggplot(aes(x=KD, y=log2_re, color = KD)) + 
  facet_wrap(~gene, scales = "free") +
  geom_point(alpha=0.8, position = position_dodge(width=0.75), shape=16) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("ctrl"= "#3f414b", "KD"= "#3b9ad9"), name = '',) +
  geom_text(data = test_df[test_df$Cell_line=="Xist KD",], aes(label = pval, y = c(0.5,-5)), 
            size=6/2.8,
            color = "black",
            nudge_x = -0.5)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme(legend.position = "top")


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_KD-XX.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)


# x-fold KD?

filtered_rel.data_long %>%
  filter(Cell_line=="Xist KD") %>% 
  filter(gene == "Xist") %>% 
  group_by(KD) %>%
  summarize(mean_re=mean(re), .groups = "drop") %>%
  mutate(fold_KD = mean_re[KD == "ctrl"] / mean_re,
         percent_KD = (mean_re[KD == "ctrl"]-mean_re)/ mean_re[KD == "ctrl"]*100)


##### Xist KD XO #####
# load data
Xist_KD <- read.csv("./fig_data/qPCR_XistKD-XO.csv")


# pre-processing

Xist_KD %<>%
  select(-Sample) 


rel.data_long <- Xist_KD


rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist",
                                                          "Tsix", "Fgf5", "Otx2", "Nanog" ,
                                                             "Esrrb"  ))

rel.data_long$KD <- if_else(rel.data_long$KD == "WT", "ctrl", "KD")
rel.data_long$KD <- factor(rel.data_long$KD, 
                           levels=c("ctrl", "KD"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist","Tsix")) 


# calculate significance
test_df <- expand.grid(Cell_line = unique(filtered_rel.data_long$Cell_line), 
                       gene = unique(filtered_rel.data_long$gene))

test_df$pval <- mapply(function(a,b)
{t.test(filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                 filtered_rel.data_long$gene==b & 
                                 filtered_rel.data_long$KD=="ctrl",]$log2_re,
        filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                 filtered_rel.data_long$gene==b & 
                                 filtered_rel.data_long$KD=="KD",]$log2_re)$p.value},
test_df$Cell_line, test_df$gene)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$KD <- "KD"

plot <- filtered_rel.data_long %>%
  ggplot(aes(x=KD, y=log2_re, color = KD)) + 
  facet_wrap(~gene, scales = "free") +
  geom_point(alpha=0.8, position = position_dodge(width=0.75), shape=16) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("ctrl"= "#3f414b", "KD"= "#3b9ad9"), name = '',) +
  geom_text(data = test_df, aes(label = pval, y = c(-4, -6.5)), 
            size=6/2.8,
            color = "black",
            nudge_x = -0.5)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme(legend.position = "top")


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_KD_XO.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)

# x-fold KD

filtered_rel.data_long %>%
  filter(gene == "Xist") %>% 
  group_by(KD) %>%
  summarize(mean_re=mean(re), .groups = "drop") %>%
  mutate(fold_KD = mean_re[KD == "ctrl"] / mean_re,
         percent_KD = (mean_re[KD == "ctrl"]-mean_re)/ mean_re[KD == "ctrl"]*100)



##### Transient Dox experiment (Supplementary Figure 3c) #####

# load data
Dox_data <- read.csv("./fig_data/qPCR_Xist_Dox.csv")


Dox_data %<>%
  select(-Sample) 


rel.data_long <- Dox_data


rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist"))
rel.data_long$Dox <- gsub("\\+", "", rel.data_long$Dox)
rel.data_long$Dox <- factor(rel.data_long$Dox, levels=c("no Dox", "24h Dox", "48h Dox"))
rel.data_long$timepoint <- if_else(rel.data_long$timepoint == 24, "Day 1", "Day 4")

# Xist relative expression
filtered_rel.data_long <- rel.data_long

filtered_rel.data_long$Dox <- if_else(filtered_rel.data_long$Dox == "24h Dox", "1d Dox", 
                                      if_else(filtered_rel.data_long$Dox == "48h Dox", "2d Dox", 
                                              filtered_rel.data_long$Dox))
filtered_rel.data_long$Dox <- factor(filtered_rel.data_long$Dox, levels=c("no Dox", "1d Dox", "2d Dox"))

# stats
stat.test <- compare_means(
  log2_re ~ Dox,
  data = filtered_rel.data_long,
  group.by = c("timepoint"),
  method = "t.test", paired=FALSE
)
stat.test %<>%
  filter(group1 == "no Dox")

stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(-2, -1, -2, -1)



plot <- filtered_rel.data_long %>%
  ggplot(aes(x=Dox, y=log2_re, color = Dox)) + 
  facet_grid(cols=vars(timepoint), scales = "fixed") +
  geom_point(alpha=0.8, position = position_dodge(width=0.75), shape=16) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Dox), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("no Dox"= "black", "1d Dox"= "#e2360e", "2d Dox" = "#F47748"), name = '',) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),)


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Transient_Dox_Xist.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)

