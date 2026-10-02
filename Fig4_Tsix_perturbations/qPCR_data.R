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

# load data
TsixKD <- read.csv("./fig_data/qPCR_TsixKD.csv")
Xist_OE <- read.csv("./fig_data/qPCR_Xist_OE.csv")
Xist_OE_Tsix_KD <- read.csv("./fig_data/qPCR_XOE-TsixKD.csv")

# pre-processing
TsixKD %<>%
  select(-Sample) 

Xist_OE %<>%
  select(-Sample) 


##### Xist OE #####
rel.data_long <- Xist_OE


rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix", "Nanog" ,
                                                             "Esrrb",   "Prdm14","Rnf12"  ))
rel.data_long$Cell_line <- if_else(rel.data_long$Cell_line == "dXIC", "XXΔXIC", "XX")
rel.data_long$Cell_line <- factor(rel.data_long$Cell_line, 
                                  levels=c("XXΔXIC", "XX"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix"))
  
filtered_rel.data_long$condition <- if_else(filtered_rel.data_long$Dox == "-Dox", "ctrl", "OE")
filtered_rel.data_long  %<>%
  filter(Cell_line == "XXΔXIC")
  
# calculate significance
  test_df <- expand.grid(Cell_line = unique(filtered_rel.data_long$Cell_line), 
                         Day = unique(filtered_rel.data_long$Day))
  
  test_df$pval <- mapply(function(a,b)
  {t.test(filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                   filtered_rel.data_long$Day==b & 
                                   filtered_rel.data_long$Dox=="-Dox",]$log2_re,
          filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                   filtered_rel.data_long$Day==b & 
                                   filtered_rel.data_long$Dox=="+Dox",]$log2_re
          )$p.value},
  test_df$Cell_line, test_df$Day)
  

  test_df$pval <- vapply(test_df$pval, format_p, character(1))
  
test_df$condition <- "OE"


plot <- filtered_rel.data_long %>%
  ggplot(aes(x=condition, y=log2_re, color = Dox)) + 
  facet_wrap(~Day, scales = "free") +
  geom_point(alpha=0.8, position = position_dodge(width=0.75), shape=16) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("-Dox"= "#3f414b", "+Dox"= "#f37748"), name = '') +
  geom_text(data = test_df, aes(label = pval, y =c(-0.5, 0.5)), 
            size=6/2.8,
            color = "black",
            nudge_x = -0.5)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.20))) +
  theme(legend.position = "top")


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_Xist-OE.pdf"), fix,
       dpi = 300, 
      device = cairo_pdf)


##### Tsix #####
rel.data_long <- TsixKD

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
                                 # filtered_rel.data_long$Day==c & 
                                 filtered_rel.data_long$KD=="ctrl",]$log2_re,
        filtered_rel.data_long[filtered_rel.data_long$Cell_line==a & 
                                 filtered_rel.data_long$gene==b & 
                                 # filtered_rel.data_long$Day==c & 
                                 filtered_rel.data_long$KD=="KD",]$log2_re)$p.value},
test_df$Cell_line, test_df$gene)

test_df$pval <- vapply(test_df$pval, format_p, character(1))


test_df$KD <- "KD"

# Tsix KD only
filtered_rel.data_long$gene <- factor(filtered_rel.data_long$gene, levels=c("Nanog", "Esrrb", "Prdm14", "Rnf12", "Tsix_spliced", "Xist"))

plot <- filtered_rel.data_long %>%
  filter(Cell_line=="Tsix KD") %>%
  ggplot(aes(x=KD, y=log2_re, color = KD)) + 
  facet_wrap(~gene, scales = "free") +
  geom_point(alpha=0.8, position = position_dodge(width=0.75), shape=16) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("ctrl"= "#3f414b", "KD"= "#3b9ad9"), name = '',) +
  geom_text(data = test_df[test_df$Cell_line=="Tsix KD",], aes(label = pval, y = c(0.5,-5)), 
            size=6/2.8,
            color = "black",
            nudge_x = -0.5)+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme(legend.position = "top")


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Tsix_KD.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)

# x-fold KD?

filtered_rel.data_long %>%
  filter(Cell_line=="Tsix KD") %>% 
  filter(gene == "Tsix_spliced") %>% 
  group_by(KD) %>%
  summarize(mean_re=mean(re), .groups = "drop") %>%
  mutate(fold_KD = mean_re[KD == "ctrl"] / mean_re,
         percent_KD = (mean_re[KD == "ctrl"]-mean_re)/ mean_re[KD == "ctrl"]*100)



##### Xist OE Tsix KD #####

rel.data_long <- Xist_OE_Tsix_KD

rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix", "Tsix_EN73_EN74",
                                                          "Tsix_ES371_ES372","Nanog" ,
                                                          "Esrrb", "Fgf5", "Otx2"))

rel.data_long$KD <- if_else(rel.data_long$dTAG == "+dTAG" & rel.data_long$Dox == "-Dox", "ctrl",
                    if_else(rel.data_long$dTAG == "-dTAG" & rel.data_long$Dox == "+Dox", "Xist OE Tsix KD",
                    if_else(rel.data_long$dTAG == "+dTAG" & rel.data_long$Dox == "+Dox", "Xist OE",
                            "Tsix KD") ))
rel.data_long$KD <- factor(rel.data_long$KD, 
                           levels=c("ctrl", "Tsix KD", "Xist OE", "Xist OE Tsix KD"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix_EN73_EN74")) # spliced Tsix


# calculate significance
stat.test_xist <- compare_means(
  log2_re ~ KD,
  data = filtered_rel.data_long[filtered_rel.data_long$gene == "Xist",],
  method = "t.test", paired=FALSE
)
stat.test_xist %<>%
  filter(group1 == "ctrl")

stat.test_xist$pval <- vapply(stat.test_xist$p, format_p, character(1))

stat.test_xist$y.position <- c(-5,-1 ,1 )


plot <- filtered_rel.data_long %>%
  filter(gene=="Xist") %>%
  ggplot(aes(x=KD, y=log2_re)) + 
  facet_grid(cols=vars(gene), rows=vars(Cell_line), scales = "free") +
  geom_point(aes(alpha=dTAG, shape=Replicate, color = KD), position = position_dodge(width=0.5)) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, 
               lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("ctrl" ="#3f414b","Tsix KD" = "#3f414b", 
                                "Xist OE"= "#f37748", "Xist OE Tsix KD"= "#f37748" ))+
  scale_alpha_manual(values=c(0.6, 1)) +

  stat_pvalue_manual(stat.test_xist, label = "pval", size = 6/2.8, vjust = -0.2)+
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_OE-KD_spliced.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)



# stats
stat.test_tsix <- compare_means(
  log2_re ~ KD,
  data = filtered_rel.data_long[filtered_rel.data_long$gene == "Tsix_EN73_EN74",],
  method = "t.test", paired=FALSE
)
stat.test_tsix %<>%
  filter(group1 == "ctrl")

stat.test_tsix$pval <- vapply(stat.test_tsix$p, format_p, character(1))

stat.test_tsix$y.position <- c(-4.5,-3.75 ,-3 )


plot <- filtered_rel.data_long %>%
  filter(gene=="Tsix_EN73_EN74") %>%
  ggplot(aes(x=KD, y=log2_re)) + 
  facet_grid(cols=vars(gene), rows=vars(Cell_line), scales = "free") +
  geom_point(aes(alpha=dTAG, shape=Replicate, color = KD), position = position_dodge(width=0.5)) +
  xlab("") + 
  ylab(expression("Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = Cell_line), geom = "crossbar", fun = "mean", width = 0.5, 
               lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("ctrl" ="#3f414b","Tsix KD" = "#3f414b", 
                                "Xist OE"= "#f37748", "Xist OE Tsix KD"= "#f37748" ))+
  scale_alpha_manual(values=c(0.6, 1)) +
  stat_pvalue_manual(stat.test_tsix, label = "pval", size = 6/2.8, vjust = -0.2)+
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_OE-KD_spliced_tsix.pdf"), fix,
       dpi = 300, 
       device = cairo_pdf)



filtered_rel.data_long %>%
  filter(gene == "Tsix_EN73_EN74") %>% 
  group_by(KD) %>%
  summarize(mean_re=mean(re), .groups = "drop") %>%
  mutate(fold_KD = mean_re[KD == "Xist OE"] / mean_re,
         percent_KD = (mean_re[KD == "Xist OE"]-mean_re)/ mean_re[KD == "Xist OE"]*100)

