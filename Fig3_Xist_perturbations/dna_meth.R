## DNA methylation Plots
## Nelly Kanata ~ OWL Schulz
## Created: 26.10.2025
## Modified: 08.04.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


setwd("./")

# load data

xist_kd_XO <- read.csv("./fig_data/Xist_KD_XO_all_methylation_table.csv")


##### XO ######

filtered_df <- xist_kd_XO %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57")

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$condition <- factor(filtered_df$condition, levels = c("+dTAG", "-dTAG"))


### mean methylation plot
RE57_amplicons <- c("RE57_YY1", "RE57_CGI")

methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 5) %>% # coverage cutoff ## THIS IS NORMALLY 10 BUT THIS LIBRARY IS UNDERSEQUENCED
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, cell_line, condition,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) 
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))
methC$condition <- if_else(methC$condition == "+dTAG", "ctrl", "KD")

# stats
stat.test <- compare_means(
  mean_meth ~ condition,
  data = methC,
  method = "t.test", paired = TRUE
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
stat.test$y.position <- 1

plot <- ggplot(methC, aes( condition, mean_meth)) +
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=condition),
    alpha=0.6,
    width = 0.7
  ) +
  
  geom_point(aes(fill = condition), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values=c(`ctrl`="#3f414b", `KD`="#3b9ad9")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.2)) +  # Extends to 1.1 but labels stop at 1
  ylab("DNA methylation (%)")+
  xlab("")+
  ggtitle(paste0("")) +
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.1, "cm"))
grid.arrange(fix)

ggsave("./Figures/mean_meth_Xist_KD_XO.pdf", fix,
       dpi = 300, useDingbats=FALSE)


#### XX ####
xist_kd_XX <- read.csv("./fig_data/Xist_KD_XX_all_methylation_table.csv")

filtered_df <- xist_kd_XX %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$condition <- factor(filtered_df$condition, levels = c("+dTAG", "-dTAG"))


# mean
methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, cell_line, condition, KD,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(KD)
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))
methC$condition <- if_else(methC$condition == "+dTAG", "ctrl", "KD")
methC$KD <- paste0(methC$KD, " KD")

methC$mean_meth <- methC$mean_meth *100



# stats
stat.test <- compare_means(
  mean_meth ~condition,
  data = methC[ methC$KD=="Xist KD",],
  method = "t.test", paired=TRUE
)


stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(87)


plot <- methC %>% filter(KD=="Xist KD") %>%
  ggplot(aes( condition, mean_meth)) +
  facet_grid(cols=vars(KD))+
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=condition),
    alpha=0.6,
    width = 0.7
  ) +
  
  geom_point(aes(fill = condition), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values=c(`ctrl`="#3f414b", `KD`="#3b9ad9")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) + 
  ylab("Mean DNA methylation (%)")+
  xlab("")+
  ggtitle(paste0("")) +
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.1, "cm"))
grid.arrange(fix)

ggsave("./Figures/mean_meth_Xist_KD_XX.pdf", fix,
       dpi = 300, useDingbats=FALSE)


