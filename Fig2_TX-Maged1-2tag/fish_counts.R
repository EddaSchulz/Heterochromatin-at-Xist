## FISH counts
## Nelly Kanata ~ OWL Schulz
## Created: 27.04.2026
## Modified: 24.06.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


setwd("./")

### sorted ####
# Read data
fish_data <- readxl::read_xlsx("fig_data/RNA-FISH-cutntag-maged1-xist.xlsx")

fish_data %<>%
  separate(Sample, c("Cell_line", "timepoint", "replicate", "gate", "Xist"), sep="_")

fish_data %<>%
  select(1:6) %>%
  pivot_wider(names_from = Xist, values_from = Percentage)

fish_data$xist_positive=fish_data$MA+fish_data$BA


fish_data$Gate <- if_else(fish_data$gate == "p1", "Xist Off",
                          if_else(fish_data$gate == "p3", "Xist-Cast", 
                                  if_else(fish_data$gate == "p4", "Xist-B6", "unsorted (mix)")))
fish_data$Gate <- factor(fish_data$Gate, levels = c("unsorted (mix)", "Xist Off", "Xist-Cast","Xist-B6"))


xist_plot <- fish_data %>%
  ggplot(aes(x = Gate, y = xist_positive)) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(width = 0.8), fill="gray70") +
  stat_summary(aes(group = Gate, y = BA), 
               geom = "bar", fun = "mean", position = position_dodge(width = 0.8),  fill = "black") +
  geom_point(aes(x = Gate, y = xist_positive),size = 1, position = position_dodge(width = 1), shape = 21) +
  scale_y_continuous(limits = c(0, 100), name = "Cells (%)") +
  scale_x_discrete(name = "")

fix <- set_panel_size(xist_plot, height = unit(2.5, "cm"), width = unit(2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/sorting_FISH.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


### comparison with TX ####

fish_data <- readxl::read_xlsx("fig_data/fish_counts_validation5.xlsx")

fish_data %<>%
  separate(Sample, c("Cell_line", "timepoint", "replicate", "Xist"), sep="_")

fish_data %<>%
  select(1:5) %>%
  pivot_wider(names_from = Xist, values_from = Percentage)

fish_data$xist_positive=fish_data$MA+fish_data$BA

fish_data$Cell_line  <- if_else(fish_data$Cell_line == "Maged1", "TX-Maged1-2tag", "TX1072")

fish_data$Cell_line <- factor(fish_data$Cell_line, levels = c("TX1072", "TX-Maged1-2tag"))
fish_data$timepoint <- gsub("d", "", fish_data$timepoint)



stat.test <- compare_means(
  xist_positive ~ Cell_line,
  data = fish_data,
  group.by = c("timepoint"),
  method = "t.test", paired=FALSE
)



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
stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(98, 68)
stat.test$Cell_line <- "TX-Maged1-2tag"



xist_plot <- fish_data %>%
  ggplot(aes(x = timepoint, y = xist_positive, group=Cell_line)) +
  stat_summary( geom = "bar", fun = "mean", aes(width=0.7, fill=Cell_line),
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("TX1072" = "gray50",
                               "TX-Maged1-2tag" = "#F47748")) +
  
  ggnewscale::new_scale_fill() +
  stat_summary(aes(y = BA, width=0.7, fill=Cell_line), 
               geom = "bar", fun = "mean", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("TX1072" = "black",
                               "TX-Maged1-2tag" = "#e2360e"))+
  geom_point(aes(x = timepoint, y = xist_positive),size = 1, position = position_dodge(width = 1), shape = 21) +
  scale_y_continuous(limits = c(0, 100), name = "Cells (%)") +
  scale_x_discrete(name = "Differentiation timepoint (days)")+
  theme(axis.text.x = element_text(angle = 0,  hjust=0.5))+
  geom_text(data = stat.test,
            aes(label = pval, y = y.position), size=6/2.8,
            nudge_x = 0.05, color="black") 
  
  
fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/TX-Maged1-2tag_FISH.pdf"), fix,
       dpi = 300, useDingbats=FALSE)

