## DNA methylation Plots
## Nelly Kanata ~ OWL Schulz
## Created: 26.10.2025
## Modified: 26.06.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  #axis.text.x = element_text(angle = 45, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


setwd("./")

# load data

timecourse <- read.csv("./fig_data/timecourse_all_methylation_table.csv")

filtered_df <- timecourse %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$timepoint <- factor(filtered_df$timepoint)
filtered_df$timepoint <- as.numeric(gsub("T", "", filtered_df$timepoint))


RE57_amplicons <- c("RE57_YY1", "RE57_CGI")

methC <- filtered_df %>%
 filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number, amplicon) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, cell_line, amplicon,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(cell_line)
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))

methC$mean_meth <- methC$mean_meth *100

methC$Day <- methC$timepoint /24


### RNA-Seq from the same experiment ####

CPM <- read.delim("./fig_data/GSE273071_CPM_RNA_timecourse.txt.gz")

xist_rnaseq <- CPM %>%
  filter(gene =="Xist") %>%
  pivot_longer(-c(1:2), names_to = "sample", values_to = "CPM") %>%
  filter(sample!="XO_R1_36h") %>% 
  separate(sample, c("cell_line", "replicate", "timepoint"), sep = "_") %>%
  mutate(timepoint = as.numeric(str_remove(timepoint, "h")))


xist_rnaseq %<>% 
  group_by(gene, timepoint, cell_line) %>% 
  na.omit() %>% 
  mutate(m=mean(CPM), s=sd(CPM))


#scale range
scale_min <- 0
scale_max <- 100

# min-max normalization
xist_rnaseq$CPM_scaled <- scale_min + (xist_rnaseq$CPM-min(xist_rnaseq$CPM)) * (scale_max - scale_min)/
  (max(xist_rnaseq$CPM)-min(xist_rnaseq$CPM))


xist_rnaseq$Day <- xist_rnaseq$timepoint / 24

# plot

plot <- methC %>%
  # filter(cell_line == "XX") %>%
  ggplot( aes(Day, mean_meth, group = cell_line)) +
  facet_wrap(~ cell_line, axes = "all")+
  #barplot
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "#84bcda",
    width = 0.3
  ) +

  # Mean points (filled vs open)
  stat_summary(fun = mean,
               aes(group = replicate),
               geom = "point",
               color = "black", shape=23,fill="#3b9ad9",
               size = 1,) +
  
  scale_y_continuous( limits = c(0,100),breaks = seq(0, 100, by = 25),
                      sec.axis = 
                        sec_axis(~ min(xist_rnaseq$CPM) + (. - 0) * (max(xist_rnaseq$CPM) - min(xist_rnaseq$CPM)) / (scale_max - 0), 
                                          name = "Xist expression (CPM)")) +
  scale_x_continuous(
    breaks = seq(0, 4, by = 1)
  )+
  ylab("DNA methylation (%)") +
  xlab("Differentiation timepoint (days)") +
  
  # rnaseq data
  # Mean line
  stat_summary(data = xist_rnaseq,
               aes(x = Day, y = CPM_scaled, group = cell_line),
               fun = mean,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.5,
               color = "black") +
  
  # Mean points
  geom_point(data = xist_rnaseq,
               aes(x = Day, y = CPM_scaled, group = cell_line),
               shape=21,
               size = 1,
               color = "black", fill="white") 


fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)

ggsave("./Figures/RE57_mean_meth_RNA_timecourse.pdf", fix,
       dpi = 300, useDingbats=FALSE)

