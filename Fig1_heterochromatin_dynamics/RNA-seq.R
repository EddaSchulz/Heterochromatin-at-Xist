#Script to plot screen Tfs expression dynamics
library(tidyverse)
library(egg)
library(gridExtra)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


setwd("./")

# load data

#Read deseq2 file and cpm
cpm_tx <- read.delim("./fig_data/GSE273071_CPM_RNA_timecourse.txt.gz")

tf_list <- c("Xist") # list of genes to plot


#Creates table with expression and only focuses on selected genes
cpm_long <- cpm_tx %>%  
  filter(gene %in% tf_list) %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("sex", "rep", "day"))


plot_df <- cpm_long %>%
  dplyr::rename(timepoint = "day") %>% 
  mutate(timepoint = as.numeric(str_remove(timepoint, "h")))

plot_df$timepoint <- plot_df$timepoint / 24


## Xist ####
plot <-  plot_df %>% 
  filter(gene == "Xist") %>% 
  ggplot(aes(x = timepoint, y = cpm, color = sex)) +
  facet_wrap(~ sex, scales = "free", nrow=1) +
    geom_point(size = 0.2) +
  geom_smooth(se = FALSE, linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), name = "Differentiation timepoint (days)") +
  ylab("RNA Expression (CPM)")+
  scale_color_manual(values = c("#676767", "#676767")) 

fix <- set_panel_size(plot, height = unit(1.8, "cm"), width = unit(1.8, "cm"))
grid.arrange(fix)
ggsave("./Figures/Timecourse_CPM_Xist-XO.pdf", fix, dpi = 300,
       useDingbats=FALSE)
