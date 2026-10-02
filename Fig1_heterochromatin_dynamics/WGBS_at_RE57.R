## WGBS at RE57 ##
## Script by Nelly Kanata, OWL Schulz ##
## Created on: 23.06.2026
## Modified on: 23.06.2026 ##

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(egg)
library(magrittr)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))



# Read all the .tab files
wdir="./"
setwd(wdir)

# Read WBGS data fro RE57
wgbs <- read.csv("fig_data/WGBS_Meth_RE57_RE58.txt", sep = "\t") %>%
  filter(RE=="RE_57")

wgbs$Day <- wgbs$Time / 24
wgbs$Meth <- wgbs$Meth *100
wgbs$Cells <- factor(wgbs$Cells, levels = c("XO-B7", "XdXIC"))

# read expression data
qPCR <- readxl::read_xlsx("fig_data/qPCR_XO_RG.xlsx")
qPCR %<>%
  pivot_longer(c("R1", "R2"), names_to = "qpcr_rep", values_to = "ct")

qPCR$ct <- if_else(qPCR$ct == "Undetermined", "40", qPCR$ct)

qPCR$ct <- as.numeric(qPCR$ct)

mean.data = qPCR %>% group_by(Sample, gene) %>% summarize(mean.ct = mean(ct), sd.ct = sd(ct))
rel.data <- mean.data %>%
  select(-sd.ct) %>%
  pivot_wider(names_from = gene, values_from = mean.ct) %>% 
  mutate(re_Xist=2^((mean(c(Arpo,Rrm2)))-Xist)) # Rnf12, Zic3 also here, but I dont need them

rel.data <- rel.data %>% select(starts_with("re_")) 


rel.data %<>%
  separate(Sample, c("cell_line", "timepoint", "replicate"), sep = " ") %>%
  mutate(timepoint = as.numeric(str_remove(timepoint, "T")))

rel.data %<>%
  ungroup() %>%
  mutate(percent_max = re_Xist/max(re_Xist)*100) 

rel.data$log2_re <- log2(rel.data$re_Xist)

#scale range
scale_min <- 0
scale_max <- 100

# min-max normalization
rel.data$re_scaled <- scale_min + (rel.data$re_Xist-min(rel.data$re_Xist)) * (scale_max - scale_min)/
  (max(rel.data$re_Xist)-min(rel.data$re_Xist))

rel.data$Cells <- factor(ifelse(rel.data$cell_line =="XO_B07", "XO-B7", "XdXIC"), levels=c("XO-B7", "XdXIC"))


rel.data$Day <- rel.data$timepoint / 24

# plot
plot <- wgbs %>%
  ggplot(aes(Day, Meth)) +
  facet_wrap(~ Cells, axes = "all")+
  
  # WGBS bars
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "#84bcda",
    width = 0.8
  ) +
  
  # Optional: show individual WGBS points
  geom_point(
    aes(Day, Meth),
    color = "black", shape=23,fill="#3b9ad9",
    size = 1,
    position = position_jitter(width = 0.1)
  ) +
  
  scale_y_continuous(
    sec.axis =
      sec_axis(~  min(rel.data$re_Xist) + (. - 0) * (max(rel.data$re_Xist) - min(rel.data$re_Xist)) / (100 - 0),
      name = "Rel. Xist expression" #,
      #breaks = seq(0, 100, by = 50)
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(wgbs$Day))
  ) +
  ylab("Mean DNA methylation (%)") +
  xlab("Differentiation timepoint (days)") +
  
  # qPCR data
  stat_summary(data = rel.data,
               aes(x = Day, y = re_scaled, group = cell_line),
               linetype="dashed",
               fun = mean,
               geom = "line",
               color = "black",
               linewidth = 0.5,
               inherit.aes = FALSE
  )+
  
  geom_point(data = rel.data,
             aes(Day, re_scaled),
             color = "black",shape=21, fill="white", inherit.aes = FALSE)

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)

ggsave("./Figures/RE57_WGBS_vs_Xist.pdf", fix,
       dpi = 300, useDingbats=FALSE)
