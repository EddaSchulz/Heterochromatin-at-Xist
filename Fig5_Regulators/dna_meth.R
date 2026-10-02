## DNA methylation Plots
## Nelly Kanata ~ OWL Schulz
## Created: 26.10.2025
## Modified: 04.02.2026

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

#### Plot rnf12_zic3 ####
# load data
rnf12_zic3 <- read.csv("./fig_data/Rnf12_Zic3_all_methylation_table.csv")

RE57_amplicons <- c("RE57_YY1", "RE57_CGI")

filtered_df <- rnf12_zic3 %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") %>%
  filter(KD %in% c("Rnf12", "Zic3")) 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$condition <- factor(if_else(filtered_df$condition=="WT", "ctrl", "KD"), levels = c("ctrl", "KD"))


### mean methylation plot
methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, condition, KD,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(KD)
methC <- as.data.frame(methC)
methC$timepoint <- "D3"

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))
methC$mean_meth <- methC$mean_meth * 100

methC$KD <- paste0(methC$KD, " KD")

# calculate significance
test_df <- expand.grid(KD = unique(methC$KD),
                       timepoint = unique(methC$timepoint))

test_df$pval <- mapply(function(a,b)
{t.test(methC[methC$KD==a &
                methC$timepoint==b &
                methC$condition=="ctrl",]$mean_meth,
        methC[methC$KD==a &
                methC$timepoint==b &
                methC$condition=="KD",]$mean_meth,
        paired=TRUE
)$p.value},
test_df$KD, test_df$timepoint)

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

test_df$condition <- "KD"


plot <- ggplot(methC, aes( condition, mean_meth)) +
  facet_grid(cols=vars(timepoint), rows=vars(KD)) +
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=condition),
    alpha=0.6,
    #fill = c("#a7a9ac", "#84bcda"),
    width = 0.7
  ) +
  
  geom_point(aes(fill = condition), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values=c(ctrl="#3f414b", KD="#F47748")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 25)) +  # Extends to 1.25 but labels stop at 1
  ylab("Mean DNA methylation (%)")+
  xlab("")+
  ggtitle(paste0("")) +
  geom_text(data = test_df,
            aes(label = pval, y = 100), size=6/2.8,
            nudge_x = -0.5, color="black")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)


ggsave("./Figures/mean_meth_Rnf12_Zic3.pdf", fix,
       dpi = 300, useDingbats=FALSE)


##### Coordinate plot ####
plot <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>%
  filter(KD=="Rnf12") %>%
  ggplot(aes(x = start, y = ratio*100, color = condition)) +
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.4, aes(fill = condition), color = NA) +
  facet_grid( cols=vars(amplicon), scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("WT" = "gray","KD" = "#F47748"))+
  scale_fill_manual(values = c("WT" = "gray","KD" = "#F47748"))+
  ylim(c(0, 100)) +
  ylab("% methylated)") + xlab("amplicon coordinates (arbitrary)")+
  theme(panel.grid = element_blank()) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(shape = timepoint),  
    size = 1.5,
    fill = "white",          # ensures open symbols look nice
    position = position_dodge(width = 0.5)
  ) 


ggsave(paste0("./Figures/DNAme_dTag_Rnf12.pdf"),  plot,
       dpi = 300, useDingbats=FALSE, height = unit(1.5, "cm"), width = unit(3.5, "cm"))



#### Plot zfp280c_zfp36l1 #####
# load data

zfp280c_zfp36l1 <- read.csv("./fig_data/Zfp280c_Zfp36l1_all_methylation_table.csv")
rex1_zfp281 <- read.csv("./fig_data/Zfp281_Myc_Rex1_all_methylation_table.csv")


filtered_df <- zfp280c_zfp36l1 %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") %>%
  filter(KD %in% c("Zfp280c", "Zfp36l1")) 

filtered_df$KD <- paste0(filtered_df$KD, " KD")

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$timepoint <- factor(filtered_df$timepoint, levels = c("D2", "D4"))
filtered_df$dTag <- factor(filtered_df$dTag, levels = c("+", "-"))
filtered_df$dTag <- factor(if_else(filtered_df$dTag=="+", "ctrl", "KD"), levels = c("ctrl", "KD"))

### mean methylation plot
RE57_amplicons <- c("RE57_YY1", "RE57_CGI")

methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, dTag, KD,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(KD)
methC <- as.data.frame(methC)
methC$timepoint <- if_else(methC$timepoint == "D2", "Day 2", "Day 4")
methC$timepoint <- factor(methC$timepoint, levels = c("Day 2", "Day 4"))

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))

methC$mean_meth <- methC$mean_meth * 100

# calculate significance
test_df <- expand.grid(KD = unique(methC$KD),
                       timepoint = unique(methC$timepoint))

test_df$pval <- mapply(function(a,b)
{t.test(methC[methC$KD==a &
                methC$timepoint==b &
                methC$dTag=="ctrl",]$mean_meth,
        methC[methC$KD==a &
                methC$timepoint==b &
                methC$dTag=="KD",]$mean_meth,
        paired=TRUE
)$p.value},
test_df$KD, test_df$timepoint)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$dTag <- "KD"


plot <- ggplot(methC, aes( dTag, mean_meth)) +
  facet_grid(cols=vars(timepoint), rows=vars(KD)) +
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=dTag),
    alpha=0.6,
    width = 0.7
  ) +
  
  geom_point(aes(fill = dTag), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values=c(ctrl="#3f414b", KD="#3b9ad9")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 25)) +  # Extends to 1.25 but labels stop at 1
  ylab("Mean DNA methylation (%)")+
  xlab("")+
  ggtitle(paste0("")) +
  geom_text(data = test_df,
            aes(label = pval, y = 105), size=6/2.8,
            nudge_x = -0.5, color="black")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave("./Figures/mean_meth_Zfp280c_Zfp36l1.pdf", fix,
       dpi = 300, useDingbats=FALSE)


#### Plot rex1 ####

filtered_df <- rex1_zfp281 %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") %>%
  filter(KD %in% c("Rex1", "NTC")) 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$timepoint <- factor(filtered_df$timepoint, levels = c("D2", "D4"))


methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, KD,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(KD)
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(6:ncol(methC)), na.rm = TRUE))

methC$mean_meth <- methC$mean_meth * 100


# Rex1 table (excluding data for individual CpGs)
methC$condition <- if_else(methC$KD == "NTC", "ctrl", "KD")
methC$KD <- "Rex1"
rex1 <- methC



#### Plot zfp281 #####

filtered_df <- rex1_zfp281 %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") %>%
  filter(KD %in% c("Zfp281", "NTC")) 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))

# mean meth plot

methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, KD,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(KD)
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(6:ncol(methC)), na.rm = TRUE))


methC$KD <- factor(methC$KD,levels = c("NTC", "Rex1", "Zfp281"))

methC$mean_meth <- methC$mean_meth * 100
methC$condition <- if_else(methC$KD == "NTC", "ctrl", "KD")
methC$KD <- "Zfp281"
methC <- rbind(methC, rex1)
methC$KD <- paste0(methC$KD, " KD")

methC$timepoint <- if_else(methC$timepoint == "D2", "Day 2", "Day 4")
methC$timepoint <- factor(methC$timepoint, levels = c("Day 2", "Day 4"))


# calculate significance
test_df <- expand.grid(KD = unique(methC$KD),
                       timepoint = unique(methC$timepoint))

test_df$pval <- mapply(function(a, b)
{t.test(methC[methC$KD==a &
                methC$timepoint==b &
                methC$condition=="ctrl",]$mean_meth,
        methC[methC$KD==a &
                methC$timepoint==b &
                methC$condition=="KD",]$mean_meth,
        paired=TRUE
)$p.value},
test_df$KD, test_df$timepoint)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$condition <- "KD"



plot <- ggplot(methC, aes( condition, mean_meth)) +
  facet_grid(cols=vars(timepoint), rows = vars(KD)) +
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=condition),
    alpha=0.6,
    width = 0.7
  ) +
  
  geom_point(aes(fill = condition), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values=c(ctrl="#3f414b", KD="#3b9ad9")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 25)) +  
  ylab("Mean DNA methylation (%)")+
  xlab("")+
  geom_text(data = test_df,
            aes(label = pval, y = 105), size=6/2.8,
            nudge_x = -0.5, color="black")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave("./Figures/mean_meth_Rex1_Zfp281.pdf", fix,
       dpi = 300, useDingbats=FALSE)


#### Plot RNF12-KO/ DKO ######
# load data

dko_data <- read.csv("./fig_data/DKO_all_methylation_table.csv")
dko_data$cell_line <- if_else(dko_data$cell_line == "RNF12", "Rnf12 KO", dko_data$cell_line)


filtered_df <- dko_data %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57")

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))


methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff ## THIS IS NORMALLY 10 BUT THIS LIBRARY IS UNDERSEQUENCED
  group_by(sample_number) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, cell_line, condition,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) 
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(6:ncol(methC)), na.rm = TRUE))
methC$cell_line <- factor(methC$cell_line, levels = c( "WT", "Rnf12 KO", "DKO"))

methC$mean_meth <- methC$mean_meth *100

my_comparisons <- list( c("WT", "Rnf12 KO"), c("WT", "DKO"), c("Rnf12 KO", "DKO") )

# stats
stat.test <- compare_means(
  mean_meth ~ cell_line,
  data = methC,
  method = "t.test"
)

stat.test$pval <- vapply(stat.test$p, format_p, character(1))

stat.test$y.position <- c(120, 140, 105)

plot <- ggplot(methC, aes( cell_line, mean_meth)) +
  # bars
  stat_summary(
    fun = mean,
    geom = "col",
    aes(fill=cell_line),
    alpha=0.6,
    width = 0.7
  ) +
  
  geom_point(aes(fill = cell_line), shape = 23, color="black", size = 1, alpha = 1)+
  scale_fill_manual(values = c("WT" ="#3f414b",`Rnf12 KO`= "#f37748", DKO="#751f58"))+
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 100, by = 20)) +  # Extends to 150 but labels stop at 1
  ylab("Mean DNA methylation (%)")+
  xlab("")+
  ggtitle(paste0("")) +
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)

fix <- set_panel_size(plot, height = unit(2.4, "cm"), width = unit(1.8, "cm"))
grid.arrange(fix)

ggsave("./Figures/mean_meth_dko_data.pdf", fix,
       dpi = 300, useDingbats=FALSE)

