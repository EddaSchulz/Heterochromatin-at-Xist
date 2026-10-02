# Allelic reads
# Nelly Kanata, OWL Schulz
# Created: 04.04.2026

library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(magrittr)
library(ggpubr)


# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

setwd("./")

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

#### Xist ####
# Read data
allelic_data <- readxl::read_xlsx("fig_data/Xist_pyroseq.xlsx", 
  col_names = TRUE)

allelic_data %<>%
  separate(Sample, c("timepoint", "Dox"))

allelic_data$timepoint <- as.factor(gsub("h", "", allelic_data$timepoint))
allelic_data$Day <- if_else(allelic_data$timepoint == "24", "Day 1",
                            if_else(allelic_data$timepoint == "48", "Day 2",
                                    if_else(allelic_data$timepoint == "72", "Day 3", 
                                            if_else(allelic_data$timepoint == "96", "Day 4", NA))))
allelic_data$Dox  <- paste0(allelic_data$Dox, " Dox")

allelic_data$Dox <- if_else(allelic_data$Dox == "24h Dox", "1d Dox", 
                                      if_else(allelic_data$Dox == "48h Dox", "2d Dox", 
                                              allelic_data$Dox))
allelic_data$Dox <- factor(allelic_data$Dox, levels=c("no Dox", "1d Dox", "2d Dox"))

allelic_data %<>%
  pivot_longer(!c(1:3, 9), names_to = "Replicate", values_to = "Cast")

allelic_data$Replicate <- gsub("_Cast", "", allelic_data$Replicate)


# Plot
filtered_allelic_data <- allelic_data 

filtered_allelic_data %<>%
  mutate(B6 = 100 - Cast)


filtered_allelic_data %<>%
  pivot_longer(cols = c(Cast, B6),
               names_to = "allele",
               values_to = "fraction")

filtered_allelic_data$allele <- factor(filtered_allelic_data$allele, levels = c("Cast", "B6"))


stat.test <- compare_means(
  fraction ~ Dox,
  data = filtered_allelic_data %>% filter(allele == "B6"),
  group.by = c("Day"),
  method = "t.test", paired=FALSE
)
stat.test %<>%
  filter(group1 == "no Dox")

stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(100, 115, 100, 115, 100, 115, 100, 115)
stat.test$allele <- "B6"

plot <- filtered_allelic_data %>%
  ggplot(aes(x = Dox, y = fraction, fill=interaction(allele,Dox))) +
  facet_grid(cols = vars(Day)) + 
  stat_summary(geom = "bar", fun.data = mean_se, 
              aes(fill = interaction(allele,Dox), color=interaction(allele,Dox)),
               position = "stack",
               size = 1.1, alpha = 0.8, color = NA, width = 0.8) +
  geom_point(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
             aes(fill = as.factor(allele)), fill="white",shape = 21, color = "grey20",
             alpha = 0.8, size = 1.5, position = position_jitter(width = 0.1, height = 0, seed = 666)) +
  stat_summary(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
               geom = "errorbar", fun.data = mean_se,
               width = 0.05, alpha = 0.8) +
  scale_y_continuous(limits =c(0,120),     breaks = seq(0, 100, 25)) +
  
  scale_fill_manual(values = c(
    "B6.no Dox" = "black",      # main signal
    "Cast.no Dox" = "gray",      # main signal
    
    "B6.1d Dox" = "#e2360e",
    "Cast.1d Dox" = "gray",      # main signal
    
    "B6.2d Dox" = "#F47748",   # background
    "Cast.2d Dox" = "gray"   # background
    
   )) +
  labs(x = "",
       y = expression("Allelic percentage " * bgroup("(",frac(B6~allele, Both~alleles),
                                                ")"))) +
  geom_hline(yintercept = 50, linetype="dashed", color="grey20")+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.border = element_blank()
  )+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2, inherit.aes = FALSE)
  
fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/Pyroseq_Dox_byDay.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



### X-linked genes ####
##### Rnf12 #####
# Read data
allelic_data <- readxl::read_xlsx("fig_data/X-linked-pyroseq.xlsx", 
                                  col_names = TRUE, sheet = "Rnf12")

# replicates
allelic_data$Replicate <- stringr::str_extract(allelic_data$`Sample ID`, "^R\\d+")
allelic_data$Replicate <- if_else(is.na(allelic_data$Replicate), "R4", allelic_data$Replicate)

#fix sample column
allelic_data$`Sample ID` <- gsub("R[0-9]_#[0-9]*_","" , allelic_data$`Sample ID`)
allelic_data$`Sample ID` <- gsub("0h Dox","no Dox" , allelic_data$`Sample ID`)


allelic_data %<>%
  separate(`Sample ID`, c("timepoint", "Dox"), sep = "_")

allelic_data$timepoint <- as.factor(gsub("h", "", allelic_data$timepoint))
allelic_data$Day <- if_else(allelic_data$timepoint == "0", "Day 0",
                            if_else(allelic_data$timepoint == "24", "Day 1",
                            if_else(allelic_data$timepoint == "48", "Day 2",
                                    if_else(allelic_data$timepoint == "72", "Day 3", 
                                            if_else(allelic_data$timepoint == "96", "Day 4", NA)))))

allelic_data$Dox <- if_else(allelic_data$Dox == "24h Dox", "1d Dox", 
                            if_else(allelic_data$Dox == "48h Dox", "2d Dox", 
                                    allelic_data$Dox))
allelic_data$Dox <- factor(allelic_data$Dox, levels=c("no Dox", "1d Dox", "2d Dox"))


allelic_data %<>%
  filter(Quality == "Passed")

colnames(allelic_data)[colnames(allelic_data) == "A  (cast)Frequency (%)"] <- "Cast"

# Plot
filtered_allelic_data <- allelic_data %>%
  filter(! is.na(Dox)) %>%
  filter(Day %in% c("Day 1", "Day 4"))

filtered_allelic_data %<>%
  mutate(B6 = 100 - Cast)


filtered_allelic_data %<>%
  pivot_longer(cols = c(Cast, B6),
               names_to = "allele",
               values_to = "fraction")

filtered_allelic_data$allele <- factor(filtered_allelic_data$allele, levels = c("Cast", "B6"))

# stats
stat.test <- compare_means(
  fraction ~ Dox,
  data = filtered_allelic_data %>% filter(allele == "B6"),
  group.by = c("Day"),
  method = "t.test", paired=FALSE
)
stat.test %<>%
  filter(group1 == "no Dox")

stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(100, 115, 100, 115)
stat.test$allele <- "B6"

plot <- filtered_allelic_data %>%
  ggplot(aes(x = Dox, y = fraction, fill=interaction(allele,Dox))) +
  facet_grid(cols = vars(Day)) + 
  stat_summary(geom = "bar", fun.data = mean_se, 
               aes(fill = interaction(allele,Dox), color=interaction(allele,Dox)),
               position = "stack",
               size = 1.1, alpha = 0.8, color = NA, width = 0.8) +
  geom_point(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
             aes(fill = as.factor(allele)), fill="white",shape = 21, color = "grey20",
             alpha = 0.8, size = 1.5, position = position_jitter(width = 0.1, height = 0, seed = 666)) +
  stat_summary(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
               geom = "errorbar", fun.data = mean_se,
               width = 0.05, alpha = 0.8) +
  scale_y_continuous(limits =c(0,120),     breaks = seq(0, 100, 25)) +
  
  scale_fill_manual(values = c(
    "B6.no Dox" = "black",      # main signal
    "Cast.no Dox" = "gray",      # main signal
    
    "B6.1d Dox" = "#e2360e",
    "Cast.1d Dox" = "gray",      # main signal
    
    "B6.2d Dox" = "#F47748",   # background
    "Cast.2d Dox" = "gray"   # background
    
  )) +
  labs(x = "",
       y = expression("Allelic percentage " * bgroup("(",frac(B6~allele, Both~alleles),
                                                     ")"))) +
  geom_hline(yintercept = 50, linetype="dashed", color="grey20")+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.border = element_blank()
  )+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2, inherit.aes = FALSE)

fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/Pyroseq_Dox_Rnf12_byDay.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


##### Fam122b #####
# Read data
allelic_data <- readxl::read_xlsx("fig_data/X-linked-pyroseq.xlsx", 
                                  col_names = TRUE, sheet = "Fam122b")

# replicates
allelic_data$Replicate <- stringr::str_extract(allelic_data$`Sample ID`, "^R\\d+")
allelic_data$Replicate <- if_else(is.na(allelic_data$Replicate), "R4", allelic_data$Replicate)

#fix sample column
allelic_data$`Sample ID` <- gsub("R[0-9]_#[0-9]*_","" , allelic_data$`Sample ID`)
allelic_data$`Sample ID` <- gsub("0h Dox","no Dox" , allelic_data$`Sample ID`)


allelic_data %<>%
  separate(`Sample ID`, c("timepoint", "Dox"), sep = "_")

allelic_data$timepoint <- as.factor(gsub("h", "", allelic_data$timepoint))
allelic_data$Day <- if_else(allelic_data$timepoint == "0", "Day 0",
                            if_else(allelic_data$timepoint == "24", "Day 1",
                                    if_else(allelic_data$timepoint == "48", "Day 2",
                                            if_else(allelic_data$timepoint == "72", "Day 3", 
                                                    if_else(allelic_data$timepoint == "96", "Day 4", NA)))))

allelic_data$Dox <- if_else(allelic_data$Dox == "24h Dox", "1d Dox", 
                            if_else(allelic_data$Dox == "48h Dox", "2d Dox", 
                                    allelic_data$Dox))
allelic_data$Dox <- factor(allelic_data$Dox, levels=c("no Dox", "1d Dox", "2d Dox"))

allelic_data %<>%
  filter(Quality == "Passed")

colnames(allelic_data)[colnames(allelic_data) == "T  (cast)Frequency (%)"] <- "Cast"

# Plot
filtered_allelic_data <- allelic_data %>%
  filter(!is.na(Dox)) %>%
  filter(Day %in% c("Day 1", "Day 4"))

filtered_allelic_data$Cast <- as.numeric(filtered_allelic_data$Cast)

filtered_allelic_data %<>%
  mutate(B6 = 100 - Cast)


filtered_allelic_data %<>%
  pivot_longer(cols = c(Cast, B6),
               names_to = "allele",
               values_to = "fraction")

filtered_allelic_data$allele <- factor(filtered_allelic_data$allele, levels = c("Cast", "B6"))

# stats
stat.test <- compare_means(
  fraction ~ Dox,
  data = filtered_allelic_data %>% filter(allele == "B6"),
  group.by = c("Day"),
  method = "t.test", paired=FALSE
)
stat.test %<>%
  filter(group1 == "no Dox")

stat.test$pval <- vapply(stat.test$p, format_p, character(1))
stat.test$y.position <- c(100, 115, 100, 115)
stat.test$allele <- "B6"

plot <- filtered_allelic_data %>%
  ggplot(aes(x = Dox, y = fraction, fill=interaction(allele,Dox))) +
  facet_grid(cols = vars(Day)) + 
  stat_summary(geom = "bar", fun.data = mean_se, 
               aes(fill = interaction(allele,Dox), color=interaction(allele,Dox)),
               position = "stack",
               size = 1.1, alpha = 0.8, color = NA, width = 0.8) +
  geom_point(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
             aes(fill = as.factor(allele)), fill="white",shape = 21, color = "grey20",
             alpha = 0.8, size = 1.5, position = position_jitter(width = 0.1, height = 0, seed = 666)) +
  stat_summary(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
               geom = "errorbar", fun.data = mean_se,
               width = 0.05, alpha = 0.8) +
  scale_y_continuous(limits =c(0,120),     breaks = seq(0, 100, 25)) +
  
  scale_fill_manual(values = c(
    "B6.no Dox" = "black",      # main signal
    "Cast.no Dox" = "gray",      # main signal
    
    "B6.1d Dox" = "#e2360e",
    "Cast.1d Dox" = "gray",      # main signal
    
    "B6.2d Dox" = "#F47748",   # background
    "Cast.2d Dox" = "gray"   # background
    
  )) +
  labs(x = "",
       y = expression("Allelic percentage " * bgroup("(",frac(B6~allele, Both~alleles),
                                                     ")"))) +
  geom_hline(yintercept = 50, linetype="dashed", color="grey20")+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.border = element_blank()
  )+
  stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2, inherit.aes = FALSE)

fix <- set_panel_size(plot, height = unit(2.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/Pyroseq_Dox_Fam122b_byDay.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
