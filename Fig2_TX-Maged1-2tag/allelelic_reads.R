# Allelic reads
# Nelly Kanata, OWL Schulz
# Created: 20.03.2026

library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(magrittr)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

setwd("./")


# Read data

unsorted <- read.csv("fig_data/SNPs/unsorted_SNPs.txt", sep="\t", header = TRUE)
P1 <- read.csv("fig_data/SNPs/P1_SNPs.txt", sep="\t", header = TRUE)
P3 <- read.csv("fig_data/SNPs/P3_SNPs.txt", sep="\t", header = TRUE)
P4 <- read.csv("fig_data/SNPs/P4_SNPs.txt", sep="\t", header = TRUE)

ref_file <- read.csv("fig_data/SNPs/RE57_SNPs.bed", sep="\t", header = FALSE)
ref_file$Region <- "RE57"

colnames(ref_file) <- c("Chr", "Start", "End", "Ref/Alt", "Region")
ref_file %<>%
  separate(`Ref/Alt`, sep="/", c("Ref", "Alt")) %>%
  select(-Start, -Chr)


combined_snp_data <- data.frame()

for (gate in c("P1", "P3", "P4", "unsorted")) {
  
  snp_data <- get(gate)
  
  snp_data %<>% pivot_longer(c("Nucleotides_R1", "Nucleotides_R2", "Nucleotides_R3"),
                             names_to = "Replicate", values_to = "Nucleotides")
  
  snp_data$Replicate <- gsub("Nucleotides_", "", snp_data$Replicate)
  
  snp_data$Nucleotides <- toupper(gsub("[^ACGTacgt]", "", snp_data$Nucleotides))
  
  snp_data$gate <- gate
  
  
  combined_snp_data <- rbind(combined_snp_data, snp_data)
  
}

for (b in c("A","C","G","T")) {
  combined_snp_data[[paste0(b, "_count")]] <- nchar(gsub(paste0("[^", b, "]"), "", combined_snp_data$Nucleotides))
}



combined_snp_data <- merge(combined_snp_data, ref_file, by.x = "Start", by.y = "End")


combined_snp_data$Reads <- nchar(combined_snp_data$Nucleotides)


combined_snp_data <- combined_snp_data %>%
  rowwise() %>%
  mutate(Cast = get(paste0(Alt, "_count")) / Reads) %>%
  ungroup()


combined_snp_data %<>%
  mutate(B6 = 1 - Cast)


combined_snp_data %<>%
  pivot_longer(cols = c(Cast, B6),
               names_to = "allele",
               values_to = "fraction")

combined_snp_data$allele <- factor(combined_snp_data$allele, levels = c("Cast", "B6"))
combined_snp_data$fraction <- combined_snp_data$fraction *100 # make it percentage


# filter SNPs with few reads
filtered_allelic_data <- combined_snp_data %>%
 filter(Reads > 3) %>%
  group_by(gate, Region, allele, Replicate) %>%
  summarise(mean_fraction = mean(fraction, na.rm=TRUE),
            sum_Reads=sum(Reads, na.rm=TRUE)) %>%
  ungroup()

summary_df <-filtered_allelic_data %>%
  group_by(gate, Region, allele) %>%
  summarize(all_mean_fraction = mean(mean_fraction),
            mean_sum_Reads = mean(sum_Reads))

filtered_allelic_data %<>%
  left_join(summary_df, by = c("gate", "Region", "allele"))



filtered_allelic_data %<>%
  filter(Region == "RE57") 

filtered_allelic_data$gate <- factor(filtered_allelic_data$gate, levels = c("unsorted", "P1", "P3", "P4"))

filtered_allelic_data$Gate <- if_else(filtered_allelic_data$gate == "P1", "Xist Off",
                          if_else(filtered_allelic_data$gate == "P3", "Xist-Cast", 
                                  if_else(filtered_allelic_data$gate == "P4", "Xist-B6", "unsorted (mix)")))
filtered_allelic_data$Gate <- factor(filtered_allelic_data$Gate, levels = c("unsorted (mix)", "Xist Off", "Xist-Cast","Xist-B6"))


plot <- filtered_allelic_data %>%
  ggplot(aes(x = Gate, y = mean_fraction)) +
  facet_grid(cols=vars(Region))+
  stat_summary(geom = "bar", fun.data = mean_se, 
               aes(fill = allele, alpha = sum_Reads),
               position = "stack",
               size = 1.1,  color = NA, width = 0.8) +
  geom_point(data=filtered_allelic_data[filtered_allelic_data$allele == "B6",],
             aes(fill = "black"), color = "black",
             alpha = 0.7, size = 1, position = position_jitter(width = 0.1, height = 0, seed = 666)) +
  stat_summary(data = filtered_allelic_data[filtered_allelic_data$allele == "B6",],
               geom = "errorbar", fun.data = mean_se,
               width = 0.05, alpha = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 50, linetype="dashed", color="grey20")+
  scale_alpha(range = c(0.3, 1))+
  scale_fill_manual(values =c(
    
  "Cast" = "#cfdaf0",
  "B6"  = "#3b9ad9"
))+
  labs(x = "",
       y="Allelic reads (%)")+
  
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
  )

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/allelic_reads_sorted_cells_averageSNPs.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
