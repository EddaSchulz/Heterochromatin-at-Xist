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


### average SNPs #####

# Read data

file_list <- list.files(path = "./fig_data/SNPs/", pattern = "\\.txt$", full.names = TRUE)


ref_file <- read.csv("../Fig2_TX-Maged1-2tag/fig_data/SNPs/RE57_SNPs.bed", sep="\t", header = FALSE)
ref_file$Region <- "RE57"


colnames(ref_file) <- c("Chr", "Start", "End", "Ref/Alt", "Region")
ref_file %<>%
  separate(`Ref/Alt`, sep="/", c("Ref", "Alt")) %>%
  select(-Start, -Chr)


combined_snp_data <- data.frame()

for (file_here in file_list) {
  
  snp_data <- read.csv(file_here, sep="\t", header = TRUE )
  
  snp_data %<>% pivot_longer(c("Nucleotides_R1", "Nucleotides_R2", "Nucleotides_R3"),
                             names_to = "Replicate", values_to = "Nucleotides")
  
  snp_data$Replicate <- gsub("Nucleotides_", "", snp_data$Replicate)
  
  snp_data$Nucleotides <- toupper(gsub("[^ACGTacgt]", "", snp_data$Nucleotides))
  
  snp_data$Sample <- sub(".*/(NK08_[^/]+)_SNPs\\.txt", "\\1", file_here)
  
  
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
  mutate(`B6 (Dox)` = 1 - Cast)


combined_snp_data %<>%
  pivot_longer(cols = c(Cast, `B6 (Dox)`),
               names_to = "allele",
               values_to = "fraction")

combined_snp_data$allele <- factor(combined_snp_data$allele, levels = c("Cast", "B6 (Dox)"))



# Plot

# filter SNPs with few reads
filtered_allelic_data <- combined_snp_data %>%
 filter(Reads > 3) %>%
  group_by(Sample, Region, allele, Replicate) %>%
  summarise(mean_fraction = mean(fraction, na.rm=TRUE),
            sum_Reads=sum(Reads, na.rm=TRUE)) %>%
  ungroup()

summary_df <-filtered_allelic_data %>%
  group_by(Sample, Region, allele) %>%
  summarize(all_mean_fraction = mean(mean_fraction),
            mean_sum_Reads = mean(sum_Reads))

filtered_allelic_data %<>%
  left_join(summary_df, by = c("Sample", "Region", "allele"))

filtered_allelic_data %<>%
  separate(col = Sample, c("Cell_line", "Dox", "dTAG"), sep = "_")

filtered_allelic_data$KD <- if_else(filtered_allelic_data$dTAG == "+dTAG" & filtered_allelic_data$Dox == "-Dox", "ctrl",
                            if_else(filtered_allelic_data$dTAG == "-dTAG" & filtered_allelic_data$Dox == "+Dox", "Xist OE Tsix KD",
                                    if_else(filtered_allelic_data$dTAG == "+dTAG" & filtered_allelic_data$Dox == "+Dox", "Xist OE",
                                            "Tsix KD") ))
filtered_allelic_data$KD <- factor(filtered_allelic_data$KD, 
                           levels=c("ctrl", "Tsix KD", "Xist OE", "Xist OE Tsix KD"))


plot <- filtered_allelic_data %>%
  ggplot(aes(x = KD, y = mean_fraction)) +
  stat_summary(geom = "bar", fun.data = mean_se, 
               aes(fill = allele, alpha = sum_Reads),
               position = "stack",
               size = 1.1,  color = NA, width = 0.7) +
  geom_point(data=filtered_allelic_data[filtered_allelic_data$allele == "B6 (Dox)",],
             aes(fill = "grey20"), shape = 21, color = "grey20",
             alpha = 0.8, size = 1, position = position_jitter(width = 0.1, height = 0, seed = 666)) +
  stat_summary(data = filtered_allelic_data[filtered_allelic_data$allele == "B6 (Dox)",],
               geom = "errorbar", fun.data = mean_se,
               width = 0.05, alpha = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  geom_hline(yintercept = 0.5, linetype="dashed", color="grey20")+
  scale_fill_manual(values =c("B6 (Dox)"="#F47748","Cast"= "gray80")) +
  labs(x = "",
       y="Allelic ratio")+
  
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.4, 'cm')
  )

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/allelic_reads_XOE_averageSNPs.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
