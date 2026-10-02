## qPCR Analysis
## Nelly Kanata ~ OWL Schulz
## Script adjusted from Alexandra Martitz
## Created: 06.04.2023
## Modified: 07.11.2025

library(tidyverse)
library(magrittr)
library(ggplot2)
library(readxl) 


wrk_dir <- "./" # change this to your working directory
setwd(wrk_dir)

# Path to your qpcr results.
excel_1 <- "2026-04-29 104240_20260429 150040.xlsx" # example file


# Load data
raw_data <- read_xlsx(excel_1, sheet = "Results", skip = 24) # the first 24 lines are metadata


# Read sample index
sample_index <- read_xlsx("exp1-2_qPCR_template.xlsx", sheet = "Sheet2",col_names = FALSE) # example file with sample metadata
colnames(sample_index) <- c("Sample_Nr", "Replicate", "timepoint", "Dox")

# pre-process metadata
sample_index$timepoint <- as.numeric(gsub("h","", sample_index$timepoint ))
sample_index$Sample_Nr <- paste0("Sample ", sample_index$Sample_Nr)

# Add metadata info to raw_data table
raw_data <- merge(raw_data, sample_index, by.x = "Sample", by.y = "Sample_Nr", all.x = TRUE)

raw_data %<>%
  filter(!is.na(Sample))

# Keep only the relevant columns
data = raw_data %>%
  select('Sample',"Replicate","timepoint", "Dox", 'Target','Cq', 'Task') %>%
  group_by(Sample, Target) %>%
  mutate(rep_qPCR = row_number()) %>% # technical replicates
  dplyr::rename(Sample = Sample, gene = Target, ct = Cq , task = Task)

# Make sure that you do not get a non-specific signal in the no-template control
data %>% filter(task=='NTC')


# convert ct values to the right class.
sel.data = data %>% filter(task!="NTC") %>% select(-task)
sel.data$ct <- if_else(sel.data$ct == "Undetermined", "40", sel.data$ct) # detection threshold
sel.data$ct <- as.numeric(sel.data$ct)


# Visualize how comparable the technical replicates are. They should ideally vary by <0.5 CT values

sel.data %>%
  ggplot(aes(x=Sample,y=ct)) + 
  geom_point() + facet_wrap(vars(gene))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("ct_values_technical_replicates.pdf",
       dpi = 300,
       width = 10,
       height = 6)


# Calculate mean and standard deviation across technical replicates.
mean.data = sel.data %>% group_by(Sample, Replicate, timepoint, Dox, gene) %>% summarize(mean.ct = mean(ct), sd.ct = sd(ct))

# Plot the standard deviation.
mean.data %>% ggplot(aes(x=Sample,y=sd.ct)) + 
  geom_point() +
  facet_wrap(vars(gene)) + 
  ylab('CT Std. Dev.' )+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("sd_technical_replicates.pdf",
       dpi = 300,
       width = 10,
       height = 6)




# Calculate gene expression relative to the house keeping genes Arpo and Rrm2.
rel.data <- mean.data %>%
  select(-sd.ct) %>%
  pivot_wider(names_from = gene, values_from = mean.ct) %>% 
  mutate(re_Xist=2^((mean(c(Arpo,Rrm2)))-Xist)) # adjust for targeted genes

rel.data <- rel.data %>% select(starts_with("re_")) 


rel.data_long <- rel.data %>%
  pivot_longer(
    cols = matches("^re_"),   # only columns that start with "re_"
    names_to = "gene",
    values_to = "re"
  )

rel.data_long$gene <- gsub("re_", "", rel.data_long$gene)



##### relative expression in log2 ##### 

rel.data_long$log2_re <- log2(rel.data_long$re)



### Save table ####

write.csv(rel.data_long, "qPCR_data.csv",
          row.names = FALSE)


# Session Info
sink("SessionInfo.txt", append = FALSE)
cat(date())
cat("\n")
cat("\n")
cat(getwd())
cat("\n")
cat("\n")
print(sessionInfo())
sink()
