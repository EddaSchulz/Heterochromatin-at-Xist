## C&T analysis: reads at REs ##
## Script by Nelly Kanata, OWL Schulz ##
## Created on: 11.08.2023
## Modified on: 03.11.2025 ##

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(egg)
library(magrittr)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  #axis.text.x = element_text(angle = 45, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))



# Read all the .tab files
wdir="./"
setwd(wdir)

file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/Gjaltema_etal/"), pattern = "\\.tab$", full.names = TRUE)

# Loop through the list of file names and read each .tab file into a dataframe
dataframe_list <- list()

for (file_path in file_list) {
  file_name <- tools::file_path_sans_ext(basename(file_path))  # Extract filename without extension
  df <- read.table(file_path, header = TRUE, sep = "\t")  # Read .tab file into a dataframe
  dataframe_list[[file_name]] <- df  # Add dataframe to the list with filename as the name
}

# Create a table for RE57 with the file names and the means 

element<-"CRE_57"
  
  result_list <- list() # Create an empty list to store the results
  
  # Loop through the list of dataframes
  for (df_name in names(dataframe_list)) {
    df <- dataframe_list[[df_name]]
    
    # Extract the row where "name" is "RE57"
    element_row <- df[df$name == element, ]
    
    # Extract the mean value from the row (assuming the "mean" column is named "mean")
    element_mean <- element_row$mean
    
    # Add the result to the result_list
    result_list[[df_name]] <- data.frame(name = df_name, mean_value = element_mean)
  }
  
  # Combine the individual dataframes in the result_list into a single dataframe
  result_df <- do.call(rbind, result_list)
  
  
  # Split the file names by Sample, Day, Replicate, Antibody, Sample Number
  result_df <- result_df %>%
    separate(name, c("Sample", "Antibody", "timepoint", "Replicate"), "_", remove = FALSE)
  
 result_df$timepoint <- as.numeric(gsub("t", "", result_df$timepoint))
  
 result_df$Day <- result_df$timepoint / 24
  

 
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
   mutate(re_Xist=2^((mean(c(Arpo,Rrm2)))-Xist)) 
 
 rel.data <- rel.data %>% select(starts_with("re_")) 
 
 
rel.data %<>%
   separate(Sample, c("cell_line", "timepoint", "replicate"), sep = " ") %>%
   mutate(timepoint = as.numeric(str_remove(timepoint, "T")))

rel.data$log2_re <- log2(rel.data$re_Xist)

#scale range
scale_min <- 0
scale_max <- max(result_df$mean_value)

# min-max normalization
rel.data$re_scaled <- scale_min + (rel.data$re_Xist-min(rel.data$re_Xist)) * (scale_max - scale_min)/
  (max(rel.data$re_Xist)-min(rel.data$re_Xist))

rel.data$Sample <- ifelse(rel.data$cell_line =="XO_B07", "XO", "XX")

rel.data$Day <- rel.data$timepoint / 24

 # plot
 plot <- result_df %>%
   ggplot( aes(Day, mean_value, group = Sample)) +
   facet_wrap(~ Sample, axes = "all")+
   # Mean line
   stat_summary(
     fun = mean,
     geom = "line",
     linewidth = 0.5,
     color = "#3b9ad9"
   ) +
   geom_point(
     aes(Day, mean_value),shape=21,
     fill = "#3b9ad9", size=1.5, color="black")+

   
   scale_y_continuous(
     sec.axis =
       sec_axis(~  min(rel.data$re_Xist) + (. - 0) * (max(rel.data$re_Xist) - min(rel.data$re_Xist)) / (scale_max - 0),
                name = "Rel. Xist expression" 
       ))+
   scale_x_continuous(
     breaks = sort(unique(result_df$Day))
   )+
   ylab("Norm. counts over RE57") +
   xlab("Differentiation timepoint (days)")+
   stat_summary(data = rel.data,
     aes(x = Day, y = re_scaled, group = Sample),
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
 
 ggsave("./Figures/RE57_H3K9me3_vs_Xist.pdf", fix,
        dpi = 300, useDingbats=FALSE)
 
  