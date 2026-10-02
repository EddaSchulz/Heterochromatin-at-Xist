## C&T analysis: reads at REs ##
## Script by Nelly Kanata, OWL Schulz ##
## Created on: 11.08.2023
## Modified on: 24.06.2026 ##

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(egg)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 45, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))



# Read all the .tab files
wdir="./"
setwd(wdir)

file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/sorting/"), pattern = "\\.tab$", full.names = TRUE)

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
    separate(name, c("Sample", "Replicate", "Antibody", "Sample_Number"), "_", remove = FALSE)
  
  result_df$Gate <- if_else(result_df$Sample == "P1", "Xist Off",
                            if_else(result_df$Sample == "P3", "Xist-Cast", 
                            if_else(result_df$Sample == "P4", "Xist-B6", "unsorted (mix)")))
 result_df$Gate <- factor(result_df$Gate, levels = c("unsorted (mix)", "Xist Off", "Xist-Cast","Xist-B6"))

 
 
 # check significance
 # 4 samples -> ANOVA
 res_aov <- aov(mean_value ~ Sample,
                data = result_df
 )
 
 summary(res_aov)

 
 
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = Gate, y = mean_value, group=Gate)) +
    geom_point(size=1, aes(color = Sample, shape= Replicate)) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.1))) +
    labs(x = "", y = "Mean norm. counts over RE57") +
    stat_summary(geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.5),
                 color = c("grey20", "#f37748", "#e2360e", "#ebc31d")) +
    scale_color_manual(values = c("P1" ="#3f414b", "P3"= "#3f414b", "P4" = "#3f414b", "unsorted" = "#3f414b"))+
    
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/sorted_H3K9me3_RE57.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  

  