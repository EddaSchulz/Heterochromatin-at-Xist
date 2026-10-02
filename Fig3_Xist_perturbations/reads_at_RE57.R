## C&T analysis: reads at REs ##
## Script by Nelly Kanata, OWL Schulz ##
## Created on: 11.08.2023
## Modified on: 03.11.2025 ##

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(egg)

# Setup
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))



# Read all the .tab files
wdir="./"

### XO ####
file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/XO/"), pattern = "\\.tab$", full.names = TRUE)

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
    separate(name, c("Sample", "Dox","dTag", "Antibody","Replicate",  "Sample_Number"), "_", remove = FALSE) 
  

  result_df$KD <- if_else(result_df$dTag == "+dTAG", "ctrl", "KD")
  result_df$KD <- factor(result_df$KD, levels = c("ctrl", "KD"))

  
  # stats
  stat.test <- compare_means(
    mean_value ~KD,
    data = result_df,
    method = "t.test", 
    paired=TRUE  )
  
  
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
  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(3100)
  
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = KD, y = mean_value, ymin=0, group=KD)) +
    geom_point(size=2, aes(color = KD, shape= Replicate)) +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",KD= "#3b9ad9"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "right")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_Xist_KD_XO.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  


  # percent change
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by( KD) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = KD, values_from = mean_across_reps) %>%
    mutate(
      percent_change = ((KD - ctrl) / ctrl) * 100
    )
  summary_df

  ### XX #### 
  
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/XX/"), pattern = "\\.tab$", full.names = TRUE)
  
  # Loop through the list of file names and read each .tab file into a dataframe
  dataframe_list <- list()
  
  for (file_path in file_list) {
    file_name <- tools::file_path_sans_ext(basename(file_path))  # Extract filename without extension
    df <- read.table(file_path, header = TRUE, sep = "\t")  # Read .tab file into a dataframe
    dataframe_list[[file_name]] <- df  # Add dataframe to the list with filename as the name
  }
  
  # Create a table for RE57 with the file names and the means 
  
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
    separate(name, c("Sample","dTag", "Antibody", "Replicate", "Sample_Number"), "_") 
  
  result_df$KD <-if_else(result_df$dTag == "+dTag", "ctrl", "KD")
  result_df$KD <- factor(result_df$KD, levels = c("ctrl", "KD"))
  result_df$Cell_line <-  "Xist KD"

  
  # stats
  stat.test <- compare_means(
    mean_value ~ KD,
    data = result_df,
    method = "t.test", paired=TRUE
  )
  
  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(2000)
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = KD, y = mean_value, ymin=0, group=KD)) +
    geom_point(size=1.5, aes(color = KD, shape= Replicate)) +
    facet_grid(cols=vars(Cell_line), scales = "free_x") +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",KD= "#3b9ad9"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "top")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_Xist_KD_XX.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  
  
  # percent change
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(Cell_line, KD) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = KD, values_from = mean_across_reps) %>%
    mutate(
      percent_change = ((KD - ctrl) / ctrl) * 100
    )
  summary_df
  colnames(summary_df) <- c("Target","ctrl", "KD", "percent_change")
  
  
  