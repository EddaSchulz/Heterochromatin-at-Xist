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


wdir="./"

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


### Xist OE ####

file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/Xist_OE/"), pattern = "\\.tab$", full.names = TRUE)

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
    separate(name, c("Sample", "Dox", "Antibody","Replicate",  "Sample_Number"), "_", remove = FALSE) 
  

  result_df$OE <- if_else(result_df$Dox == "2i+Dox" | result_df$Dox == "diff+Dox", "OE", "ctrl")
  result_df$OE <- factor(result_df$OE, levels = c("ctrl", "OE"))
  
  result_df$timepoint <- if_else(result_df$Dox == "diff+Dox" | result_df$Dox == "diff-Dox", "Day 2", "ESCs")
  result_df$timepoint <- factor(result_df$timepoint, levels = c("ESCs", "Day 2"))
  
  
  result_df %<>%
    #filter(Sample == "dXIC") %>%
    arrange(timepoint, OE, Replicate)
  
  

  # stats
  stat.test <- compare_means(
    mean_value ~ OE,
    data = result_df[result_df$timepoint == "Day 2",],
    method = "t.test", 
    paired=TRUE  )
  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(3100)
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    filter(timepoint == "Day 2") %>%
    ggplot(aes(x = OE, y = mean_value, ymin=0, group=OE)) +
    geom_point(size=1.5, aes(color = OE, shape= Replicate)) +
    facet_grid(cols=vars(Sample), scales = "free_x") +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",OE= "#f37748"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "right")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_XistOE_diff.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  # stats
  stat.test <- compare_means(
    mean_value ~ OE,
    data = result_df[result_df$timepoint == "ESCs",],
    method = "t.test", 
    paired=TRUE  )
  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(3100)
  
  plot <- result_df %>%
    filter(timepoint == "ESCs") %>%
    ggplot(aes(x = OE, y = mean_value, ymin=0, group=OE)) +
    geom_point(size=1.5, aes(color = OE, shape= Replicate)) +
    facet_grid(cols=vars(Sample), scales = "free_x") +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",OE= "#f37748"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "right")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_XistOE_ESCs.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  

  
  # percent change
  # (mean of OE - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(OE, timepoint) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = OE, values_from = mean_across_reps) %>%
    mutate(
      percent_change = ((OE - ctrl) / ctrl) * 100
    )
  summary_df

  
  ### Tsix KDs #####
  
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/Tsix_KD/"), pattern = "\\.tab$", full.names = TRUE)
  
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
  result_df$Cell_line <- "Tsix KD"

  
# stats
  stat.test <- compare_means(
    mean_value ~ KD,
    data = result_df_tsix,
    method = "t.test",paired=TRUE
  )

  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(2600)

  
  
  
  # Create the plot
  
  plot <- result_df %>%
    filter(Cell_line == "Tsix KD") %>%
    ggplot(aes(x = KD, y = mean_value, ymin=0, group=KD)) +
    geom_point(size=1.5, aes(color = KD, shape= Replicate)) +
    facet_grid(cols=vars(Cell_line), scales = "free_x") +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",KD= "#3b9ad9"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "top")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_Tsix_KD.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  

 
 ### Xist OE Tsix dTag #####

  
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/XistOE_TsixKD/"), pattern = "\\.tab$", full.names = TRUE)
  
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
    separate(name, c("Sample","Dox", "dTAG", "Antibody", "Replicate", "Sample_Number"), "_") 
  
  
  result_df$KD <- if_else(result_df$Dox == "-Dox" & result_df$dTAG == "+dTAG", "ctrl",
                    if_else( result_df$Dox == "-Dox" & result_df$dTAG == "-dTAG", "Tsix KD",
                       if_else(result_df$Dox == "+Dox" & result_df$dTAG == "+dTAG", "Xist OE",
                          if_else(result_df$Dox == "+Dox" & result_df$dTAG == "-dTAG", "Xist OE Tsix KD",
                                  ""
                                        ))
                                  ))
  result_df$KD <- factor(result_df$KD, levels = c("ctrl", "Tsix KD", "Xist OE", "Xist OE Tsix KD"))
  
  #  result_df$Sample_Number <- gsub(".noAmpR", "", result_df$Sample_Number)
  
  result_df$Cell_line <- "Xist OE Tsix KD"
  
  result_df %<>%
    arrange(KD, Replicate)
  
  # stats
  stat.test <- compare_means(
    mean_value ~ KD,
    data = result_df,
    method = "t.test", paired=TRUE  )
  stat.test %<>%
    filter(group1 == "ctrl")
  
  stat.test$pval <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(800,1150 ,1400 )
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = KD, y = mean_value, ymin=0, group=KD)) +
    geom_point(size=1.5, aes(color = KD, alpha=dTAG, shape= Replicate)) +
    facet_grid(cols=vars(Cell_line), scales = "free_x") +
    stat_pvalue_manual(stat.test, label = "pval", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean signal over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b","Tsix KD" = "#3f414b", 
                                  "Xist OE"= "#f37748", "Xist OE Tsix KD"= "#f37748" ))+
    scale_alpha_manual(values=c(0.6, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/H3K9me3_RE57_XistOE_TsixKD.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  

  # percent change
  
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(KD, Cell_line) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    ungroup() %>%
    mutate(ctrl_value = mean_across_reps[ KD == "ctrl" ]) %>%
    # Calculate percent change relative to ctrl
    mutate(percent_change_vs_ctrl = ((mean_across_reps - ctrl_value) / ctrl_value) * 100)
  
  summary_df <- summary_df %>%
    filter(KD !="ctrl") %>%
    select(-Cell_line) %>%
    select(1,3,2,4)
  
  summary_df

  