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


setwd("./")

# Read all the .tab files
wdir="./"

  ### repressors dTag #####
  
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/repressors_dTag/"), pattern = "\\.tab$", full.names = TRUE)
  
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
    separate(name, c("Sample", "Replicate", "Antibody", "Sample_Number"), "_") %>%
    
    mutate(dTAG = substr(Sample, nchar(Sample), nchar(Sample)),
           KD_guides = substr(Sample, 1, nchar(Sample) - 1))
  
  result_df$KD <- if_else(result_df$KD_guides == "TS35", "Zfp280c KD", 
                          if_else(result_df$KD_guides == "TS37", "Setdb1 KD",
                                  if_else(result_df$KD_guides == "TS38", "Zfp36l1 KD","")))
  result_df$Sample_Number <- gsub(".noAmpR", "", result_df$Sample_Number)
  
  result_df$Cell_line <- if_else(result_df$dTAG == "+", "ctrl", "KD")
  result_df$Cell_line <- factor(result_df$Cell_line, levels = c("ctrl", "KD"))
  
  
  # calculate significance
  test_df <- expand.grid(KD = unique(result_df$KD))
  
  test_df$pval <- mapply(function(a,b)
  {t.test(result_df[result_df$KD==a &
                      result_df$Cell_line=="ctrl",]$mean_value,
          result_df[result_df$KD==a &
                      result_df$Cell_line=="KD",]$mean_value,
          paired=TRUE
  )$p.value},
  test_df$KD)
  
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
  
  test_df$Cell_line <- "KD"
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = Cell_line, y = mean_value, ymin=0, group=Cell_line)) +
    geom_point(size=1.5, aes(color = Cell_line, shape= Replicate)) +
    facet_grid(rows=vars(KD), scales = "free") +
    geom_text(data = test_df,
              aes(label = pval, y = 5350), size=6/2.8,
              nudge_x = -0.5, color="black")+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",KD= "#3b9ad9"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "top")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/repressors_dTag_H3K9me3_RE57.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  

  # percent change
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(KD, dTAG) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = dTAG, values_from = mean_across_reps) %>%
    mutate(
      percent_change = ((`-` - `+`) / `+`) * 100
    )
  summary_df
  
  
 ### repressors ABA #####
  
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/repressors_ABA/"), pattern = "\\.tab$", full.names = TRUE)
  
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
    separate(name, c("KD_guides", "Replicate", "Antibody", "Sample_Number"), "_") 
  
  
  result_df$KD <- if_else(result_df$KD_guides == "LR15", "NTC", 
                          if_else(result_df$KD_guides == "NK06", "Rex1 KD",
                                  if_else(result_df$KD_guides == "NK07", "Myc KD",
                                          if_else(result_df$KD_guides == "NK11", "Zfp281 KD",""))))
  result_df$KD <- factor(result_df$KD, levels = c("NTC", "Myc KD", "Rex1 KD", "Zfp281 KD"))
  

  result_df$Cell_line <- if_else(result_df$KD == "NTC", "NTC", "KD")
  
  # rearrange so that i plot rex1 and Zfp281 separately
  temp_df <- result_df[result_df$Cell_line == "NTC",]
  temp_df$KD <-"Rex1 KD"
  
  result_df$KD[result_df$Cell_line == "NTC"]<-"Zfp281 KD"
  result_df <- rbind(result_df, temp_df)
  

    # calculate significance
  test_df <- expand.grid(KD = unique(result_df$KD))
  
  test_df$pval <- mapply(function(a,b)
  {t.test(result_df[result_df$KD==a &
                      result_df$Cell_line=="NTC",]$mean_value,
          result_df[result_df$KD==a &
                      result_df$Cell_line=="KD",]$mean_value,
          paired=TRUE
  )$p.value},
  test_df$KD)
  
  test_df$pval <- vapply(test_df$pval, format_p, character(1))
  
  test_df$Cell_line <- "KD"
  
  result_df$Cell_line <- if_else(result_df$Cell_line == "NTC", "ctrl", result_df$Cell_line)
  result_df$Cell_line <- factor(result_df$Cell_line, levels=c("ctrl", "KD"))
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = Cell_line, y = mean_value, ymin=0, group=KD)) +
    geom_point(size=1.5, aes(color = Cell_line, shape= Replicate)) +
    facet_grid(rows=vars(KD), scales = "free") +
    geom_text(data = test_df,
              aes(label = pval, y = 5350), size=6/2.8,
              nudge_x = -0.5, color="black")+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b","KD"= "#3b9ad9"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/repressors_ABA_H3K9me3_RE57.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  

  # percent change
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(KD, Cell_line) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    mutate(ntc_value = mean_across_reps[ Cell_line == "ctrl" ]) %>%
    ungroup() %>%
    # Calculate percent change relative to NTC
    mutate(percent_change_vs_NTC = ((mean_across_reps - ntc_value) / ntc_value) * 100)
  
  summary_df 
  
  # activators ####
  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/activators"), pattern = "\\.tab$", full.names = TRUE)
  
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
    separate(name, c("Sample", "Day", "Replicate", "Antibody", "Sample_Number"), "_", remove = FALSE) %>%
    mutate(Target = case_when(grepl("T", name) ~ "Zic3 KD", TRUE ~ case_when(grepl("G", name) ~ "Rnf12 KD", TRUE ~ "NTC")),
           dTag = str_match(name, "(?<=[TGL])[0-9]{1,3}"))
  
  result_df$dTag <- factor(result_df$dTag, levels= c("0", "500"))
  result_df$KD <- if_else(result_df$dTag == "500", "ctrl", "KD")
  result_df$KD <- factor(result_df$KD, levels = c("ctrl", "KD"))
  
  # calculate significance
  test_df <- expand.grid(Target = unique(result_df$Target))
  
  test_df$pval <- mapply(function(a,b)
  {t.test(result_df[result_df$Target==a &
                      result_df$KD=="ctrl",]$mean_value,
          result_df[result_df$Target==a &
                      result_df$KD=="KD",]$mean_value,
          paired=TRUE
  )$p.value},
  test_df$Target)
  
  test_df$pval <- vapply(test_df$pval, format_p, character(1))
  
  test_df$KD <- "KD"
  
  
  
  # Create the plot
  plot <- result_df %>%
    ggplot(aes(x = KD, y = mean_value,  ymin=0, group=KD)) +
    geom_point(size=1.5, aes(color = KD, shape= Replicate)) +
    facet_grid(rows=vars(Target), scales = "free") +
    geom_text(data = test_df,
              aes(label = pval, y = 3350), size=6/2.8,
              nudge_x = -0.5, color="black")+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("ctrl" ="#3f414b",KD= "#f37748"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(panel.grid = element_blank()) +
    theme(aspect.ratio = 2)+
    theme(legend.position = "right")
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/activators_H3K9me3_RE57.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  
  # percent change
  # (mean of KD - mean of ctrl ) mean of ctrl x100%
  
  summary_df <- result_df %>%
    group_by(Target, KD) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = KD, values_from = mean_across_reps) %>%
    mutate(
      percent_change = ((KD - ctrl) / ctrl) * 100
    )
  summary_df
  
  
  
  # Rnf12 KO- DKO ####

  file_list <- list.files(path = paste0(wdir, "fig_data/RE57_counts_tab_files/dko/"), pattern = "\\.tab$", full.names = TRUE)
  
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
    separate(name, c("Sample","Condition", "Day", "Antibody", "Replicate", "Sample_Number"), "_") 
  
  result_df$Sample <- if_else(result_df$Sample == "RNF12", "Rnf12 KO", result_df$Sample)
  result_df$Sample <- factor(result_df$Sample, levels = c("WT" , "Rnf12 KO" ,"DKO"))
  
  # stats
  stat.test <- compare_means(
    mean_value ~ Sample,
    data = result_df,
    method = "t.test"
  )

  stat.test$label <- vapply(stat.test$p, format_p, character(1))
  
  stat.test$y.position <- c(4500, 6500, 5500)
  
  
  
  # Create the plot
  
  plot <- result_df %>%
    ggplot(aes(x = Sample, y = mean_value, ymin=0, group=Sample)) +
    geom_point(size=1.5, aes(color = Sample, shape= Replicate)) +
    stat_pvalue_manual(stat.test, label = "label", size = 6/2.8, vjust = -0.2)+
    labs(x = "", y = "Mean norm. counts over RE57") +
    scale_color_manual(values = c("WT" ="#3f414b",`Rnf12 KO`= "#f37748", DKO="#751f58"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "right")
  
  
  fix <- set_panel_size(plot, height = unit(2.4, "cm"), width = unit(1.8, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/DKO_H3K9me3_RE57.pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  
  
  # percent change
  # (mean of KD - mean of WT ) mean of WT x100%
  
  result_df %>%
    group_by(Sample) %>%
    summarize(mean_across_reps = mean(mean_value)) %>%
    pivot_wider(names_from = Sample, values_from = mean_across_reps) %>%
    mutate(
      percent_change_RNF12 = (`Rnf12 KO` / WT) * 100,
      percent_change_DKO = (DKO / WT) * 100
    )
  
  
  