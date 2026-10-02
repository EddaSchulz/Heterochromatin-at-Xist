## DNA methylation WGBS vs amplicon
## Nelly Kanata ~ OWL Schulz
## Created: 21.05.2026
## Modified: 21.05.2026

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  #axis.text.x = element_text(angle = 45, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


setwd("./")

# load data

timecourse <- read.csv("./fig_data/timecourse_all_methylation_table.csv")

filtered_df <- timecourse %>%
  filter(region %in% c("RE57")) %>%
  filter(amplicon != "other", amplicon !="intermediate_RE57") 

filtered_df$amplicon <- factor(filtered_df$amplicon, levels=c("RE57_YY1","RE57_CGI"))
filtered_df$timepoint <- factor(filtered_df$timepoint)
filtered_df$timepoint <- as.numeric(gsub("T", "", filtered_df$timepoint))


### mean methylation plot ####
RE57_amplicons <- c("RE57_YY1", "RE57_CGI")

methC <- filtered_df %>%
  filter(amplicon %in% RE57_amplicons) %>% # I combine RE57 amplicons in one metric
  filter(totalC > 10) %>% # coverage cutoff
  group_by(sample_number, amplicon) %>%
  mutate(totalC_mean = round(mean(totalC), 0)) %>%
  select(start, ratio, cell_line, amplicon,timepoint,totalC_mean,
         replicate, sample_number) %>%
  pivot_wider(names_from = start, values_from = ratio, values_fill = NA) %>%
  arrange(cell_line)
methC <- as.data.frame(methC)

methC %<>%
  rowwise() %>%
  mutate(mean_meth = mean(c_across(7:ncol(methC)), na.rm = TRUE))


methC$mean_meth <- methC$mean_meth *100




# read WGBS analysis data
wgbs <- list.files("fig_data/WGBS_quantified", full.names = TRUE)

#(these are already average of two replicates)
wgbs_data_all <- map_dfr(wgbs, function(x) {
  
  fname <- basename(x)
  fname <- tools::file_path_sans_ext(fname)
  
  parts <- strsplit(fname, "_")[[1]]
  
  cell_time <- strsplit(parts[3], "-")[[1]]
  
  read.table(x, stringsAsFactors = FALSE, skip = 1) %>%
    mutate(
      cell_line = cell_time[1],
      amplicon  = parts[2],
      timepoint = tail(cell_time, 1), 
      
    )
})


wgbs_data_all$amplicon <- gsub("-", "_", wgbs_data_all$amplicon)

for (amplicon_i in RE57_amplicons) {
  
  
  # average across replicates
  methC_selected <- methC %>%
    filter(amplicon == amplicon_i) %>%
    
    select(-sample_number, -replicate, -amplicon, -totalC_mean) %>%
    group_by(cell_line, timepoint) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    select(where(~ any(!is.na(.)))) # to exclude columns from other amplicons
  
  
  # filter out CpGs that are not sequenced (middle of amplicon)
  if (amplicon_i == "RE57_CGI") {
    coordinates_to_exclude <- c(103482178, 103482243)
    wgbs_data <- wgbs_data_all %>%
      filter(amplicon==amplicon_i) %>%
      filter(V2 < coordinates_to_exclude[1] | V2 > coordinates_to_exclude[2])
  } else {
    wgbs_data <- wgbs_data_all %>%
      filter(amplicon==amplicon_i)
  }
  
  
  # assign amplicon CpG position to coordinates: !doublecheck!
  
  # this is to confirm that we have the same number of positions in both experiments.
  # There are 4 non-position columns in methC_selected
 length(colnames(methC_selected))-3 == length(unique(wgbs_data$V2))
  
 
  offset <- unique(unique(wgbs_data$V2) - as.numeric(colnames(methC_selected[, grepl("^[0-9]", names(methC_selected))])))
 
   wgbs_data$Position <- wgbs_data$V2 - offset

     wgbs_data %<>%
    select(V4, cell_line, timepoint, Position) %>%
    pivot_wider(names_from = Position, values_from = V4, values_fill = NA) %>%
    arrange(timepoint, cell_line)
  
  #mean methylation across all positions
  wgbs_data %<>%
    rowwise() %>%
    mutate(mean_meth = mean(c_across(3:ncol(wgbs_data)), na.rm = TRUE))
  
  methC_selected$experiment <- "Amplicon Seq"
  wgbs_data$experiment <- "WGBS"
  
  wgbs_data$timepoint <- as.numeric(gsub("T", "", wgbs_data$timepoint))
  wgbs_data$cell_line <- if_else(wgbs_data$cell_line == "XOB7", "XO", "XX")
  

  
  
  merged_df <- rbind(methC_selected, wgbs_data)
  merged_df$experiment <- factor(merged_df$experiment, levels=c("WGBS", "Amplicon Seq"))
  
  merged_df$Day <- merged_df$timepoint / 24
  
  plot <- merged_df%>%
       filter(timepoint %in% c(0, 48, 96)) %>%
       filter(cell_line == "XO") %>%
  ggplot( aes(Day, mean_meth)) +
    facet_grid(cols = vars(cell_line)) +
    geom_point(aes(colour = experiment, shape=experiment), size = 1.5, alpha = 0.9)+
    scale_color_manual(values=c(`Amplicon Seq`="black",WGBS="#f37748")) +
    scale_y_continuous( limits = c(0,110),breaks = seq(0, 100, by = 25)) +
    scale_x_continuous(
      breaks = c(0,2, 4))+
    ylab("DNA methylation (%)")+
    xlab("Differentiation timepoint (days)")
    ggtitle(amplicon_i)
  
  
  
  fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
  grid.arrange(fix)
  
  ggsave(paste0("./Figures/WGBS_vs_amplicon_", amplicon_i, ".pdf"), fix,
         dpi = 300, useDingbats=FALSE)
  
  
}

