## epiallele Analysis
## Nelly Kanata ~ OWL Schulz
## Created: 04.02.2026
## Modified: 04.02.2026


library("epialleleR")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(magrittr)
library(egg)

setwd("./Fig1_heterochromatin_dynamics")


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  plot.title = element_text(size = 8),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

# Read data
bam_dir <- "../Upstream_analyses/Amplicon_bis_seq/bam" # this is the directory where the bam files from bsmapz are


# define amplicon regions
bed_YY1 <- GRanges(
  seqnames = Rle("FJ14-FJ15.FJ10-FJ11::chrX:103481327-103482781"),
  ranges = IRanges(start =400, end = 679))

bed_CpG <- GRanges(
  seqnames = Rle("FJ14-FJ15.FJ10-FJ11::chrX:103481327-103482781"),
  ranges = IRanges(start =692, end = 1054))

bed_RE58_1 <- GRanges(
  seqnames = Rle("FJ28-FJ29.FJ26-FJ27::chrX:103482701-103483891"),
  ranges = IRanges(start =588, end = 790))

bed_RE58_2 <- GRanges(
  seqnames = Rle("FJ28-FJ29.FJ26-FJ27::chrX:103482701-103483891"),
  ranges = IRanges(start =400, end = 643))

# sample reference genome
genome <- preprocessGenome("../Upstream_analyses/Amplicon_bis_seq/ref/mm10_RE57-58.fa")



beta_list_YY1 <- list()
beta_list_CGI <- list()
for (rep in c("Rep1", "Rep2", "Rep3")) {
  
  file_list <- list.files(
    path = bam_dir,
    pattern = rep,
    full.names = TRUE
  )
  
  for (chosen_file in file_list) {
    
    sample <- stringr::str_extract(basename(chosen_file), "^[^_]+_[^_]+")
    
    
    input.bam <- chosen_file
    # as the bam files were created with bsmapz, they have to be modified
    # resulting BAM with XG/XM tags
    output.bam <- tempfile(pattern="output-", fileext=".bam")
    
    # calls cytosine methylation and stores it in the output BAM
    callMethylation(input.bam, output.bam, genome)
    
    # First, let's extract base methylation information for sequencing reads
    patterns_YY1 <- extractPatterns(
      bam=output.bam,
      bed=bed_YY1)
    nrow(patterns_YY1)
    
    
    patterns_CpG <- extractPatterns(
      bam=output.bam,
      bed=bed_CpG)
    nrow(patterns_CpG)
    
    
    beta_list_YY1[[length(beta_list_YY1) + 1]] <- tibble(
      beta = patterns_YY1$beta,
      sample = sample,
      replicate = rep)
    
    beta_list_CGI[[length(beta_list_CGI) + 1]] <- tibble(
      beta = patterns_CpG$beta,
      sample = sample,
      replicate = rep)
    
  }
}

combined_beta_list_YY1 <- dplyr::bind_rows(beta_list_YY1)
combined_beta_list_CGI <- dplyr::bind_rows(beta_list_CGI)


combined_beta_list_YY1 %<>%
  separate(sample, c("cell_line", "timepoint"), sep="_")
combined_beta_list_CGI %<>%
  separate(sample, c("cell_line", "timepoint"), sep="_")


plot <- ggplot(combined_beta_list_YY1, aes(x = timepoint, y = beta, fill = timepoint)) +
  facet_grid(rows = vars(cell_line)) +
  
  geom_violin(trim = TRUE, color = NA, scale="area") +
  scale_fill_manual(values = colorRampPalette(
    c("#a3c9a8", "#699d73", "#2d6a4f", "#1b4332")
  )(length(unique(combined_beta_list_YY1$timepoint)))) +
  scale_y_continuous(
    limits = c(-0, 1),
    breaks = seq(0, 1, by = 0.25)
  )+
  xlab("Differentiation timepoint (days)") +
  ylab("Average DNA methylation per read (%)") +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  ggtitle("YY1 amplicon")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3.5, "cm"))
grid.arrange(fix)

ggsave(paste0("Figures/YY1_epialleles_timecourse.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



plot <- ggplot(combined_beta_list_CGI, aes(x = timepoint, y = beta, fill = timepoint)) +
  facet_grid(rows = vars(cell_line)) +
  
  geom_violin(trim = TRUE, color = NA, scale="area") +
  scale_fill_manual(values = colorRampPalette(
    c("#a3c9a8", "#699d73", "#2d6a4f", "#1b4332")
  )(length(unique(combined_beta_list_CGI$timepoint)))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  )+
  xlab("Differentiation timepoint (days)") +
  ylab("Average DNA methylation per read (%)") +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  ggtitle("RE57_CGI amplicon")

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3.5, "cm"))
grid.arrange(fix)

ggsave(paste0("Figures/CGI_epialleles_timecourse.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
