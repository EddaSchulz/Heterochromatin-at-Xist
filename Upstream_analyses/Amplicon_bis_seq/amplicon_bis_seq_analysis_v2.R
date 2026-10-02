## Amplicon Bisulfite seq Analysis
## Nelly Kanata ~ OWL Schulz
## Created: 27.05.2024
## Modified: 05.01.2026

library(tidyverse)
library(magrittr)

analysis_date <- gsub("-", "", Sys.Date())
setwd("./")


# Read methylation data from MOABS mcall

bedfiles <- list.files("bed", full.names = TRUE) # folder "bed" with all bedfiles from MOABS


all_methylation <- do.call(
  rbind,
  lapply(bedfiles, function(x) {
    
    read.table(x, stringsAsFactors = FALSE, skip = 1) %>%
      mutate(filename = str_remove(x, "^bed/"),
             filename = str_remove(filename, "\\.G\\.bed$")) %>%
      separate(
        filename,
        into = c("cell_line", "condition", "replicate", "sample_number"),
        sep = "_",
        remove = TRUE
      )
  })
)


colnames(all_methylation) <- c("chrom",	"start",	"end",	"ratio",	"totalC",	"methC",	"strand",	"next",	
                               "Plus",	"totalC_plus",
                               "methC_plus",	"Minus",	"totalC_minus",	"methC_minus",	"localSeq", "cell_line", "condition", 
                               "replicate", "sample_number")



# assign regions to amplicons
all_methylation$region <- if_else(grepl("FJ10", all_methylation$chrom), "RE57", 
                                  if_else(grepl("FJ28", all_methylation$chrom), "RE58",
                                          "MSR"))

# The coordinates below represent the position of the amplicons in the fasta file used for mapping. 
all_methylation %<>%
  mutate(amplicon = if_else(region == "RE57" & 
                              (start > 400 & start < (400+279)),
                            "RE57_YY1",
                            if_else(all_methylation$region == "RE57" &
                                      (all_methylation$start > (400+279+13) & all_methylation$start < (400+279+13+362)),
                                    "RE57_CGI",
                                    
                                    if_else(region == "RE58" & 
                                              (start > 400 & start < (400+243)),
                                            "RE58_CGI_2",
                                            if_else(all_methylation$region == "RE58" &
                                                      (all_methylation$start > (588) & all_methylation$start < (588+202)),
                                                    "RE58_CGI_1",
                                                    
                                                    "")))))
all_methylation %<>%
  mutate(amplicon = if_else(region == "RE57" & start %in% c(856, 863, 871, 899, 904), "intermediate_RE57",
                            if_else(amplicon == "", "other", amplicon)))

# order sample number
all_methylation$sample_number <- factor(all_methylation$sample_number,
                                   levels = unique(all_methylation$sample_number[order(as.numeric(sub("S", "", all_methylation$sample_number)))]))


write.csv(all_methylation, paste0(analysis_date, "_", "all_methylation_table.csv"))



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

