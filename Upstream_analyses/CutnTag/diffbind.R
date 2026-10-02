# diffbind peak analysis
# Nelly Kanata, OWL Schulz
# Created: 14.01.2026
# Modified: 14.01.2026


library(DiffBind)
library(magrittr)
library(dplyr)

setwd("/path/to/your/directory/Analysis/epic2/") # replace this with your working directory
analysis_date <- gsub("-", "", Sys.Date())

# Load data and metadata

dataset <- read.csv("../libraries.csv", sep = ";") # this is the metadata file


#fix colnames to be compatible with what diffbind expects
colnames(dataset) <- c("SampleID","Condition", "Replicate")
  

# location of bam/bed files

dataset$bamReads <- paste0("../bam/", dataset$Condition, "_D4_Rep", 
                           dataset$Replicate, "_S", dataset$SampleID, 
                           ".mapped.sort.blacklisted.bam")

dataset$Peaks <- paste0("../epic2/", dataset$Condition, "_D4_Rep", 
                        dataset$Replicate, "_S", dataset$SampleID,
                           "_epic2.bed")
dataset$PeakCaller <- "bed"

# re-assign SampleIDs because they cannot be duplicated
dataset$SampleID <- seq(1,nrow(dataset))

peaksets.dba <- dba(sampleSheet=dataset, attributes = c(DBA_CONDITION))
peaksets.dba


peaksets.count <- dba.count(peaksets.dba, 
                            minOverlap = 3, 
                            filter=5,
                            bUseSummarizeOverlaps=TRUE, 
                            summits=FALSE, 
                            score=DBA_SCORE_RPKM)
# as i am using summits =FALSE, not all intervals are of the same size. Therefore, I should normalize
# the count of reads per interval size (DBA_SCORE_RPKM)

peaksets.count


# separate analysis for RNF12KO and DKO

for (ko in c("RNF12-KO", "DKO")) {
  
  if (sum(peaksets.count[["samples"]][["Condition"]] == ko) < 3) {
    next
  }
  
  # Subset peaksets.count
  
  subset.peaksets.count <- dba.mask(peaksets.count, attribute = DBA_CONDITION,
                          c("WT", ko), combine='or', bApply=TRUE)
  
  
  peaksets.norm <- dba.normalize(subset.peaksets.count, method=DBA_ALL_METHODS,  
                                 normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

    
  
  # model setup
  peaksets.norm <- dba.contrast(peaksets.norm,  categories=DBA_CONDITION, 
                                reorderMeta=list(Condition="WT"), bNot = FALSE)
  
  peaksets.norm 
  
  # differential analysis
  peaksets.analyzed <- dba.analyze(peaksets.norm, method=DBA_ALL_METHODS)
  peaksets.analyzed
  
  dba.show(peaksets.analyzed, bContrasts = TRUE)
  
  peaksets.DESEQ2 <- dba.report(peaksets.analyzed,
                                contrast=1, th=1, bCounts=TRUE, bNormalized = TRUE,
                                bCalled = TRUE,  method=DBA_DESEQ2)
  peaksets.DESEQ2
     
  
  gr_sorted <- peaksets.DESEQ2[order(mcols(peaksets.DESEQ2)$FDR), ]
  gr_sorted_df <- as.data.frame(gr_sorted)
  
  
  gr_sorted_df$direction <- if_else(gr_sorted_df$Fold > 0.5 & gr_sorted_df$FDR < 0.05, "Up", 
                                    if_else(gr_sorted_df$Fold < -0.5 & gr_sorted_df$FDR < 0.05, "Down","Unchanged"))
  

  # save results 
  
  write.csv(gr_sorted_df, file = paste0("diffbind/", analysis_date, "_",ko, "vsWT_deseq2.csv"),
            row.names = FALSE)

  # save peaks as bedfile 
  
  write_tsv(gr_sorted_df[,c(1:3)], 
            file = paste0("diffbind/", analysis_date, "_",ko, "vsWT_deseq2.bed"), col_names = FALSE
            )
  
  
}


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
