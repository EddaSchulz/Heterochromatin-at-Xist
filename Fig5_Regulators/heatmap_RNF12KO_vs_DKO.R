# heatmap
# Nelly Kanata, OWL Schulz
# Created: 14.01.2026
# Modified: 16.02.2026


library(GenomicRanges)
library(ComplexHeatmap)
library(rtracklayer)
library(dplyr)
library(grid)

setwd("./")


# Read data
Rnf12KO <- read.csv(file = 
                      "fig_data/deseq2/RNF12vsWT_deseq2.csv")
Rnf12KO_GR <- GRanges(
  seqnames = Rnf12KO[[1]],
  ranges   = IRanges(start = Rnf12KO[[2]], end = Rnf12KO[[3]]),
  strand   = "*",
  mcols    = Rnf12KO[, -(1:3)]
)

# Read bw tracks
RNF12_vs_WT_log2FC <- BigWigFile("") # INSERT PATH TO THE BIGWIG FILES (eg download from GEO)
DKO_vs_WT_log2FC <- BigWigFile("") # SAME

# get mean signal in defined peaks
mean_signal <- summary(DKO_vs_WT_log2FC, Rnf12KO_GR, type="mean")
Rnf12KO_GR$DKO_vs_WT_log2FC <- sapply(mean_signal, function(x) x$score)

mean_signal <- summary(RNF12_vs_WT_log2FC, Rnf12KO_GR, type="mean")
Rnf12KO_GR$RNF12_vs_WT_log2FC <- sapply(mean_signal, function(x) x$score)

Rnf12KO_GR_df <- as.data.frame(Rnf12KO_GR)


Rnf12KO_GR_df$mcols.direction <- factor(Rnf12KO_GR_df$mcols.direction, levels = c("Up", "Unchanged", "Down"))

Rnf12KO_GR_df <- Rnf12KO_GR_df[order(Rnf12KO_GR_df$mcols.direction, Rnf12KO_GR_df$mcols.FDR),]


Rnf12KO_GR_df$peakID <- rownames(Rnf12KO_GR_df)
Rnf12KO_GR_df$chrX <- factor(if_else(Rnf12KO_GR_df$seqnames == "chrX", "chrX", "autosome"))


heatmap_matrix <- as.matrix(Rnf12KO_GR_df[,c("RNF12_vs_WT_log2FC", "DKO_vs_WT_log2FC")])
rownames(heatmap_matrix) <- Rnf12KO_GR_df$peakID

heatmap_annot <- Rnf12KO_GR_df[,c("chrX", "mcols.direction")]


gaps <- cumsum(table(heatmap_annot$mcols.direction))

table(heatmap_annot$mcols.direction)
heatmap_annot$mcols.direction <- factor(
  heatmap_annot$mcols.direction,
  levels = c("Up", "Unchanged", "Down"),
  labels = c(paste0("Up (n=", gaps[1], ")"), paste0("Unchanged (n=",
                                                    gaps[2]-gaps[1], ")"), paste0("Down (n=", gaps[3]-gaps[2], ")"))
)
colnames(heatmap_annot) <- gsub("mcols.", "", colnames(heatmap_annot))

annotation_colors <- list(
  direction = c(
    "Up (n=547)" = "#f37748",
    "Unchanged (n=9357)" = "#ebebeb",
    "Down (n=835)" = "#107bc0"))

colnames(heatmap_matrix) <- gsub("_", " ", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub("RNF12 ", "RNF12 KO ", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub(" vs WT log2FC", "", colnames(heatmap_matrix))


# set 99th percentile as limit to exclude outliers in coloring the heatmap
thresh99 <- quantile(as.vector(heatmap_matrix), 0.99, na.rm = TRUE)


pdf(file= paste0("Figures/RNF12KO-DKO_heatmap_log2fc.bw.pdf"))

Heatmap(
  heatmap_matrix,
  name="Average log2FC\n(rel. to WT) H3K9me3\nacross peak",
  col = circlize::colorRamp2(
    seq(-thresh99, thresh99, length.out = 100),
    colorRampPalette(c("#107bc0", "white", "#f37748"))(100)
  ),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  left_annotation = rowAnnotation(
    df = heatmap_annot[, 2, drop = FALSE],
    col = annotation_colors, show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm")
    ),
  width = unit(1.75, "cm"),
  height = unit(3.5, "cm"),
  row_split = c(factor(rep(c( "Up\n(n=547)", "Unchanged\n(n=9357)", "Down\n(n=835)"), c(547, 9357, 835)),
                       levels=c("Up\n(n=547)", "Unchanged\n(n=9357)", "Down\n(n=835)"))),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  )
  
)

dev.off()