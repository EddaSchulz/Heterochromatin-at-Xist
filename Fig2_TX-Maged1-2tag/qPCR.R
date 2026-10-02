### qPCR analysis

library(tidyverse)
library(egg)
library(readxl)
library(gridExtra)
library(EnvStats)


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


wd <- "./"
setwd(wd)

## Load data and clean it
raw <- read_xlsx(path = "./fig_data/qPCR5_val5_rawdata.xlsx", sheet = "Results", skip = 25) %>% 
  select("Sample", "Target", "Cq") %>%
  dplyr::rename(CT = "Cq", Gene = "Target") %>%
  filter(Sample != "NTC") %>%
  filter(CT != "UNDETERMINED") %>%
  transform(CT = as.numeric(CT)) %>%
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarize(CT = geoMean(CT)) %>%
  na.omit() %>% 
  filter(Gene != "Arpo-")


#Calculate relative expression to Arpo and Rrm2
analysis <- raw %>% 
  mutate(Sample = sub("^\\d+_", "", Sample)) %>% 
  pivot_wider(names_from = Gene, values_from = CT) %>% 
  mutate(Cont = (Arpo + Rrm2) /2, Arpo = NULL, Rrm2 = NULL) %>%
  select(Sample, Cont, everything()) 

for( i in 3:length(analysis) ) {
  analysis[i] <- analysis[i] - analysis[2]
  analysis[i] <- 2^ - analysis[i]
  analysis[i] <- log2(analysis[i])
}

output <- analysis %>%
  select(-Cont) 

qpcr <- output %>%
  separate(Sample, c("cell_line", "day", "rep"), sep = "_") %>% 
  pivot_longer(-c(cell_line, day, rep), names_to = "gene", values_to = "rel_exp") %>% 
  na.omit()

qpcr$gene <- factor(qpcr$gene, levels = c("Xist", "Nanog", "Esrrb"))
qpcr$cell_line <- factor(qpcr$cell_line , levels = c("tx1072", "maged1"))
qpcr$day <- as.numeric(gsub("d", "", qpcr$day))

#Calculate pvalues
test_df <- expand.grid(day = c(0, 2, 4), gene = unique(qpcr$gene))

test_df$pval <- mapply(function(a, b) {
  t.test(qpcr[qpcr$day==a & qpcr$gene==b & qpcr$cell_line=="maged1",]$rel_exp,
         qpcr[qpcr$day==a & qpcr$gene==b & qpcr$cell_line=="tx1072",]$rel_exp
         )$p.value
}, test_df$day, test_df$gene)

# format pvalue function
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
test_df$cell_line <- "maged1"

#plot relative expression

plot <- qpcr %>% 
  ggplot(aes(x = day, y = rel_exp, fill = cell_line, group = cell_line)) +
  facet_wrap(~gene, scales = "free", ncol=1) +
  geom_point(position = position_dodge(width = 1), shape = 21, stroke = 0.01) +
  stat_summary(fun = "mean", geom = "crossbar", color = "black", position = position_dodge(width = 1),
               lwd = 0.25, width =0.75) +
  scale_fill_manual(values = c("#666666", "#F47748")) +
  geom_text(data = test_df,
            aes(label = pval, y = c(-0.5, -2.5, -1.8,
                                    1, -1, 0, 
                                  -6, 1, 2.2)), size=6/2.8,
            nudge_x = 0, color="black") +
  ylab(expression("Rel. expression  (log"[2]*")")) +
  xlab("Differentiation timepoint (days)") +
  scale_x_continuous( breaks = seq(0, 4, by = 2), expand = expansion(mult = c(0.1, 0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
  

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(fix)

ggsave("./Figures/TX-Maged1-2tag_qPCR.pdf", fix, dpi = 300,
       useDingbats=FALSE)
