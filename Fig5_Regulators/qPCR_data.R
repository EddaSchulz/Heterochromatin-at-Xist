## qPCR Analysis
## Nelly Kanata ~ OWL Schulz
## Created: 26.10.2025
## Modified: 26.10.2025

library(tidyverse)
library(magrittr)
library(ggpubr)
library(egg)

setwd("./")


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  plot.title = element_text(size = 8),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

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

## repressors ####

# load data
zfp280c_zfp36l1 <- read.csv("./fig_data/qPCR_Zfp280c_Zfp36l1_2.csv")
rex1_zfp281 <- read.csv("./fig_data/qPCR_Rex1_Zfp281.csv")

# pre-processing
zfp280c_zfp36l1 %<>%
  select(-Sample) %>%
  mutate(Cell_line = if_else(Cell_line == "WT", "ctrl", "KD"))

rex1_zfp281 %<>%
  select(-Sample) %>%
  mutate(Cell_line = if_else(KD == "NTC", "ctrl", "KD"))
  

# combine data
rel.data_long <- bind_rows(zfp280c_zfp36l1, rex1_zfp281)

rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix","Dnmt3b", "Fgf5",  "Otx2", "Setdb1",
                                                             "Zfp280c", "Zfp36l1", "Zfp281",  "Rex1",   "Myc",  "Nanog" ,
                                                             "Esrrb",   "Prdm14" , "Oct4",
                                                          "Rnf12","Zic3"))

rel.data_long$Cell_line <- factor(rel.data_long$Cell_line, levels=c("ctrl", "KD"))

# Xist and Tsix relative expression
# Castuner plot
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix")) %>%
  filter(!KD %in% c("Myc","Setdb1")) %>%
  filter(KD %in% c("Zfp280c", "Zfp36l1", "Rnf12", "Zic3")) %>%
  filter( Day != "D0" )


Zfp281 <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix")) %>%
  filter(!KD %in% c("Myc","Setdb1")) %>%
  filter(KD %in% c("NTC", "Zfp281")) %>%
  filter(Day != "D0")

Zfp281$KD <- "Zfp281"

Rex1 <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix")) %>%
  filter(!KD %in% c("Myc","Setdb1")) %>%
  filter(KD %in% c("NTC", "Rex1")) %>%
  filter(Day != "D0")

Rex1$KD <- "Rex1"

filtered_rel.data_long <- rbind(filtered_rel.data_long, Rex1, Zfp281)

filtered_rel.data_long$KD <- paste0(filtered_rel.data_long$KD, " KD")

filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c("Rnf12 KD", "Zic3 KD", 
                                                                        "Zfp280c KD", "Zfp36l1 KD",
                                                                        "Rex1 KD", "Zfp281 KD") )

### effect on Xist/Tsix #####
test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD),
                         gene = unique(filtered_rel.data_long$gene),
                       Day= unique(filtered_rel.data_long$Day))

test_df$pval <- mapply(function(a,b,c)
{t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                 filtered_rel.data_long$gene==b &
                                 filtered_rel.data_long$Day==c &
                                 filtered_rel.data_long$Cell_line=="ctrl",]$log2_re,
        filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                 filtered_rel.data_long$gene==b &
                                 filtered_rel.data_long$Day==c &
                                 filtered_rel.data_long$Cell_line=="KD",]$log2_re,
       paired=FALSE
      # var.equal = TRUE
        )$p.value},
test_df$KD, test_df$gene, test_df$Day)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$Cell_line <- "KD"

# significance position
label_pos <- filtered_rel.data_long %>%
  group_by(gene, KD, Cell_line) %>%
  summarise(y = max(log2_re, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y + 1)   # small offset above the highest point


test_df %<>%
  left_join(label_pos, by = c("gene", "KD", "Cell_line"))


plot <- filtered_rel.data_long %>%
  filter(KD %in% c("Zfp280c KD", "Zfp36l1 KD", "Zfp281 KD", "Rex1 KD")) %>%
  filter(gene=="Xist") %>%
  ggplot(aes(x=Cell_line, y=log2_re, color = Cell_line)) +
  facet_grid(cols=vars(Day), rows=vars(KD), scales = "free") +
  geom_point(position = position_dodge(width=0.75), shape=16) +
  xlab("") + ylab(expression("Xist Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = KD), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("ctrl"= "#3f414b", "KD"= "#3b9ad9"), name = '') +
  geom_text(data = test_df[test_df$gene == "Xist",],
            aes(label = pval, y = y), size=6/2.8,
            nudge_x = -0.5, color="black")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme(legend.position = "none")


fix <- set_panel_size(plot, height = unit(1.8, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_repressorKDs_D4.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


### KD efficiencies ####

filtered_rel.data_long <- rel.data_long %>%  
  filter(gene == KD) %>%
  mutate(perc_of_NTC = if_else(!is.na(re_NTC), re/re_NTC*100,
                               re_minus/re_plus*100)) %>%
  filter(!KD %in% c("Myc","Setdb1")) %>%
  filter(Day == "D4" | is.na(Day) | Day == "D2") %>%
  filter(Cell_line == "KD")

filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c("Zfp280c", "Zfp36l1", "Rex1", "Zfp281", "Rnf12", "Zic3"))


# t-tests
test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD),
                       Day = unique(filtered_rel.data_long$Day))

test_df$pval <- mapply(function(a, b)
{
  if (is.na(filtered_rel.data_long$re_NTC[filtered_rel.data_long$KD==a])[1]) {
    
    t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$Day==b
                                  ,]$re_minus,
           filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$Day==b,]$re_plus
    )$p.value
  } else {
    t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$Day==b,]$re,
           filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$Day==b,]$re_NTC

    )$p.value
    
          }
  },
test_df$KD, test_df$Day)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$Cell_line <- "KD"


plot <- filtered_rel.data_long %>%
  filter(KD %in% c("Zfp280c", "Zfp281", "Zfp36l1", "Rex1")) %>%
  ggplot(aes(x=KD, y=perc_of_NTC)) + 
  geom_point(aes(color=Day),
             position = position_dodge(width=0.9), shape=16) +

  stat_summary(aes(group = Day), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, 
               position = position_dodge(width = 0.9),
               color = "black", show.legend = FALSE)+
  
  geom_hline(yintercept =100, linetype="dashed", color="gray") +
  scale_color_manual(values=c("D4"= "black",
                              "D2"= "#a5c4ab"), name = '')+
  geom_text(data = test_df %>% filter(KD %in% c("Zfp280c", "Zfp281", "Zfp36l1", "Rex1")), 
            aes(label = pval, 
            y = c(15,45, 60,15, 15, 45,55, 15),
            group = Day), 
            size=6/2.8, angle=90, hjust=-0.1, vjust=0.5,
            position = position_dodge(width = 0.9), show.legend = FALSE)+
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, by = 20)) +
  xlab("Knock-down") + ylab(expression("% control expression")) +
  theme(legend.position = "none",
        )


fix <- set_panel_size(plot, height = unit(2.1, "cm"), width = unit(3, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_KD_repressors_efficiencies_D2_D4.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



### effect on differentiation #####

filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Nanog", "Esrrb", "Fgf5", "Otx2")) %>%
  mutate(log2FC_NTC = log2(if_else(!is.na(re_NTC), re/re_NTC,
                               re_minus/re_plus))) %>%
  filter(!KD %in% c("Myc","Setdb1", "NTC")) %>%
  filter(!KD %in% c("Rnf12", "Zic3")) %>%
  filter(Cell_line == "KD") %>%
  filter(Day != "D0")

filtered_rel.data_long$Day <- factor(filtered_rel.data_long$Day, levels=c( "D2", "D4"))

filtered_rel.data_long$KD <- paste0(filtered_rel.data_long$KD, " KD")
filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c("Zfp280c KD", "Zfp36l1 KD", "Rex1 KD", "Zfp281 KD"))


# t-tests
test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD),
                       gene= unique(filtered_rel.data_long$gene),
                       Day= unique(filtered_rel.data_long$Day))

test_df$pval <- mapply(function(a, b,c)
{
  if (is.na(filtered_rel.data_long$re_NTC[filtered_rel.data_long$KD==a])[1]) {
    
    t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                  filtered_rel.data_long$gene==b &
                                  filtered_rel.data_long$Day==c,]$re_minus,
           filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$gene==b &
                                    filtered_rel.data_long$Day==c,]$re_plus,
           paired=FALSE
    )$p.value
  } else {
    t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$gene==b &
                                    filtered_rel.data_long$Day==c,]$re,
           filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                    filtered_rel.data_long$gene==b &
                                    filtered_rel.data_long$Day==c,]$re_NTC,
           
           paired=FALSE
           
    )$p.value
    
  }
},
test_df$KD, test_df$gene, test_df$Day)

test_df$pval <- vapply(test_df$pval, format_p, character(1))


plot <- filtered_rel.data_long %>%
  ggplot(aes(x=Day, y=log2FC_NTC)) + 
  facet_grid(cols=vars(gene), rows=vars(KD), scales="free") +
  geom_hline(yintercept =0, linetype="dashed", color="gray")+
  geom_point( position = position_dodge(width=0.75), shape=16, color="#3b9ad9") +
  geom_text(data = test_df, aes(label = pval, y = 0.9, group = Day), 
    size=6/2.8, angle=90, hjust=-0.1, vjust=0.5,
    position = position_dodge(width = 0.9), show.legend = FALSE)+
  
  stat_summary(aes(group = KD), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, 
               position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE)+
  scale_y_continuous(limits = c(-2.3,4.6), breaks = seq(-2, 2, by = 2)) +
  xlab("Differentiation Timepoint") + ylab(expression("log"[2]*"FC")) 

fix <- set_panel_size(plot, height = unit(1.8, "cm"), width = unit(1.2, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_repressors_KD_differentiation.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


## activators ####

# load data
rnf12 <- read.csv("./fig_data/qPCR_Rnf12_KD.csv")
zic3 <- read.csv("./fig_data/qPCR_Zic3.csv") 

# pre-processing
rnf12$KD <- "Rnf12"
colnames(rnf12)[1] <- "Biological_Replicate"

zic3 %<>%
  select(-ID)
colnames(zic3)[1] <- "Biological_Replicate"

# combine data
rel.data_long <- bind_rows(rnf12, zic3)

rel.data_long$gene <- factor(rel.data_long$gene, levels=c("Xist", "Tsix","Dnmt3b", "Fgf5",  "Otx2", "Setdb1",
                                                          "Zfp280c", "Zfp36l1", "Zfp281",  "Rex1",   "Myc",  "Nanog" ,
                                                          "Esrrb",   "Prdm14" , "Oct4",
                                                          "Rnf12","Zic3"))

rel.data_long$Cell_line <- factor(rel.data_long$Cell_line, levels=c("ctrl", "KD"))

# Xist and Tsix relative expression
filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Xist", "Tsix")) 

filtered_rel.data_long$KD <- paste0(filtered_rel.data_long$KD, " KD")

filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c("Rnf12 KD", "Zic3 KD") )

### effect on Xist/Tsix #####
filtered_rel.data_long %<>%
  arrange(Biological_Replicate)

test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD),
                       gene = unique(filtered_rel.data_long$gene))

test_df$pval <- mapply(function(a,b)
{t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                 filtered_rel.data_long$gene==b &
                                 filtered_rel.data_long$Cell_line=="ctrl",]$log2_re,
        filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                 filtered_rel.data_long$gene==b &
                                 filtered_rel.data_long$Cell_line=="KD",]$log2_re,
        paired=FALSE
        # var.equal = TRUE
)$p.value},
test_df$KD, test_df$gene)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$Cell_line <- "KD"

plot <- filtered_rel.data_long %>%
  filter(gene=="Xist") %>%
  ggplot(aes(x=Cell_line, y=log2_re, color = Cell_line)) +
  facet_grid(cols=vars(gene), rows=vars(KD), scales = "free") +
  geom_point(position = position_dodge(width=0.75), shape=16) +
  xlab("") + ylab(expression("Xist Rel. expression  (log"[2]*")")) +
  stat_summary(aes(group = KD), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE) +
  scale_color_manual(values=c("ctrl"= "#3f414b", "KD"= "#f37748"), name = '') +
  geom_text(data = test_df[test_df$gene == "Xist",],
            aes(label = pval, y = -2.5), size=6/2.8,
            nudge_x = -0.5, color="black")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme(legend.position = "none")


fix <- set_panel_size(plot, height = unit(1.8, "cm"), width = unit(1, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_Xist_activatorsKDs.pdf"), fix,
       dpi = 300, useDingbats=FALSE)


### KD efficiencies ####

filtered_rel.data_long <- rel.data_long %>%  
  filter(gene == KD) %>%
  mutate(perc_of_NTC = re_minus/re_plus*100) %>%
  filter(Cell_line == "KD")

filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c( "Rnf12", "Zic3"))


# t-tests
test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD))

test_df$pval <- mapply(function(a)
{
  
  t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a 
                                ,]$re_minus,
         filtered_rel.data_long[filtered_rel.data_long$KD==a 
                                ,]$re_plus
  )$p.value
  
},
test_df$KD)

test_df$pval <- vapply(test_df$pval, format_p, character(1))

test_df$Cell_line <- "KD"


plot <- filtered_rel.data_long %>%
  ggplot(aes(x=KD, y=perc_of_NTC)) + 
  geom_point(position = position_dodge(width=0.9), shape=16) +
  stat_summary( geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, 
                position = position_dodge(width = 0.9),
                color = "black", show.legend = FALSE)+
  
  geom_hline(yintercept =100, linetype="dashed", color="gray") +
  geom_text(data = test_df, 
            aes(label = pval, 
                y = 30), 
            size=6/2.8, angle=90, hjust=-0.1, vjust=0.5,
            position = position_dodge(width = 0.9), show.legend = FALSE)+
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, by = 20)) +
  xlab("Knock-down") + ylab(expression("% control expression")) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


fix <- set_panel_size(plot, height = unit(2.1, "cm"), width = unit(1.3, "cm"))
grid.arrange(fix)

ggsave(paste0("./Figures/qPCR_KD_activators_efficiencies.pdf"), fix,
       dpi = 300, useDingbats=FALSE)



### effect on differentiation #####

filtered_rel.data_long <- rel.data_long %>%  
  filter(gene %in% c("Nanog", "Esrrb", "Fgf5", "Otx2")) %>%
  mutate(log2FC_NTC = log2(re_minus/re_plus)) %>%
  filter(Cell_line == "KD") 


filtered_rel.data_long$KD <- paste0(filtered_rel.data_long$KD, " KD")
filtered_rel.data_long$KD <- factor(filtered_rel.data_long$KD, levels=c("Rnf12 KD", "Zic3 KD"))

# t-tests
test_df <- expand.grid(KD = unique(filtered_rel.data_long$KD),
                       gene= unique(filtered_rel.data_long$gene))

test_df$pval <- mapply(function(a, b)
{
  t.test(filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                  filtered_rel.data_long$gene==b ,]$re_minus,
         filtered_rel.data_long[filtered_rel.data_long$KD==a &
                                  filtered_rel.data_long$gene==b ,]$re_plus,
         paired=FALSE
  )$p.value
},
test_df$KD, test_df$gene)

test_df$pval <- vapply(test_df$pval, format_p, character(1))
test_df$Cell_line <- "KD"

plot <- filtered_rel.data_long %>%
  ggplot(aes(x=Cell_line, y=log2FC_NTC)) + 
  facet_grid(cols=vars(gene), rows=vars(KD), scales="free") +
  geom_hline(yintercept =0, linetype="dashed", color="gray")+
  geom_point( position = position_dodge(width=0.75), shape=16, color="#f37748") +
  geom_text(data = test_df, aes(label = pval, y = 1.6), 
            size=6/2.8, 
            position = position_dodge(width = 0.9), show.legend = FALSE)+
  
  stat_summary(aes(group = KD), geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, 
               position = position_dodge(width = 0.75),
               color = "black", show.legend = FALSE)+
  scale_y_continuous(limits = c(-1.6,1.8), breaks = seq(-2, 2, by = 1)) +
  xlab("Differentiation Timepoint") + ylab(expression("log"[2]*"FC")) 

fix <- set_panel_size(plot, height = unit(1.8, "cm"), width = unit(1, "cm"))
grid.arrange(fix)


ggsave(paste0("./Figures/qPCR_activators_KD_differentiation.pdf"), fix,
       dpi = 300, useDingbats=FALSE)
