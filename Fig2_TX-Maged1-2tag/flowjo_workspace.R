## importing gating from flowjo and plotting
## Nelly Kanata ~ OWL Schulz
## 01.05.2026

library(CytoML)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(egg)

# theme
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6),
                  axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

setwd("./")

# mESCs ####

# read data
wsfile <- "2026-02-02_sorting-2i-clones.wsp"
ws <- open_flowjo_xml(wsfile)

# extract groups
fj_ws_get_sample_groups(ws)
# extract samples
fj_ws_get_samples(ws)

gs <- flowjo_to_gatingset(ws, name="Maged1 E6", subset = "Day0_Maged1_E6_004.fcs", path = "./")
flowWorkspace::plot(gs)

 p <- autoplot(gs, gate = c("1","2", "3", "4"), bins = 256) +
  ggcyto_par_set(limits=list(x=c(0,1.3), y=c(-0.45,1.3)))+ # works as scaling factor!
   ylab("mScarlet-I3 (615/20) (a.u.)") + xlab("mStayGold (510/20) (a.u.)")+labs(title = NULL)+
   theme(legend.position = "right", legend.key.height =unit(0.2,"cm"), legend.key.width =unit(0.4,"cm"))

 p
 
 ggsave(paste0("./Figures/d0_scatterplot_autoplot.pdf"), p ,height = unit(2, "cm"), width = unit(2.5, "cm"),
        dpi = 300, useDingbats=FALSE)
 
 
 # sorted day 4 ####
 
 # read data
 dataDir <- system.file("./")
 wsfile <- "2026-03-18_cutntag_maged1_d4.wsp"
 ws <- open_flowjo_xml(wsfile)
 
 # extract groups
 fj_ws_get_sample_groups(ws)
 # extract samples
 fj_ws_get_samples(ws)
 
 gs <- flowjo_to_gatingset(ws, name="rep1", subset = "Day4_Maged1_E6_001.fcs", path = "./")
 flowWorkspace::plot(gs)
 
 
 # To see the names of the gates
 colnames(gs[[1]])
 
 p <- autoplot(gs, gate = c("1", "3", "4"), bins = 256) +
   ggcyto_par_set(limits=list(x=c(0,0.9), y=c(0.04,0.9)))+ # works as scaling factor!
   ylab("mScarlet-I3 (615/20) (a.u.)") + xlab("mStayGold (510/20) (a.u.)")+labs(title = NULL)+
   theme(legend.position = "right", legend.key.height =unit(0.2,"cm"), legend.key.width =unit(0.4,"cm"))
 
 p
 
 ggsave(paste0("./Figures/rep1_scatterplot_autoplot.pdf"), p ,height = unit(2, "cm"), width = unit(2.5, "cm"),
        dpi = 300, useDingbats=FALSE)
 
 
 # upstream gates ####
 
 #### Live-dead ####
 
 # SSC-A vs FSC-A
 # check autoplot first
 autoplot(gs, gate = c("Cells"), bins = 64) +
   ggcyto_par_set(limits="data")+
   theme_classic()
 
 
 labels_df <- data.frame(
   pop = c("Cells"),
   label = c("Live cells"),
   x = c(1.8e05), 
   y = c(0.3e05)
 )
 
 p <- ggcyto(gs, aes(x = `FSC-A`, y = `SSC-A`),
             subset = "root") + 
   
   geom_hex(bins = 128) +   
   geom_gate(c("Cells"), colour="red", size=0.65) +
   geom_stats(adjust = 0.05, size =6/2.8) +  
   ggcyto_par_set(limits = "data") +
   geom_label(data = labels_df, aes(x = x, y = y, label = label), fill = c("white"),
              size=6/2.8) +
   labs(title = NULL)+
   theme(legend.position = "right", legend.key.height =unit(0.2,"cm"), legend.key.width =unit(0.4,"cm"))
 
 p
 
 ggsave(paste0("./Figures/rep1_SSC-A_FSC-A.pdf"), p ,height = unit(2, "cm"), width = unit(2.5, "cm"),
        dpi = 300, useDingbats=FALSE)

 
 #### FSC Singlets ####
 # FSC-W vs FSC-H
 autoplot(gs, gate = c("P2"), bins = 64) +
   ggcyto_par_set(limits="data")+
   theme_classic()
 
 
 labels_df <- data.frame(
   pop = c("P2"),
   label = c("P2"),
   x = c(200000), 
   y = c(85000)
 )
 
 p <- ggcyto(gs, aes(x = `FSC-H`, y = `FSC-W`),
             subset = "Cells") + 
   
   geom_hex(bins = 128) +
   geom_gate(c("P2"), colour="red", size=0.65) +
   geom_stats(adjust = 0.05, size =6/2.8) +  
   ggcyto_par_set(limits = "data") +
   geom_label(data = labels_df, aes(x = x, y = y, label = label), fill = c("white"),
              size=6/2.8) +
   labs(title = NULL)+
   theme(legend.position = "right", legend.key.height =unit(0.2,"cm"), legend.key.width =unit(0.4,"cm"))
 
 p
 
 ggsave(paste0("./Figures/rep1_FSC_singlets.pdf"), p ,height = unit(2, "cm"), width = unit(2.5, "cm"),
        dpi = 300, useDingbats=FALSE)
 
 #### SSC Singlets ####
 
 # SSC-W vs SSC-H
 autoplot(gs, gate = c("Singlets"), bins = 64) +
   ggcyto_par_set(limits="data")+
theme_classic()
 
 
 labels_df <- data.frame(
   pop = c("Singlets"),
   label = c("Singlets"),
   x = c(150000),
   y = c(120000)
 )
 
 p <- ggcyto(gs, aes(x = `SSC-H`, y = `SSC-W`),
             subset = "P2") + 
   
   geom_hex(bins = 128) +   
   geom_gate(c("Singlets"), colour="red", size=0.65) +
   geom_stats(adjust = 0.05, size =6/2.8) + 
  ggcyto_par_set(limits = "data") +
   geom_label(data = labels_df, aes(x = x, y = y, label = label), fill = c("white"),
              size=6/2.8) +
   labs(title = NULL)+
   theme(legend.position = "right", legend.key.height =unit(0.2,"cm"), legend.key.width =unit(0.4,"cm"))
 
 p
 
 ggsave(paste0("./Figures/rep1_SSC_singlets.pdf"), p ,height = unit(2, "cm"), width = unit(2.5, "cm"),
        dpi = 300, useDingbats=FALSE)
 

 # Session Info
 sink("SessionInfo.txt", append=FALSE)
 cat(date())
 cat("\n")
 cat("\n")
 cat(getwd())
 cat("\n")
 cat("\n")
 print(sessionInfo())
 sink()
 