#Compensation ideas
library(ggplot2)
library(data.table)
library(dplyr)

rm(list = ls())
THEMEMAIN <- function() {
  theme_bw() +
    theme(legend.text=element_text(size=20), axis.text = element_text(size = 20), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

load("~/Documents/Output/Projects/Pleiotropy/DEworkspace/separatenooutlierspostcontrast.RData")
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Compensation" #Directory for figure output
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Compensation" #Directory for file output

cissamplesvec <- c("C","B","GG","M","MM")
paralogsamplesvec <- c("F","FF","C","K")
transvec <- samplesvec[!samplesvec %in% c(cissamplesvec,"BB","TT","A","F","FF","K")] #Taking out cis-reg mutants, controls, and paralogs

#Expression levels of paralogs in the TDH3 mutant
Paralogsindel <- data.frame(C[c("YGR192C","YJR009C","YJL052W"),])
Paralogsindel$Gene <- c("TDH3","TDH2","TDH1")

ggplot(data = Paralogsindel, aes(x = Gene, y = 2^log2FoldChange)) +
  geom_col() +
  THEMEMAIN() +
  geom_errorbar(aes(ymin = 2^(log2FoldChange - lfcSE), ymax = 2^(log2FoldChange + lfcSE))) +
  xlab("") +
  ylab("Fold Change in TDH3 Deletion\n(Relative to Wild Type)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1)
ggsave("TDH3delexp.pdf", plot = last_plot(), path = figdir, width = 4, height = 7)

#Expression using CPM values
# TDH1cpms <- c((countsmatrix["YJL052W","X120000"]/sum(countsmatrix[,"X120000"])*1000000))
# TDH2cpms <- c((countsmatrix["YJR009C","X120000"]/sum(countsmatrix[,"X120000"])*1000000))
# TDH3cpms <- c((countsmatrix["YGR192C","X120000"]/sum(countsmatrix[,"X120000"])*1000000))
# 
# TDH2cpms <- c((countsallauram["YJR009C",]/colSums(countsallauram)*1000000))
# TDH3cpms <- c((countsallauram["YGR192C",]/colSums(countsallauram)*1000000))
# TDH1cpms <- c((countsallauram["YJL052W",]/colSums(countsallauram)*1000000))
# 
# TDH1cpmsdf <- samples[samples$condition %in% c("C","TT"),]
# TDH1cpmsdf$SampleID <- as.character(TDH1cpmsdf$SampleID)
# for (i in 1:nrow(TDH1cpmsdf)) {
#   TDH1cpmsdf[i,"TDH1"] <- TDH1cpms[TDH1cpmsdf[i,"SampleID"]]
#   TDH1cpmsdf[i,"TDH2"] <- TDH2cpms[TDH1cpmsdf[i,"SampleID"]]
#   TDH1cpmsdf[i,"TDH3"] <- TDH3cpms[TDH1cpmsdf[i,"SampleID"]]
# }
# TDH1cpmsdf$condition <- factor(TDH1cpmsdf$condition, levels = c("TT","C"))
# ggplot(data = TDH1cpmsdf, aes(x = condition, y = TDH3)) +
#   geom_point(size = 5, color = "#de6600") +
#   THEMEMAIN() +
#   xlab("") +
#   ylab("") +
#   theme(axis.text.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm"))
# ggsave("TDH3cpms.pdf", plot = last_plot(), path = figdir, width = 2, height = 7)
# 
# ggplot(data = TDH1cpmsdf, aes(x = condition, y = TDH2)) +
#   geom_point(size = 5, color = "#003f5a") +
#   THEMEMAIN() +
#   xlab("") +
#   ylab("") +
#   theme(axis.text.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm"))
# ggsave("TDH2cpms.pdf", plot = last_plot(), path = figdir, width = 2, height = 7)
# 
# ggplot(data = TDH1cpmsdf, aes(x = condition, y = TDH1)) +
#   geom_point(size = 5, color = "#007a7a") +
#   THEMEMAIN() +
#   xlab("") +
#   ylab("") +
#   theme(axis.text.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm"))
# ggsave("TDH1cpms.pdf", plot = last_plot(), path = figdir, width = 2, height = 7)

#Fitness
# growthdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Growthrate"
# growthstats <- read.table(paste0(growthdir,"/Combstats.txt"), sep = "\t", header = 1) 
# dir2 <- "/Users/petravandezande/Documents" #Directory containing the document with the Nameskey
# Nameskey <- read.table(paste0(dir2,"/Nameskey.txt"), sep = "\t", header = 1, stringsAsFactors = FALSE)
# Nameskey$COLLECTION <- as.character(Nameskey$COLLECTION)
# #getting rid of the empty control columns
# growthstats <- growthstats[!growthstats$COLLECTION %in% c(30161,30162,30163),]
# for (i in 1:nrow(growthstats)) {
#   growthstats[i,"condition"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "CONDITION"]
#   growthstats[i,"nickname"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "NICKNAME"]
# }
# growth <- growthstats[growthstats$condition %in% paralogsamplesvec,]
# growth$Deletion <- c("TDH3del","TDH1del","TDH2/3del","TDH2del")
# growth$Deletion <- factor(growth$Deletion, levels = c("TDH1del","TDH2del","TDH3del","TDH2/3del"))
# ggplot(data = growth, aes(x = Deletion, y = REL.r.u)) +
#   geom_point(size = 5) +
#   geom_errorbar(aes(ymin = REL.r.lower, ymax = REL.r.upper)) +
#   THEMEMAIN() +
#   ylab("Growth Rate\n(Relative to Wild Type)") +
#   xlab("Paralog Deletion Strain")
# ggsave("DelGrowth.pdf", plot = last_plot(), path = figdir, width = 6, height = 8)

#Looking at dosage changes over different levels of TDH3 change
inallsets <- rownames(C)[which(rownames(C) %in% rownames(MM))]
Cislogfcs <- cbind(C[inallsets,"log2FoldChange"], B[inallsets,"log2FoldChange"], GG[inallsets,"log2FoldChange"], M[inallsets,"log2FoldChange"], MM[inallsets,"log2FoldChange"])
rownames(Cislogfcs) <- inallsets
TDH1 <- Cislogfcs["YJL052W",]
TDH2 <- Cislogfcs["YJR009C",]
YFP <- Cislogfcs["YFP",]
TDH3 <- Cislogfcs["YGR192C",]
Paralogs <- data.frame(rbind(TDH1, TDH2, TDH3, YFP))
Paralogs$Gene <- c("TDH1","TDH2","TDH3","pTDH3-YFP")
Paralogs$Gene <- factor(Paralogs$Gene, levels = c("TDH1","TDH2","TDH3","pTDH3-YFP"))
colnames(Paralogs) <- c(TDH3,"Gene")
Paralogsm <- melt(Paralogs)

Cisses <- cbind(C[inallsets,"lfcSE"], B[inallsets,"lfcSE"], GG[inallsets,"lfcSE"], M[inallsets,"lfcSE"], MM[inallsets,"lfcSE"])
rownames(Cisses) <- inallsets
TDH1 <- Cisses["YJL052W",]
TDH2 <- Cisses["YJR009C",]
YFP <- Cisses["YFP",]
TDH3 <- Cisses["YGR192C",]
Parases <- data.frame(rbind(TDH1, TDH2, TDH3, YFP))
Parases$Gene <- c("TDH1","TDH2","TDH3","pTDH3-YFP")
Parases$Gene <- factor(Parases$Gene, levels = c("TDH1","TDH2","TDH3","pTDH3-YFP"))
colnames(Parases) <- c(0,20,50,85,135,"Gene")
Parasesm <- melt(Parases)
Paraplot <- cbind(Paralogsm, sevalue = Parasesm$value)
Paraplot$variable <- as.numeric(as.character(Paraplot$variable))

ggplot(data = Paraplot, aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_point(aes(colour = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue), color = Gene)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  scale_color_manual(values = c("#007a7a","#003f5a","#de6600","#fea02f")) +
  labs(color = "Gene")
ggsave("Paralogsplot.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = Paraplot[Paraplot$Gene != "pTDH3-YFP",], aes(x = variable, y = 2^value, group = Gene)) +
  geom_point(aes(colour = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue), color = Gene)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  scale_color_manual(values = c("#007a7a","#003f5a","#de6600","#fea02f")) +
  labs(color = "Gene")
ggsave("ParalogsplotnoYFP.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Doing these separately

ggplot(data = Paraplot[Paraplot$Gene %in% c("TDH1", "TDH3"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(shape = Gene), size = 6) +
  geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue), shape = Gene, color = Gene), size = 1.5) +
  scale_shape_manual(values = c(16,1)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Fold Change\n(Relative to Wild Type)") +
  #ylim(0,2.75) +
  geom_hline(yintercept = 1) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"), legend.position = c(0.8,0.85)) +
  geom_line(aes(color = Gene)) +
  scale_color_manual(values = c("#007a7a","black")) + # ,"#de6600"
  labs(shape = "Gene", color = "Gene")
ggsave("TDH1titr.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Paraplot[Paraplot$Gene %in% c("TDH2","TDH3"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(color = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue), shape = Gene, color = Gene), size = 1.5) +
  scale_shape_manual(values = c(16,1)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Fold Change\n(Relative to Wild Type)") +
  #ylim(0,2.75) +
  geom_hline(yintercept = 1) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"), legend.position = c(0.8,0.85)) +
  geom_line(aes(color = Gene)) +
  scale_color_manual(values = c("#003f5a","black")) + #"#de6600"
  labs(shape = "Gene", color = "Gene")
ggsave("TDH2titr.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Paraplot[Paraplot$Gene %in% c("TDH3","pTDH3-YFP"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(color = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue), shape = Gene), size = 1.5) +
  scale_shape_manual(values = c(1,16)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Fold Change\n(Relative to Wild Type)") +
  #ylim(0,2.75) +
  geom_hline(yintercept = 1) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"), legend.position = c(0.75,0.15)) +
  geom_line() +
  scale_color_manual(values = c("#de6600","#fea02f")) +
  labs(shape = "Gene")
ggsave("YFPtitr.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

# ggplot(data = Paraplot[Paraplot$Gene == "TDH3",], aes(x = variable, y = 2^value, group = Gene)) +
#   geom_point(color = "#de6600", size = 5, alpha = 0.8) +
#   geom_pointrange(aes(ymin = 2^(value - sevalue), ymax = 2^(value + sevalue)), color = "#de6600") +
#   THEMEMAIN() +
#   xlab("TDH3 Expression Mutant") +
#   ylab("Fold Change\n(Relative to Wild Type)") +
#   #ylim(0,2.75) +
#   geom_hline(yintercept = 1) +
#   theme(legend.text=element_text(size=15), axis.text = element_text(size = 15), plot.margin = unit(c(0,0,0,0),"cm")) +
#   geom_line(colour = "#de6600")
# #scale_color_manual(values = c("#007a7a","#003f5a","#de6600","#fea02f")) +
# #labs(color = "Gene")
# ggsave("TDH3titr.pdf", plot = last_plot(), path = figdir, width = 5, height = 6)

#TFs
RAP1 <- Cislogfcs["YNL216W",]
GCR1 <- Cislogfcs["YPL075W",]
TDH3 <- Cislogfcs["YGR192C",]
TFs <- data.frame(rbind(RAP1, GCR1, TDH3))
TFs$Gene <- c("RAP1","GCR1","TDH3")
TFs$Gene <- factor(TFs$Gene, levels = c("RAP1","GCR1","TDH3"))
colnames(TFs) <- c(TDH3,"Gene")
TFsm <- melt(TFs)

RAP1 <- Cisses["YNL216W",]
GCR1 <- Cisses["YPL075W",]
TDH3 <- Cisses["YGR192C",]
TFses <- data.frame(rbind(RAP1, GCR1, TDH3))
TFses$Gene <- c("RAP1","GCR1","TDH3")
TFses$Gene <- factor(TFses$Gene, levels = c("RAP1","GCR1","TDH3"))
colnames(TFses) <- c(0,20,50,85,135,"Gene")
TFsesm <- melt(TFses)
Tfplot <- cbind(TFsm, sevalues = TFsesm$value)
Tfplot$variable <- as.numeric(as.character(Tfplot$variable))

ggplot(data = Tfplot, aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(colour = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Gene), size = 1.5) +
  scale_shape_manual(values = c(17,15,1)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c("#994c7d","#4a0438","black")) + #"#de6600"
  labs(color = "Gene", shape = "Gene") +
  theme(plot.margin = unit(c(0,0,0,0),"cm"), legend.position = c(0.75,0.15))
ggsave("TFsplot.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)


#Looking at how the paralogs are expressed in the RAP1 and GCR1 mutants
RAP1strains <- c("QQ","SS","W","Q")
GCR1strains <- c("PP","DD","KK","V","CC")
TFs <- c(RAP1strains, GCR1strains)
TFdf <- data.frame(Strain = TFs)
TFdf$Strain <- as.character(TFdf$Strain)
for (i in 1:nrow(TFdf)) {
  TFdf[i,"TDH3"] <- get(TFdf[i,"Strain"])["YGR192C","log2FoldChange"]
  TFdf[i,"TDH2"] <- get(TFdf[i,"Strain"])["YJR009C","log2FoldChange"]
  TFdf[i,"TDH1"] <- get(TFdf[i,"Strain"])["YJL052W","log2FoldChange"]
  TFdf[i,"pTDH3-YFP"] <- get(TFdf[i,"Strain"])["YFP","log2FoldChange"]
}

TFses <- data.frame(Strain = TFs)
TFses$Strain <- as.character(TFses$Strain)
for (i in 1:nrow(TFses)) {
  TFses[i,"TDH3se"] <- get(TFses[i,"Strain"])["YGR192C","lfcSE"]
  TFses[i,"TDH2se"] <- get(TFses[i,"Strain"])["YJR009C","lfcSE"]
  TFses[i,"TDH1se"] <- get(TFses[i,"Strain"])["YJL052W","lfcSE"]
  TFses[i,"pTDH3-YFPse"] <- get(TFses[i,"Strain"])["YFP","lfcSE"]
}
TFplot2 <- inner_join(TFdf, TFses, by = "Strain")
TFplot2$TFID <- ifelse(TFplot2$Strain %in% RAP1strains, "RAP1", "GCR1")

#Bar plot of TDH3 expression levels
#setting factor levels for easy ordering
TFplot2$Strain <- factor(TFplot2$Strain, levels = c("PP","V","DD","KK","CC","SS","QQ","W","Q"))
TFplot2 <- TFplot2[order(TFplot2$Strain),]
TFplot2$Nickname <- c("GCR1.162","GCR1.339","GCR1.281","GCR1.37","GCR1.241","RAP1.238","RAP1.54","RAP1.484","RAP1.357")
TFplot2$Nickname <- factor(TFplot2$Nickname, levels = c("GCR1.162","GCR1.339","GCR1.281","GCR1.37","GCR1.241","RAP1.238","RAP1.54","RAP1.484","RAP1.357"))

ggplot(data = TFplot2, aes(x = Nickname, y = TDH3)) +
  geom_col(aes(fill = TFID)) +
  scale_fill_manual(values = c("dimgrey","darkgrey")) +
  THEMEMAIN() +
  theme(axis.text.x = element_text(angle = 90, size = 15)) +
  ylab("TDH3 Log2 Fold Change\nRelative to Wild Type") +
  xlab("") +
  theme(legend.position = "none") +
  labs(fill = "") +
  geom_errorbar(aes(ymin = (TDH3 - TDH3se), ymax = (TDH3 + TDH3se)))
ggsave("TDH3inGRmuts.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = TFplot2, aes(x = 2^TDH3, y = 2^TDH2)) +
  #geom_point(aes(x = 2^TDH3, y = 2^TDH3, shape = TFID), color = "#de6600", size = 5, alpha = 0.7) +
  #geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH3 - TDH3se), ymax = 2^(TDH3 + TDH3se)),color = "#de6600") +
  #geom_point(aes(shape = TFID), color = "#003f5a", size = 5, alpha = 0.7) +
  geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH2 - TDH2se), ymax = 2^(TDH2 + TDH2se), shape = TFID), color = "#003f5a", size = 1.5) +
  scale_shape_manual(values = c(2,15)) +
  THEMEMAIN() +
  ylab("TDH2 Fold Change\nRelative to Wild Type") +
  xlab("TDH3 Expression\nin RAP1/GCR1 mutants") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position = c(0.2,0.8), legend.title = element_text(size = 20)) +
  labs(shape = "Mutant")
ggsave("TDH2inGRmuts.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = TFplot2, aes(x = 2^TDH3, y = 2^TDH1)) +
  #geom_point(aes(x = 2^TDH3, y = 2^TDH3, shape = TFID), color = "#de6600", size = 5, alpha = 0.7) +
  #geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH3 - TDH3se), ymax = 2^(TDH3 + TDH3se)),color = "#de6600") +
  #geom_point(aes(shape = TFID), color = "#007a7a", size = 5, alpha = 0.7) +
  geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH1 - TDH1se), ymax = 2^(TDH1 + TDH1se), shape = TFID), color = "#007a7a", size = 1.5) +
  scale_shape_manual(values = c(2,15)) +
  THEMEMAIN() +
  ylab("TDH1 Fold Change\nRelative to Wild Type") +
  xlab("TDH3 Expression\nin RAP1/GCR1 mutants") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position = c(0.5,0.8), legend.title = element_text(size = 20)) +
  labs(shape = "Mutant")
ggsave("TDH1inGRmuts.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = TFplot2, aes(x = 2^TDH3, y = 2^TFplot2[,5])) +
  #geom_point(aes(x = 2^TDH3, y = 2^TFplot2[,5], shape = TFID), color = "#fea02f", size = 5, alpha = 0.7) +
  geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TFplot2[,5] - TFplot2[,9]), ymax = 2^(TFplot2[,5] + TFplot2[,9]), shape = TFID), size = 1) +
  #geom_point(aes(x = 2^TDH3, y = 2^TDH3, shape = TFID), color = "#de6600", size = 5, alpha = 0.7) +
  #geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH3 - TDH3se), ymax = 2^(TDH3 + TDH3se)),color = "#de6600") +
  scale_shape_manual(values = c(2,15)) +
  THEMEMAIN() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position = c(0.2,0.8), legend.title = element_text(size = 20)) +
  labs(shape = "Mutant") +
  ylab("pTDH3-YFP Fold Change\nRelative to Wild Type") +
  xlab("TDH3 Expression\nin RAP1/GCR1 mutants")
ggsave("YFPinGRmuts.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

#The above plotted according to growth rate to see if TDH1 expression follows (had to re-read in growthstats)
growthdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Growthrate"
growthstats <- read.table(paste0(growthdir,"/Combstats.txt"), sep = "\t", header = 1, stringsAsFactors = FALSE)
dir2 <- "/Users/petravandezande/Documents" #Directory containing the document with the Nameskey
Nameskey <- read.table(paste0(dir2,"/Nameskey.txt"), sep = "\t", header = 1, stringsAsFactors = FALSE)
Nameskey$COLLECTION <- as.character(Nameskey$COLLECTION)
#getting rid of the empty control columns
growthstats <- growthstats[!growthstats$COLLECTION %in% c(30161,30162,30163),]
for (i in 1:nrow(growthstats)) {
  growthstats[i,"condition"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "CONDITION"]
  growthstats[i,"nickname"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "NICKNAME"]
}

TFplot2$Strain <- as.character(TFplot2$Strain)
for (i in 1:nrow(TFplot2)) {
  TFplot2[i,"REL.r.u"] <- growthstats[growthstats$condition == TFplot2[i,"Strain"],"REL.r.u"]
  TFplot2[i,"REL.r.upper"] <- growthstats[growthstats$condition == TFplot2[i,"Strain"],"REL.r.upper"]
  TFplot2[i,"REL.r.lower"] <- growthstats[growthstats$condition == TFplot2[i,"Strain"],"REL.r.lower"]
}
ggplot(data = TFplot2, aes(x = REL.r.u, y = 2^TDH1)) +
  geom_point(aes(shape = TFID), color = "#007a7a", size = 5) +
  geom_pointrange(aes(ymax = 2^(TDH1 + TDH1se), ymin = 2^(TDH1-TDH1se)), color = "#007a7a") +
  geom_errorbarh(aes(xmax = REL.r.upper, xmin = REL.r.lower), color = "#007a7a") +
  THEMEMAIN() +
  ylab("TDH1 Fold Change\n(Relative to Wild Type)") +
  xlab("Growth Rate\n(Relative to Wild Type)") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  labs(shape = "")
ggsave("TDH1TFsgrowth.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)
# 
# #RAP1 separately
# ggplot(data = TFdf[TFdf$Strain %in% RAP1strains,], aes(x = REL.r.u, y = TDH2, label = Strain)) +
#   geom_point(color = "green", size = 5) +
#   geom_point(aes(x = REL.r.u, y = TDH1), color = "pink", size = 5)+
#   geom_point(aes(x = REL.r.u, y = TFdf[TFdf$Strain %in% RAP1strains,5]), color = "purple", size = 5) +
#   geom_point(aes(x = REL.r.u, y = TDH3), color = "blue", size = 5) +
#   THEMEMAIN() +
#   ylab("Paralog Expression") +
#   xlab("Relative Growth Rate")
# ggsave("ParalogsRAP1growth.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)
# 
# ggplot(data = TFdf[TFdf$Strain %in% RAP1strains,], aes(x = TDH3, y = TDH2, label = Strain)) +
#   geom_point(color = "green", size = 5) +
#   geom_point(aes(x = TDH3, y = TDH1), color = "pink", size = 5)+
#   geom_point(aes(x = TDH3, y = TFdf[TFdf$Strain %in% RAP1strains,5]), color = "purple", size = 5) +
#   geom_point(aes(x = TDH3, y = TDH3), color = "blue", size = 5) +
#   THEMEMAIN() +
#   ylab("Paralog Expression") +
#   xlab("TDH3 Expression")
# ggsave("ParalogsRAP1.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)
# 
# #GCR1 separately
# ggplot(data = TFdf[TFdf$Strain %in% GCR1strains,], aes(x = REL.r.u, y = TDH2, label = Strain)) +
#   geom_point(color = "green", size = 5) +
#   geom_point(aes(x = REL.r.u, y = TDH1), color = "pink", size = 5)+
#   geom_point(aes(x = REL.r.u, y = TFdf[TFdf$Strain %in% GCR1strains,5]), color = "purple", size = 5) +
#   geom_point(aes(x = REL.r.u, y = TDH3), color = "blue", size = 5) +
#   THEMEMAIN() +
#   ylab("Paralog Expression") +
#   xlab("Relative Growth Rate")
# ggsave("ParalogsGCR1growth.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)
# 
# ggplot(data = TFdf[TFdf$Strain %in% GCR1strains,], aes(x = TDH3, y = TDH2, label = Strain)) +
#   geom_point(color = "green", size = 5) +
#   geom_point(aes(x = TDH3, y = TDH1), color = "pink", size = 5)+
#   geom_point(aes(x = TDH3, y = TFdf[TFdf$Strain %in% GCR1strains,5]), color = "purple", size = 5) +
#   geom_point(aes(x = TDH3, y = TDH3), color = "blue", size = 5) +
#   THEMEMAIN() +
#   ylab("Paralog Expression") +
#   xlab("TDH3 Expression")
# ggsave("ParalogsGCR1.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)

#Looking at other metabolic genes
PGK1 <- Cislogfcs["YCR012W",]
PFK2 <- Cislogfcs["YMR205C",]
ENO1 <- Cislogfcs["YGR254W",]
TDH3 <- Cislogfcs["YGR192C",]
TFs <- data.frame(rbind(PGK1, PFK2, ENO1, TDH3))
TFs$Gene <- c("PGK1","PFK2","ENO1","TDH3")
TFs$Gene <- factor(TFs$Gene, levels = c("PGK1","PFK2","ENO1","TDH3"))
colnames(TFs) <- c(TDH3,"Gene")
TFsm <- melt(TFs)
TFsm$variable <- as.numeric(as.character(TFsm$variable))

PGK1 <- Cisses["YCR012W",]
PFK2 <- Cisses["YMR205C",]
ENO1 <- Cisses["YGR254W",]
TDH3 <- Cisses["YGR192C",]
TFses <- data.frame(rbind(PGK1, PFK2, ENO1, TDH3))
TFses$Gene <- c("PGK1","PFK2","ENO1","TDH3")
TFses$Gene <- factor(TFses$Gene, levels = c("PGK1","PFK2","ENO1","TDH3"))
colnames(TFses) <- c(0,20,50,85,135,"Gene")
TFsesm <- melt(TFses)
Tfplot <- cbind(TFsm, sevalues = TFsesm$value)

ggplot(data = Tfplot, aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(colour = Gene, shape = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Gene), size = 1.5) +
  scale_shape_manual(values = c(15,16,17,1)) +
  THEMEMAIN() +
  xlab("TDH3 Expression Mutant") +
  ylab("Expression Fold Change\nRelative to Wild Type") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c("deepskyblue4","deepskyblue3","deepskyblue","black")) +
  labs(color = "Gene", shape = 'Gene')
ggsave("glycoplot.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Now looking at these expression in GCR1/RAP1 mutants
TFs <- c(RAP1strains, GCR1strains)
TFdf <- data.frame(Strain = TFs)
TFdf$Strain <- as.character(TFdf$Strain)
for (i in 1:nrow(TFdf)) {
  TFdf[i,"TDH3"] <- get(TFdf[i,"Strain"])["YGR192C","log2FoldChange"]
  TFdf[i,"PGK1"] <- get(TFdf[i,"Strain"])["YCR012W","log2FoldChange"]
  TFdf[i,"ENO1"] <- get(TFdf[i,"Strain"])["YGR254W","log2FoldChange"]
  TFdf[i,"PFK2"] <- get(TFdf[i,"Strain"])["YMR205C","log2FoldChange"]
}

TFses <- data.frame(Strain = TFs)
TFses$Strain <- as.character(TFses$Strain)
for (i in 1:nrow(TFses)) {
  TFses[i,"TDH3se"] <- get(TFses[i,"Strain"])["YGR192C","lfcSE"]
  TFses[i,"PGK1se"] <- get(TFses[i,"Strain"])["YCR012W","lfcSE"]
  TFses[i,"ENO1se"] <- get(TFses[i,"Strain"])["YGR254W","lfcSE"]
  TFses[i,"PFK2se"] <- get(TFses[i,"Strain"])["YMR205C","lfcSE"]
}
TFplot2 <- inner_join(TFdf, TFses, by = "Strain")
TFplot2$TFID <- c(rep("RAP1",4), rep("GCR1",5))

ggplot(data = TFplot2, aes(x = 2^TDH3, y = 2^PGK1)) +
  #geom_point(aes(x = 2^TDH3, y = 2^TFplot2[,5]), color = "deepskyblue3", size = 5, alpha = 0.7) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^PFK2, ymin = 2^(PFK2 - PFK2se), ymax = 2^(PFK2 + PFK2se), shape = TFID),color = "deepskyblue3", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^ENO1, ymin = 2^(ENO1 - ENO1se), ymax = 2^(ENO1 + ENO1se), shape = TFID),color = "deepskyblue", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^PGK1, ymin = 2^(PGK1 - PGK1se), ymax = 2^(PGK1 + PGK1se), shape = TFID),color = "deepskyblue4", size = 1.5) +
  #geom_point(aes(x = 2^TDH3, y = 2^TDH3), color = "#de6600", size = 5, alpha = 0.7) +
  #geom_pointrange(aes(x = 2^TDH3, ymin = 2^(TDH3 - TDH3se), ymax = 2^(TDH3 + TDH3se)),color = "#de6600") +
  scale_shape_manual(values = c(2,15)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 1) +
  geom_hline(yintercept = 1) +
  THEMEMAIN() +
  theme(legend.position = c(0.2,0.8), legend.title = element_text(size = 20)) +
  labs(shape = "Mutant") +
  ylab("Expression Fold Change\nRelative to Wild Type") +
  xlab("TDH3 Expression\nin RAP1/GCR1 mutants")
ggsave("glycoinTFs.pdf", plot = last_plot(), path = figdir, width = 6, height = 8)

#Glycolytic genes that do not show this pattern
PFK1 <- Cislogfcs["YGR240C",]
ENO2 <- Cislogfcs["YHR174W",]
TPI1 <- Cislogfcs["YDR050C",]
FBA1 <- Cislogfcs["YKL060C",]
GPM1 <- Cislogfcs["YKL152C",]
TDH3 <- Cislogfcs["YGR192C",]
TFs <- data.frame(rbind(PFK1, ENO2, TPI1, FBA1, GPM1, TDH3))
TFs$Gene <- c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3")
TFs$Gene <- factor(TFs$Gene, levels = c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3"))
colnames(TFs) <- c(TDH3,"Gene")
TFsm <- melt(TFs)
TFsm$variable <- as.numeric(as.character(TFsm$variable))

PFK1 <- Cisses["YGR240C",]
ENO2 <- Cisses["YHR174W",]
TPI1 <- Cisses["YDR050C",]
FBA1 <- Cisses["YKL060C",]
GPM1 <- Cisses["YKL152C",]
TDH3 <- Cisses["YGR192C",]
TFses <- data.frame(rbind(PFK1, ENO2, TPI1, FBA1, GPM1, TDH3))
TFses$Gene <- c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3")
TFses$Gene <- factor(TFses$Gene, levels = c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3"))
colnames(TFses) <- c(0,20,50,85,135,"Gene")
TFsesm <- melt(TFses)
Tfplot <- cbind(TFsm, sevalues = TFsesm$value)

ggplot(data = Tfplot[Tfplot$Gene !="TDH3",], aes(x = 2^variable, y = 2^value, group = Gene)) +
  #geom_point(aes(colour = Gene), size = 5, alpha = 0.8) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  geom_line(aes(colour = Gene)) +
  geom_pointrange(data = Tfplot[Tfplot$Gene == "TDH3",], aes(x = 2^variable, y = 2^value, ymin = 2^(value - sevalues), ymax = 2^(value + sevalues)), color = "black", size = 1.5, shape = 1) +
  geom_hline(yintercept = 1) +
  geom_abline(intercept = 0, slope = 1) +
  THEMEMAIN() +
  xlab("TDH3 Expression Mutant") +
  ylab("Expression Fold Change\nRelative to Wild Type") +
  scale_color_manual(values = c("orange","yellowgreen","darkgreen","blue","plum3")) +
  labs(color = "Gene")
ggsave("glycoplot2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Now looking at these expression in GCR1/RAP1 mutants
TFs <- c(RAP1strains, GCR1strains)
TFdf <- data.frame(Strain = TFs)
TFdf$Strain <- as.character(TFdf$Strain)
for (i in 1:nrow(TFdf)) {
  TFdf[i,"PFK1"] <- get(TFdf[i,"Strain"])["YGR240C","log2FoldChange"]
  TFdf[i,"ENO2"] <- get(TFdf[i,"Strain"])["YHR174W","log2FoldChange"]
  TFdf[i,"TDH3"] <- get(TFdf[i,"Strain"])["YGR192C","log2FoldChange"]
  TFdf[i,"TPI1"] <- get(TFdf[i,"Strain"])["YDR050C","log2FoldChange"]
  TFdf[i,"FBA1"] <- get(TFdf[i,"Strain"])["YKL060C","log2FoldChange"]
  TFdf[i,"GPM1"] <- get(TFdf[i,"Strain"])["YKL152C","log2FoldChange"]
}

TFses <- data.frame(Strain = TFs)
TFses$Strain <- as.character(TFses$Strain)
for (i in 1:nrow(TFses)) {
  TFses[i,"PFK1se"] <- get(TFses[i,"Strain"])["YGR240C","lfcSE"]
  TFses[i,"ENO2se"] <- get(TFses[i,"Strain"])["YHR174W","lfcSE"]
  TFses[i,"TDH3se"] <- get(TFses[i,"Strain"])["YGR192C","lfcSE"]
  TFses[i,"TPI1se"] <- get(TFses[i,"Strain"])["YDR050C","lfcSE"]
  TFses[i,"FBA1se"] <- get(TFses[i,"Strain"])["YKL060C","lfcSE"]
  TFses[i,"GPM1se"] <- get(TFses[i,"Strain"])["YKL152C","lfcSE"]
}
TFplot2 <- inner_join(TFdf, TFses, by = "Strain")
TFplot2$TFID <- c(rep("RAP1",4), rep("GCR1",5))

ggplot(data = TFplot2, aes(x = 2^TDH3, y = 2^PFK1)) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^PFK1, ymin = 2^(PFK1 - PFK1se), ymax = 2^(PFK1 + PFK1se), shape = TFID),color = "orange", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^ENO2, ymin = 2^(ENO2 - ENO2se), ymax = 2^(ENO2 + ENO2se), shape = TFID),color = "yellowgreen", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^TPI1, ymin = 2^(TPI1 - TPI1se), ymax = 2^(TPI1 + TPI1se), shape = TFID),color = "darkgreen", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^FBA1, ymin = 2^(FBA1 - FBA1se), ymax = 2^(FBA1 + FBA1se), shape = TFID),color = "blue", size = 1.5) +
  geom_pointrange(aes(x = 2^TDH3, y = 2^GPM1, ymin = 2^(GPM1 - GPM1se), ymax = 2^(GPM1 + GPM1se), shape = TFID),color = "plum3", size = 1.5) +
  scale_shape_manual(values = c(2,15)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 1) +
  geom_hline(yintercept = 1) +
  THEMEMAIN() +
  theme(legend.position = c(0.2,0.8), legend.title = element_text(size = 20)) + 
  labs(shape = "Mutant") +
  ylab("Expression Fold Change\nRelative to Wild Type") +
  xlab("TDH3 Expression\nin RAP1/GCR1 mutants")
ggsave("glycoinTFs2.pdf", plot = last_plot(), path = figdir, width = 6, height = 8)
