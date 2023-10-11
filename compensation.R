#Script generating figures for compensation for TDH3 by TDH2
#Libraries needed
library(ggplot2)
library(data.table)
library(dplyr)
#clearing environment
rm(list = ls())
#ggplot theme
THEMEMAIN <- function() {
  theme_bw() +
    theme(legend.text=element_text(size=20), axis.text = element_text(size = 20), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
#Data already analyzed as published in Vande Zande, Hill and Wittkopp (2022)
load("/DEworkspace/separatenooutlierspostcontrast.RData")

figdir <- "" #Directory for figure output
outputdir <- "" #Directory for file output

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


#Looking at other glycolytic genes
PGK1 <- Cislogfcs["YCR012W",]
PFK2 <- Cislogfcs["YMR205C",]
ENO1 <- Cislogfcs["YGR254W",]
TDH3 <- Cislogfcs["YGR192C",]
PFK1 <- Cislogfcs["YGR240C",]
ENO2 <- Cislogfcs["YHR174W",]
TPI1 <- Cislogfcs["YDR050C",]
FBA1 <- Cislogfcs["YKL060C",]
GPM1 <- Cislogfcs["YKL152C",]
TFs <- data.frame(rbind(PFK1, ENO2, TPI1, FBA1, GPM1, TDH3, PGK1, PFK2, ENO1))
TFs$Gene <- c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3","PGK1","PFK2","ENO1")
TFs$Gene <- factor(TFs$Gene)
colnames(TFs) <- c(TDH3,"Gene")
TFsm <- melt(TFs)
TFsm$variable <- as.numeric(as.character(TFsm$variable))

PGK1 <- Cisses["YCR012W",]
PFK2 <- Cisses["YMR205C",]
ENO1 <- Cisses["YGR254W",]
TDH3 <- Cisses["YGR192C",]
PFK1 <- Cisses["YGR240C",]
ENO2 <- Cisses["YHR174W",]
TPI1 <- Cisses["YDR050C",]
FBA1 <- Cisses["YKL060C",]
GPM1 <- Cisses["YKL152C",]
TFses <- data.frame(rbind(PFK1, ENO2, TPI1, FBA1, GPM1, TDH3, PGK1, PFK2, ENO1))
TFses$Gene <- c("PFK1","ENO2","TPI1","FBA1","GPM1","TDH3","PGK1","PFK2","ENO1")
TFses$Gene <- factor(TFses$Gene)
colnames(TFses) <- c(0,20,50,85,135,"Gene")
TFsesm <- melt(TFses)
Tfplot <- cbind(TFsm, sevalues = TFsesm$value)
Tfplot$Type <- rep("Cis",nrow(Tfplot))


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
  TFdf[i,"TDH3"] <- get(TFdf[i,"Strain"])["YGR192C","log2FoldChange"]
  TFdf[i,"PGK1"] <- get(TFdf[i,"Strain"])["YCR012W","log2FoldChange"]
  TFdf[i,"ENO1"] <- get(TFdf[i,"Strain"])["YGR254W","log2FoldChange"]
  TFdf[i,"PFK2"] <- get(TFdf[i,"Strain"])["YMR205C","log2FoldChange"]
}

TFdft <- t(TFdf[,2:10])
TFdft <- data.frame(TFdft)
colnames(TFdft) <- TFdft["TDH3",]
TFdft$Gene <- rownames(TFdft)
TFm <- melt(TFdft)

TFses <- data.frame(Strain = TFs)
TFses$Strain <- as.character(TFses$Strain)
for (i in 1:nrow(TFses)) {
  TFses[i,"PFK1"] <- get(TFses[i,"Strain"])["YGR240C","lfcSE"]
  TFses[i,"ENO2"] <- get(TFses[i,"Strain"])["YHR174W","lfcSE"]
  TFses[i,"TDH3"] <- get(TFses[i,"Strain"])["YGR192C","lfcSE"]
  TFses[i,"TPI1"] <- get(TFses[i,"Strain"])["YDR050C","lfcSE"]
  TFses[i,"FBA1"] <- get(TFses[i,"Strain"])["YKL060C","lfcSE"]
  TFses[i,"GPM1"] <- get(TFses[i,"Strain"])["YKL152C","lfcSE"]
  TFses[i,"TDH3"] <- get(TFses[i,"Strain"])["YGR192C","lfcSE"]
  TFses[i,"PGK1"] <- get(TFses[i,"Strain"])["YCR012W","lfcSE"]
  TFses[i,"ENO1"] <- get(TFses[i,"Strain"])["YGR254W","lfcSE"]
  TFses[i,"PFK2"] <- get(TFses[i,"Strain"])["YMR205C","lfcSE"]
}

TFsedft <- t(TFses[,2:10])
TFsedft <- data.frame(TFsedft)
TFsedft$Gene <- rownames(TFsedft)
TFsemelt <- melt(TFsedft)

Tfplot3 <- cbind(TFm, sevalues = TFsemelt$value)
Tfplot3$Type <- rep("Trans",nrow(Tfplot3))

Plottogether <- rbind(Tfplot, Tfplot3)
Plottogether$variable <- as.numeric(as.character(Plottogether$variable))
#Adding info on which strains are from RAP1, GCR1, or TDH3 mutants (they are currently in order)
Plottogether$Mutant <- c(rep("TDH3", 45), rep("RAP1",36), rep("GCR1",45))


ggplot(data = Plottogether[Plottogether$Type == "Cis" & Plottogether$Gene %in% c("PFK1","PFK2"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(data = Plottogether[Plottogether$Type == "Trans" & Plottogether$Gene %in% c("PFK1","PFK2"),], aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Mutant), size = 1.5) +
  scale_shape_manual(values = c(2,15)) +
  labs(color = "Gene") +
  scale_color_manual(values = c("grey","deepskyblue3")) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("PFKs.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Plottogether[Plottogether$Type == "Cis" & Plottogether$Gene %in% c("FBA1","TPI1"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(data = Plottogether[Plottogether$Type == "Trans" & Plottogether$Gene %in% c("FBA1","TPI1"),], aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Mutant), size = 1.5) + #shape = TF
  scale_shape_manual(values = c(2,15)) +
  labs(color = "Gene") +#shape = "TF"
  scale_color_manual(values = c("darkgrey","lightgrey")) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("FBATIP.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Plottogether[Plottogether$Type == "Cis" & Plottogether$Gene %in% c("PGK1"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(data = Plottogether[Plottogether$Type == "Trans" & Plottogether$Gene %in% c("PGK1"),], aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Mutant), size = 1.5) + #shape = TF
  scale_shape_manual(values = c(2,15)) +
  labs(color = "Gene") +#shape = "TF"
  scale_color_manual(values = c("deepskyblue4")) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("PGK.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Plottogether[Plottogether$Type == "Cis" & Plottogether$Gene %in% c("GPM1"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(data = Plottogether[Plottogether$Type == "Trans" & Plottogether$Gene %in% c("GPM1"),], aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Mutant), size = 1.5) + #shape = TF
  scale_shape_manual(values = c(2,15)) +
  labs(color = "Gene") +#shape = "TF"
  scale_color_manual(values = c("darkgrey")) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("GPM.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

ggplot(data = Plottogether[Plottogether$Type == "Cis" & Plottogether$Gene %in% c("ENO1","ENO2"),], aes(x = 2^variable, y = 2^value, group = Gene)) +
  geom_pointrange(aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene), size = 1.5) +
  THEMEMAIN() +
  xlab("TDH3 Expression Level") +
  ylab("Expression Fold Change\n(Relative to Wild Type)") +
  geom_line(aes(colour = Gene)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(data = Plottogether[Plottogether$Type == "Trans" & Plottogether$Gene %in% c("ENO1","ENO2"),], aes(ymin = 2^(value - sevalues), ymax = 2^(value + sevalues), color = Gene, shape = Mutant), size = 1.5) + #shape = TF
  scale_shape_manual(values = c(2,15)) +
  labs(color = "Gene") +#shape = "TF"
  scale_color_manual(values = c("deepskyblue3","grey")) +
  theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("ENOs.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)
