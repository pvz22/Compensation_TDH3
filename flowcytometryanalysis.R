#####################################
##Data analysis and figure generation for flow cytometry data contained in "Active compensation 
#for changes in TDH3 expression mediated by direct regulators of TDH3 in Saccharomyces cerevisiae"
#Vande Zande et al, 2024
####################################


#Loading required libraries
library(flowCore)
library(flowClust)
library(flowViz)
library(ggplot2)
library(ggcyto)
library(flowStats)
library(reshape2)
library(dplyr)
#Setting ggplot2 theme for figures
THEMEMAIN <- function() {
  theme_bw() +
    theme(legend.text=element_text(size=20), axis.text = element_text(size = 20), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
#Inputting paths to directories containing raw FCS files and directory where figures will be generated. (Must be edited)
data.dir <- ""
fig.dir <- ""
#Reading in FCS data
FILENAMES <- list.files(data.dir,pattern=".fcs",recursive=TRUE,include.dirs=TRUE,full.names=TRUE)
frames <- list()
for (i in c(1:80)) {
  frames[i] <- read.FCS(FILENAMES[i],transformation=FALSE,alter.names=TRUE)
}
fs <- as(frames, "flowSet") #Doesn't work if I also include the unmixed samples

#Taking a look at FSC vs SSC to get the single cells
ggplot(data = frames[[1]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled() +
  geom_vline(xintercept = 5.5) +
  geom_vline(xintercept = 6.5) +
  geom_hline(yintercept = 5.5) +
  geom_hline(yintercept = 6.5)

ggplot(data = frames[[1]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled() +
  xlim(6.05,6.5) +
  ylim(5.95,6.5)

ggplot(data = frames[[80]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled() +
  xlim(6.05,6.5) +
  ylim(5.95,6.5)

ggplot(data = frames[[1]], aes(x = log10(FSC.A), y = log10(FSC.H))) +
  geom_density2d_filled() +
  xlim(6.05,6.75) +
  ylim(5.75,6.25)

rectgate <- rectangleGate(filterId = "Single cells", "FSC.A" = c(10^5.5,10^6.5), "SSC.A" = c(10^5.5,10^6.5))
fsrect <- Subset(fs, rectgate)

ggplot(data = fsrect[[5]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled()
ggplot(data = fsrect[[5]], aes(x = FSC.A, y =SSC.A)) +
  geom_density2d_filled()

ggplot(data = fsrect[[5]], aes(x = log10(FSC.A), y = log10(FSC.H))) +
  geom_density2d_filled()

#there is still a small population of small cells. I am going to try removing it with a cluster
Size.filter <- flowClust(fsrect[[1]],varNames=c("FSC.A"),K=2,min.count=10)
Size.filter <- flowClust(fsrect[[61]],varNames=c("FSC.A"),K=2,min.count=10)

plot(Size.filter, data = fsrect[[1]])

#Which ones will this not work for? - After running, I think this should work for all of them
for (i in 1:length(fsrect)) {
  ggplot(data = fsrect[[i]], aes(x = log10(FSC.A), y = log10(FSC.H))) +
    geom_density2d_filled() +
    ggtitle(paste0(i))
  ggsave(paste0(i,"AHplot.png"),plot = last_plot(), path = fig.dir)
}

#For loop to do this with all of them
Subsettedls1 <- list()
Subsettedls2 <- list()

for (i in 1:length(fsrect)) {
    Size.filter <- flowClust(fsrect[[i]],varNames=c("FSC.A"),K=2,min.count=10)
    plot(Size.filter, data = fsrect[[i]])
    Subsetted <- split(fsrect[[i]], Size.filter)
    Subsettedls1[[i]] <- Subsetted[[1]]
    Subsettedls2[[i]] <- Subsetted[[2]]
}

#Stopped here and saved R workspace, but got rid of frames first because it is so large.
rm(frames)
save.image("~/Documents/Documents/Minnesota/Data/Paper3flow/ADanalysisworkspace.RData")

#Lets see how well this worked
ggplot(data = Subsettedls2[[5]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled()

ggplot(data = Subsettedls2[[5]], aes(x = log10(FSC.A), y = log10(FSC.H))) +
  geom_density2d_filled()

fsSub <- as(Subsettedls2, "flowSet")

#Gating for singlets

sg <- gate_singlet(fsSub[[1]], area = "FSC.A", height = "FSC.H")
ggcyto(fsSub[[1]], aes(x = "FSC.A", y = "FSC.H")) +
  geom_hex() +
  geom_gate(sg)

Singlets <- Subset(fsrect, sg)

ggplot(data = Singlets[[1]], aes(x = log10(FSC.A), y = log10(SSC.A))) +
  geom_density2d_filled()

ggplot(data = Singlets[[1]], aes(x = log10(FSC.A), y = log10(FSC.H))) +
  geom_density2d_filled()

#Okay, pretty happy with that. Now to plot some fluorescence values.

#Now getting these labelled
x <- strsplit(FILENAMES[1:80],split = " ")
Strains <- data.frame("Strain" = c(unlist(x)[seq(from = 2, to = 160,by = 2)]))

#Reading in the key for the wells
Wellkey <- read.table(paste0(data.dir,"/Wellkey.txt"), header = 1) #This has the info about TDH3 status, etc.
Wellkey$Strain <- as.character(Wellkey$Strain)
#Just need to get rid of the Y and .fcs (and those other two wierdos)
Strains$Strain <- gsub("Y","", Strains$Strain)
Strains$Strain <- gsub(".fcs","", Strains$Strain)
Strains$Strain <- gsub("_1","", Strains$Strain)

Strains <- left_join(Strains, Wellkey, by = "Strain")

min(fsApply(Singlets, nrow)) #16844
Singlets <- transform(Singlets, 'logB2.A' = log10(`B2.A`), 'logFSC.A' = log10(`FSC.A`))
Singlets <- transform(Singlets, 'NormFluor' = `logB2.A`/`logFSC.A`)
#Just want to pull out medians to try doing stats on those
for (i in 1:nrow(Strains)) {
  Strains[i,"Median"] <- median(exprs(Singlets[[i]])[,"NormFluor"], na.rm = TRUE)
}

t.test(Strains[Strains$Strain == "3824","Median"], Strains[Strains$Strain == "3324","Median"], alternative = c("greater")) #WT, pval = 8.7e-06, FC = 1.013
t.test(Strains[Strains$Strain == "3998","Median"], Strains[Strains$Strain == "3990","Median"], alternative = c("greater")) #GCR1.1mut, pval = 0.7841, FC = 0.9982367
t.test(Strains[Strains$Strain == "4000","Median"], Strains[Strains$Strain == "3992","Median"], alternative = c("greater")) #GCR1.1mut2, pval = 0.00013, FC = 1.047
t.test(Strains[Strains$Strain == "3999","Median"], Strains[Strains$Strain == "3991","Median"], alternative = c("greater")) #GCR1.2mut, pval = 0.00067, FC = 1.013
t.test(Strains[Strains$Strain == "3997","Median"], Strains[Strains$Strain == "3989","Median"], alternative = c("greater")) #RAP1mut, pval = 0.014, FC = 1.01

t.test(Strains[Strains$Strain == "3876","Median"], Strains[Strains$Strain == "3761","Median"], alternative = c("greater")) #WT, p-value = 5.17e-06, FC = 1.086
t.test(Strains[Strains$Strain == "4087","Median"], Strains[Strains$Strain == "4085","Median"], alternative = c("greater")) #GCRabol,p-value = 0.1802, FC = 1.0054
t.test(Strains[Strains$Strain == "4086","Median"], Strains[Strains$Strain == "4084","Median"], alternative = c("greater")) #RAP1mut, p-value 0.027, FC = 1.0060

t.test(Strains[Strains$Strain == "3857","Median"], Strains[Strains$Strain == "3733","Median"], alternative = c("greater")) #p-value = 7.13e-06, FC = 1.087

i <- 1
x <- as.data.frame(exprs(Singlets[[i]]))
Events <- c(sample(x$NormFluor, size = 15000, replace = FALSE))
Strain <- rep(Strains[i,"Strain"], 15000)
PlotDF <- data.frame("Events" = Events, "Strain" = Strain)

for (i in 2:length(Singlets)) {
  x <- as.data.frame(exprs(Singlets[[i]]))
  Events <- c(sample(x$NormFluor, size = 15000, replace = FALSE))
  Strain <- rep(Strains[i,"Strain"], 15000)
  DF <- data.frame(Events = Events, Strain = Strain)
  PlotDF <- rbind(PlotDF, DF)
}

PlotDF$Replicate <- c(rep("One",20*15000),rep("Two",20*15000),rep("Three",20*15000),rep("Four",20*15000))

Boxplotmelt <- left_join(PlotDF, Wellkey, by = "Strain")

Boxplotmelt$Strain <- factor(Boxplotmelt$Strain)
Boxplotmelt$Promgeno <- factor(Boxplotmelt$Promgeno, levels = c("WT","GCR1.1mut","GCR1.1mut2","GCR1.2mut","RAP1mut","GCRmut","GCRabol"))
Strains$Promgeno <- factor(Strains$Promgeno, levels = c("WT","GCR1.1mut","GCR1.1mut2","GCR1.2mut","RAP1mut","GCRmut","GCRabol"))

ggplot(data = Boxplotmelt[Boxplotmelt$Reporter == "ptdh3:YFP" & Boxplotmelt$Strain != "3872",], aes(x = TDH3geno, y = Events)) +
  geom_violin(fill = "grey") +
  #stat_summary(fun.y = "median", geom = "point") +
  geom_point(data = Strains[Strains$Reporter == "ptdh3:YFP" & Strains$Strain != "3872",], aes(x = TDH3geno, y = Median)) +
  facet_grid(~ Promgeno) +
  THEMEMAIN() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fluorescence") +
  xlab("Background")
ggsave("TDH3reporter.pdf", plot = last_plot(), path = fig.dir, width = 7.5, height = 7)


#Performing a T-test for all events grouped together to see how that compares.
t.test(Boxplotmelt[Boxplotmelt$Strain == "3824","Events"], Boxplotmelt[Boxplotmelt$Strain == "3324","Events"], alternative = c("greater")) #WT, pval = < 2.2e-16, FC = 1.013019
t.test(Boxplotmelt[Boxplotmelt$Strain == "3998","Events"], Boxplotmelt[Boxplotmelt$Strain == "3990","Events"], alternative = c("greater")) #GCR1.1mut, pval = 1, FC = 0.9982367
t.test(Boxplotmelt[Boxplotmelt$Strain == "4000","Events"], Boxplotmelt[Boxplotmelt$Strain == "3992","Events"], alternative = c("greater")) #GCR1.1mut2, pval = < 2.2e-16, FC = 1.047735
t.test(Boxplotmelt[Boxplotmelt$Strain == "3999","Events"], Boxplotmelt[Boxplotmelt$Strain == "3991","Events"], alternative = c("greater")) #GCR1.2mut, pval = < 2.2e-16, FC = 1.022197
t.test(Boxplotmelt[Boxplotmelt$Strain == "3997","Events"], Boxplotmelt[Boxplotmelt$Strain == "3989","Events"], alternative = c("greater")) #RAP1mut, pval = < 2.2e-16, FC = 1.011505

ggplot(data = Boxplotmelt[Boxplotmelt$Reporter == "ptdh2:YFP",], aes(x = TDH3geno, y = Events)) +
  geom_violin(fill = "#003f5a") +
  #stat_summary(fun.y = "median", geom = "point") +
  geom_point(data = Strains[Strains$Reporter == "ptdh2:YFP",], aes(x = TDH3geno, y = Median), color = "white") +
  facet_grid(~ Promgeno) +
  THEMEMAIN() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fluorescence \n (arbitary units)") +
  xlab("Background")
ggsave("TDH2reporter.pdf", plot = last_plot(), path = fig.dir, width = 8, height = 7)

#T-tests for these guys
t.test(Boxplotmelt[Boxplotmelt$Strain == "3876","Events"], Boxplotmelt[Boxplotmelt$Strain == "3761","Events"], alternative = c("greater")) #WT, p-value < 2.2x10-16, FC = 1.086913
t.test(Boxplotmelt[Boxplotmelt$Strain == "4087","Events"], Boxplotmelt[Boxplotmelt$Strain == "4085","Events"], alternative = c("greater")) #GCRabol,p-value < 2.2e-16, FC = 1.005897
t.test(Boxplotmelt[Boxplotmelt$Strain == "4086","Events"], Boxplotmelt[Boxplotmelt$Strain == "4084","Events"], alternative = c("greater")) #RAP1mut, p-value < 2.2e-16, FC = 1.00627

ggplot(data = Boxplotmelt[Boxplotmelt$Reporter == "ptdh2::CFP",], aes(x = TDH3geno, y = Events)) +
  geom_violin(fill = "#003f5a") +
  #stat_summary(fun.y = "median", geom = "point") +
  geom_point(data = Strains[Strains$Reporter == "ptdh2::CFP",], aes(x = TDH3geno, y = Median), color = "white") +
  THEMEMAIN() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fluorescence \n (arbitary units)") +
  xlab("Background")
ggsave("TDH2CFPreporter.pdf", plot = last_plot(), path = fig.dir, width = 5, height = 7)

#T-tests for these guys
t.test(Boxplotmelt[Boxplotmelt$Strain == "3857","Events"], Boxplotmelt[Boxplotmelt$Strain == "3733","Events"], alternative = c("greater")) #p-value < 2.2e-16, FC = 1.089647

