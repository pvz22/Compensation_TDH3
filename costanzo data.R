#Looking up Costanzo (2010) data on TDH1,2, and 3
library(data.table)

data <- fread("/Downloads/sgadata_costanzo2009_rawdata_101120.txt", sep = "\t", header = FALSE)

colnames(data) <- c("QueryORF","QueryName","ArrayORF","ArrayName","GIS","SD","p-value","QuerySMF","QSMFSD","ArraySMF","ASMFSD","DMF","DMFSD")

"TDH3" %in% data$QueryName
"TDH2" %in% data$ArrayName
"TDH1" %in% data$QueryName
data[data$QueryName == "TDH3" & data$ArrayName == "TDH2",]
data[data$QueryName == "TDH2" & data$ArrayName == "TDH3",]

data[data$QueryName == "TDH3" & data$ArrayName == "TDH1",] #Does not exist
data[data$QueryName == "TDH1" & data$ArrayName == "TDH3",]

data[data$QueryName == "TDH1" & data$ArrayName == "TDH2",] #They have basically the same single mutant growth.
