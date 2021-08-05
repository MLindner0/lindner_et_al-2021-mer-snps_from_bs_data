
### Get csv file for python dictionary
### for
### RBC WGBS data (ERC selection lines)
###
### BM Lindner v2020

### -------------- load required packages & set pathes -----------------------------

library(dplyr)
library(tidyr)
library(stringr)

p_in <- "/home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data"
p_out <- "/home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/src"


### -------------- load & format data ----------------------------------------------

setwd(p_in)
FileNames_seq2018 = list.files(path="seq2018_full", pattern = "*fastq$")
FileNames_reseq2020 = list.files(path="reseq2020_full", pattern = "*fastq.gz$")
FileNames_reseq2020_2 = list.files(path="reseq2020_SP_full", pattern = "*fastq.gz$")

# set up loop to get input file for python dictionary
NewData <- NULL
In <- FileNames_seq2018
for(f in seq(1, length(In), by=1)) {
  Parts <- str_split_fixed(In[f],"_",8)
  if(Parts[,8]==""){
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], sep="_")
    Barcode <- Parts[,4]
    Lane <- str_split_fixed(Parts[,5],"00",2)[2]
  }
  else {
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], Parts[,4], sep="_")
    Barcode <- Parts[,5]
    Lane <- str_split_fixed(Parts[,6],"00",2)[2]
  }
Data <- data.frame(Sample=Sample, Barcode=Barcode, Lane=Lane)
NewData <- rbind(NewData, Data)
}
NewData$Library <- rep(1, nrow(NewData))
NewData$Flowcell <- rep("H3FCMDSXX", nrow(NewData))
Run_seq2018 <- gather(NewData, names(NewData[2:5]), key = "Key", value = "Value")

NewData <- NULL
In <- FileNames_reseq2020
for(f in seq(1, length(In), by=1)) {
  Parts <- str_split_fixed(In[f],"_",8)
  if(Parts[,8]==""){
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], sep="_")
    Barcode <- Parts[,4]
    Lane <- str_split_fixed(Parts[,5],"00",2)[2]
  }
  else {
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], Parts[,4], sep="_")
    Barcode <- Parts[,5]
    Lane <- str_split_fixed(Parts[,6],"00",2)[2]
  }
Data <- data.frame(Sample=Sample, Barcode=Barcode, Lane=Lane)
NewData <- rbind(NewData, Data)
}
NewData$Library <- rep(2, nrow(NewData))
NewData$Flowcell <- rep("H2MY3DSXY", nrow(NewData))
Run_reseq2020 <- gather(NewData, names(NewData[2:5]), key = "Key", value = "Value")

NewData <- NULL
In <- FileNames_reseq2020_2
for(f in seq(1, length(In), by=1)) {
  Parts <- str_split_fixed(In[f],"_",8)
  if(Parts[,8]==""){
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], sep="_")
    Barcode <- Parts[,4]
    Lane <- str_split_fixed(Parts[,5],"00",2)[2]
  }
  else {
    Sample <- paste(Parts[,1], Parts[,2], Parts[,3], Parts[,4], sep="_")
    Barcode <- Parts[,5]
    Lane <- str_split_fixed(Parts[,6],"00",2)[2]
  }
Data <- data.frame(Sample=Sample, Barcode=Barcode, Lane=Lane)
NewData <- rbind(NewData, Data)
}
NewData$Library <- rep(3, nrow(NewData))
NewData$Flowcell <- rep("HVL2MDRXX", nrow(NewData))
Run_reseq2020_2 <- gather(NewData, names(NewData[2:5]), key = "Key", value = "Value")

# save data
setwd(p_out)
write.table(Run_seq2018, "Run_seq2018.csv", quote=F, sep=",", row.names=F)
write.table(Run_reseq2020, "Run_reseq2020.csv", quote=F, sep=",", row.names=F)
write.table(Run_reseq2020_2, "Run_reseq2020_SP.csv", quote=F, sep=",", row.names=F)
