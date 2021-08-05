
### Make meta data csv-file for gemBS
### for
### WGBS data (ERC selection lines)
###
### M Lindner v2020

### -------------- load required packages & set pathes -----------------------------

library(dplyr)
library(tidyr)
library(stringr)

p_in <- "/home/NIOO.INT/melaniel/projects/WGBS_Snakemake_SnpCalling/raw_data/full_sample_name/"
p_out <- "/home/NIOO.INT/melaniel/projects/WGBS_Snakemake_SnpCalling/gemBS/"
options(width=200)


### -------------- load & format data ----------------------------------------------

setwd(p_in)
FileNames_seq2018 = list.files(pattern = "*fastq$")

# set up loop to get input file for python dictionary
NewData <- NULL
In <- FileNames_seq2018
for(f in 1:length(In)) {
  Parts <- str_split_fixed(In[f],"_",8)
  #Barcode <- Parts[,5]
  Sample <- paste(Parts[,1], Parts[,2], Parts[,3], Parts[,4], sep="_")
  Sample_short <- paste(Parts[,3], Parts[,4], sep="_")
  Run1 <- data.frame(Barcode=Sample, Library="KapaBiosystems-1", file_id=paste("seq2018", Sample, sep="/"),
  end_1=paste("seq2018", Sample, "R1_val_1.fq.gz", sep="/"), end_2=paste("seq2018", Sample, "R2_val_2.fq.gz", sep="/"))
  Run2 <- data.frame(Barcode=Sample, Library="KapaBiosystems-2", file_id=paste("reseq2020", Sample, sep="/"),
  end_1=paste("reseq2020", Sample, "R1_val_1.fq.gz", sep="/"), end_2=paste("reseq2020", Sample, "R2_val_2.fq.gz", sep="/"))
  NewData <- rbind(NewData, Run1, Run2)
}

# save data
setwd(p_out)
write.table(NewData, "Metadata.csv", quote=F, sep=",", row.names=F)
