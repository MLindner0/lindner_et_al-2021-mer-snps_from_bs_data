
# R script adjusted from: https://evodify.com/gatk-in-non-model-organism/
library(ggplot2)
library(cowplot)
library(stringr)


## ----- bayes dynamicP -----
FileLocation <- "SNP_calls"
FileNames <- paste(FileLocation, list.files(path = FileLocation, pattern = "*.bayes.dynamicP_SNPs.table"), sep="/")

FileList <- lapply(FileNames, function(x) read.csv(x, header=T, na.strings=c("","NA"), sep="\t"))

lapply(FileList, function(x) dim(x))
lapply(FileList, function(x) head(x))

sample_temp <- str_split_fixed(FileNames, "/", 2)[,2]
sample_help <- str_split_fixed(sample_temp, "\\.", 3)[,1]

Data_plots <- NULL
for(i in 1:4) {
  data <- FileList[[i]]
  data$DP <- NULL
  sample <- sample_help[i]

  # remove sample name from headers
  names(data)[grep("NA00001", names(data))] <- str_split_fixed(names(data)[grep("NA00001", names(data))], "\\.", 2)[,2]
  data$SAMPLE <- rep(sample, nrow(data))
  Data_plots <- rbind(Data_plots, data)
}

# Make plots: here DP and GQ
colfunc <- colorRampPalette(c("darkcyan", "darkgoldenrod1"))
col <- colfunc(4)

## DP - combined depth per SNP across samples (4 samples here)
(Quan <- quantile(Data_plots$DP, probs=c(0.05,0.99)))
DP_up <- as.numeric(Quan[2])
VCF_f <- Data_plots[Data_plots$DP<DP_up,] # remove snps with DP higher than 99th percentile to aid visualization

DP <- ggplot(VCF_f, aes(x=DP, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "bottom")
legend <- get_legend(DP)

DP <- ggplot(VCF_f, aes(x=DP, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  geom_vline(xintercept=10, col="black") +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

quantile(Data_plots$GQ, probs=c(0.05,0.99))
GQ <- ggplot(Data_plots, aes(x=GQ, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

(Quan <- quantile(Data_plots$QUAL, probs=c(0.05,0.99)))
QUAL <- ggplot(Data_plots, aes(x=QUAL, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

Plot <- plot_grid(DP, GQ, QUAL, labels = "AUTO", label_size = 10, scale = 0.9, ncol=3)
Plot_legend <- plot_grid(Plot, legend, rel_heights=c(10, 1), nrow=2)
save_plot("SNP_calls/Diagnostics_bayes.dynamicP_DP_GQ.pdf", Plot_legend, base_height = 2, base_width = 6)


## ----- binom -----
FileLocation <- "SNP_calls"
FileNames <- paste(FileLocation, list.files(path = FileLocation, pattern = "*.binom_SNPs.table"), sep="/")

FileList <- lapply(FileNames, function(x) read.csv(x, header=T, na.strings=c("","NA"), sep="\t"))

lapply(FileList, function(x) dim(x))
lapply(FileList, function(x) head(x))

sample_temp <- str_split_fixed(FileNames, "/", 2)[,2]
sample_help <- str_split_fixed(sample_temp, "\\.", 3)[,1]

Data_plots <- NULL
for(i in 1:4) {
  data <- FileList[[i]]
  data$DP <- NULL
  sample <- sample_help[i]

  # remove sample name from headers
  names(data)[grep("NA00001", names(data))] <- str_split_fixed(names(data)[grep("NA00001", names(data))], "\\.", 2)[,2]
  data$SAMPLE <- rep(sample, nrow(data))
  Data_plots <- rbind(Data_plots, data)
}

# Make plots: here DP and GQ
colfunc <- colorRampPalette(c("darkcyan", "darkgoldenrod1"))
col <- colfunc(4)

## DP - combined depth per SNP across samples (4 samples here)
(Quan <- quantile(Data_plots$DP, probs=c(0.05,0.99)))
DP_up <- as.numeric(Quan[2])
VCF_f <- Data_plots[Data_plots$DP<DP_up,] # remove snps with DP higher than 99th percentile to aid visualization

DP <- ggplot(VCF_f, aes(x=DP, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "bottom")
legend <- get_legend(DP)

DP <- ggplot(VCF_f, aes(x=DP, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  geom_vline(xintercept=10, col="black") +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

quantile(Data_plots$GQ, probs=c(0.05,0.99))
GQ <- ggplot(Data_plots, aes(x=GQ, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

(Quan <- quantile(Data_plots$QUAL, probs=c(0.05,0.99)))
QUAL <- ggplot(Data_plots, aes(x=QUAL, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

Plot <- plot_grid(DP, GQ, QUAL, labels = "AUTO", label_size = 10, scale = 0.9, ncol=3)
Plot_legend <- plot_grid(Plot, legend, rel_heights=c(10, 1), nrow=2)
save_plot("SNP_calls/Diagnostics_binom_DP_GQ.pdf", Plot_legend, base_height = 2, base_width = 6)

