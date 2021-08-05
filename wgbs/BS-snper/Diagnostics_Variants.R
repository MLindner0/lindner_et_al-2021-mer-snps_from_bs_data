
# R script adjusted from: https://evodify.com/gatk-in-non-model-organism/
library(ggplot2)
library(cowplot)
library(stringr)


## ----- new Bis-SNP version (with Base Quality Score Recalibration) -----
FileLocation <- "SNP_calls"
FileNames <- paste(FileLocation, list.files(path = FileLocation, pattern = "*.table"), sep="/")

FileList <- lapply(FileNames, function(x) read.csv(x, header=T, na.strings=c("","NA"), sep="\t"))

lapply(FileList, function(x) dim(x))
lapply(FileList, function(x) head(x))

Data_plots <- NULL
for(i in 1:4) {
  data <- FileList[[i]]
  data$DP <- NULL
  sample <- str_split_fixed(names(data)[ncol(data)], "\\.", 9)[,6]

  # remove sample name from headers
  names(data)[grep(sample, names(data))] <- str_split_fixed(names(data)[grep(sample, names(data))], "\\.", 9)[,9]
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
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  geom_vline(xintercept=10, col="black") +
  theme_classic() + theme(legend.position = "none")

(Quan <- quantile(Data_plots$QUAL, probs=c(0.05,0.99)))
help <- Data_plots[Data_plots$QUAL<as.numeric(Quan)[2],]
QUAL <- ggplot(help, aes(x=QUAL, fill=SAMPLE)) + geom_density(alpha=.5, size=0.1) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none")

Plot <- plot_grid(DP, QUAL, labels = "AUTO", label_size = 10, scale = 0.9, ncol=2)
Plot_legend <- plot_grid(Plot, legend, rel_heights=c(10, 1), scale = c(1,0.5), nrow=2)
save_plot("SNP_calls/Diagnostics_DP.pdf", Plot_legend, base_height = 2.5, base_width = 6)

