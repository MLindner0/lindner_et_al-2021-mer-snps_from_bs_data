# R script adjusted from: https://evodify.com/gatk-in-non-model-organism/
library(ggplot2)
library(cowplot)
library(stringr)

VCFsnps <- read.csv('variants/GVCFall_SNPs.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(VCFsnps)
VCF <- VCFsnps

snps <- "darkcyan"


## DP - combined depth per SNP across samples (4 samples here)
quantile(VCFsnps$DP, probs = seq(0,1,by=0.01))
range(VCF$DP) # up to 25,809

Quan <- as.numeric(quantile(VCF$DP, probs = c(0.01,0.99))) # 1th and 99th percentile: 36,319

VCF_f <- VCF[VCF$DP<=Quan[2],] # remove snps with DP higher than 99th percentile to aid visualization
DP <- ggplot(VCF_f, aes(x=DP)) + geom_density(fill=snps) +
  theme_classic() + scale_y_continuous(expand=c(0,0))


## QD - variant confidence standardized by depth
VCF_f <- VCF[!is.na(VCF$QD),] # remove NA
QD <- ggplot(VCF_f, aes(x=QD)) + geom_density(fill=snps) +
  geom_vline(xintercept=2, col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0))


## FS - strand bias in support for REF vs ALT allele calls
## More bias is indicative of false positive calls (Auwera et al, 2013)
## Threshold: 60.00 (Auwera et al, 2013)
range(VCF$FS)

FS <- ggplot(VCF, aes(x=FS)) + geom_density(fill=snps) +
  geom_vline(xintercept=60, col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0))

## MQ - Mapping quality of a SNP
## Threshold: 40.00 (Auwera et al, 2013)
# Quan <- as.numeric(quantile(VCFsnps$MQ, probs = 0.01))

MQ <- ggplot(VCF, aes(x=MQ)) + geom_density(fill=snps) +
  geom_vline(xintercept=40, col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0))

## MQRankSum - Rank sum test for mapping qualities of REF vs. ALT reads
## Threshold: 12.50 (Auwera et al, 2013)

VCF_f <- VCF[!is.na(VCF$MQRankSum),] # remove NA
VCF_f <- VCF_f[VCF_f$MQRankSum>=(-20)&VCF_f$MQRankSum<=20,] # remove snps with MQRankSum higer/smaller than +/- 20 to aid visualization
MQRankSum <- ggplot(VCF_f, aes(x=MQRankSum)) + geom_density(fill=snps) +
  geom_vline(xintercept=c(-12.5,12.5), col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0))


## SOR - sequencing bias in which one DNA strand is favored over the other
## Threshold: 4.00
SOR <- ggplot(VCF, aes(x=SOR)) + geom_density(fill=snps) +
  geom_vline(xintercept=4, col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0))


## ReadPosRankSum - do all the reads support a SNP call tend to be near the end of a read
## Distance from the end of the read for reads with the alternate allele
## If the alternate allele is only seen near the ends of reads; this is indicative of error
VCF_f <- VCF[!is.na(VCF$ReadPosRankSum),] # remove NA
ReadPosRankSum <- ggplot(VCF_f, aes(x=ReadPosRankSum)) + geom_density(fill=snps) +
  geom_vline(xintercept=c(-10,10), col="black") +
  theme_classic() + scale_y_continuous(expand=c(0,0))

Plot <- plot_grid(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, labels = "AUTO", label_size = 10, scale = 0.9, ncol=2)
save_plot("variants/Diagnostics_VariantScores.pdf", Plot, base_height = 8, base_width = 6)


## per sample/GT coverage
VCFsnps_GT <- VCFsnps[,-c(3:10)]
Percentiles_5th_99th <- NULL
for (i in 1:4) {
  s <- 2+2*i
  sample <- str_split_fixed(names(VCFsnps_GT)[s], "\\.", 2)[,1]

  DP <- VCFsnps_GT[,s]
  Percentiles <- as.numeric((quantile(DP[!is.na(DP)], probs=c(0.05,0.99))))
  Data <- data.frame(Sample=sample, Percentile_5th=Percentiles[1], Percentile_99th=Percentiles[2])
  Percentiles_5th_99th <- rbind(Percentiles_5th_99th, Data)
}
DP_low <- floor(mean(Percentiles_5th_99th$Percentile_5th))
DP_up <- ceiling(mean(Percentiles_5th_99th$Percentile_99th))

# plot DP distribution
col <- "darkcyan"

Plots_DP <- NULL
for (i in 1:4) {
  s <- 2+2*i
  sample <- str_split_fixed(names(VCFsnps_GT)[s], "\\.", 2)[,1]

  DP <- VCFsnps_GT[,s]
  DP_noNA <- DP[!is.na(DP)]
  Data_f <- data.frame(DP=DP_noNA[DP_noNA<=150])

  Plot <- ggplot(Data_f, aes(x=DP)) + geom_density(fill=col, size=0.1) +
    geom_vline(xintercept=c(DP_low-1, DP_up), col="black") +
    theme_classic() + scale_y_continuous(expand=c(0,0))

  assign(sample, Plot)
  Plots_DP <- append(Plots_DP, sample, after = length(Plots_DP))
}

Plot <- plot_grid(get(Plots_DP[1]), get(Plots_DP[2]), get(Plots_DP[3]), get(Plots_DP[4]), labels = "AUTO", label_size = 10, scale = 0.9, ncol=2)
save_plot("variants/Diagnostics_DP.pdf", Plot, base_height = 5, base_width = 6)
