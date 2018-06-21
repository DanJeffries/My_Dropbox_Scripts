#### BAYESCAN #####

source("~/RAD_programs/BayeScan2.1/R functions//plot_R.r")

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYESCAN/")

par(mfrow = c(1,1))



data <- read.delim("pr_odds_1000/Batch_1_bayescan_fst_formatted.txt", sep = "\t")
names(data) <- c("prob", "log10_PO", "qval", "alpha", "fst")

hist(data$fst, breaks = (max(data$fst) - min(data$fst))*1000)


plot_bayescan("../../ref_aligned/Ho_filtered_baysecan_fst.txt", FDR =0.001)

plot_bayescan("pr_odds_1000/Batch_1_bayescan_fst.txt", FDR =0.001)

plot_bayescan("ref_aligned//NEU_only/P17_bayescan_format_fst.txt", FDR =0.001)


plot_bayescan("BF-BOR/prior_odds_100/BF-BOR.baye_fst.txt", FDR =0.98)


plot_bayescan("BF-BOR/prior_odds_100/BF-BOR.baye_fst.txt", FDR =0.1)



plot_bayescan("bayescan_before_loc_filter/BAYESCAN_OUTPUT_fst.txt", FDR = 0.1)

head(outliers)

sorted_fgsts <- read.delim("fsts_sorted.txt", header = F)

summary(sorted_fgsts)

plot(sorted_fgsts$V1, sorted_fgsts$V2)

datafile <- read.delim("BAYESCAN_OUTPUT_fst.txt", header = T, sep = ' ')

head(datafile)

hist(datafile$log10.PO., main = "Histogram of Fst", breaks = 100)

hist(sumstats$Obs.Het)
hist(sumstats$Fis)


fst_BF_BOR <- read.delim("fst_only.txt")
length(datafile$fst)


length(fst_BF_BOR$Fst)

summary(fst_BF_BOR)
hist(fst_BF_BOR$Fst, breaks = 114)



mean(datafile[,7])

outliers_sumstats <- read.delim("outlier_only_sumstats_NOHEADERS.tsv", header = T)
non_outliers_sumstats <- read.delim("non_outlier_only_sumstats_NOHEADERS.tsv", header = T)

sumstats <- read.delim("batch_1.sumstats_no_headers.tsv", header = T)

head(sumstats)

hist(sumstats$P, breaks = 100)  ## most common allele freqs
hist(sumstats$Obs.Het, breaks = 100) ## heterozygosities
hist(sumstats$Exp.Het, breaks = 100) 
hist(sumstats$Fis, breaks = 100) ## Fis values. 
hist(sumstats$Pi, breaks = 100) ## Nuc diversity

hist(outliers$P, breaks = 100)  ## most common allele freqs
hist(outliers$Obs.Het, breaks = 100) ## heterozygosities
hist(outliers$Exp.Het, breaks = 100) 
hist(outliers$Fis, breaks = 100) ## Fis values. 
hist(outliers$Pi, breaks = 100) ## Nuc diversity




sumstats$P[1:10:2]


head(sumstats$V10[4:nrow(sumstats)])

sumstats2 <- (tail(sumstats, nrow(sumstats)-3))

levels(sumstats2$V10)
