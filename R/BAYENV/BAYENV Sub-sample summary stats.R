#### Sub sample summary stats ###


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV_CLOUD_COMP/Sub-samples")

SNPS <- read.delim("SNP_numbers.txt", header = F)
temps <- read.delim("Sub_samp_temp_ranges.txt", header = F)
Fsts <- read.delim("Sub_sample_fsts.txt", header = F)
names(SNPS) <- c("N_bottled", "Rep", "N_SNPs")
head(SNPS)

par(mar=c(5,5,2,2))
barplot(SNPS$N_SNPs, names.arg = paste(SNPS$N_bottled,"/", SNPS$Rep), las = 2)
title(xlab = "N Bottled pops/Rep", ylab = "Number of SNPs", main = "Number of polymorphic positions in each sub-sample")

par(mar=c(10,1,1,1))
barplot(temps$V2, names.arg = temps$V1, las = 2)
title(main = "Standardised temperature range in Sub samples")
my_line <-line(seq(1,55,1), temps$V2)
abline(coef(my_line))


par(mar=c(13,5,2,1))
barplot(Fsts$V2, names.arg = Fsts$V1, las = 2, ylab = "Fst", cex.names = 0.7)
title(main = "Average PW Fst in Sub samples")

?barplot

mean_bottle_stats <- read.delim("Mean_bottleneck_stats_per_subsample.tsv", header = T)

head(mean_bottle_stats)

mean_bottle_stats$Micro_AR <- (mean_bottle_stats$Micro_AR)/10

### Are they correlated? ## 


mean_bottle_stats.melted <- melt(mean_bottle_stats)

head(mean_bottle_stats.melted)

ggplot(mean_bottle_stats.melted, aes(X.Pop, value, fill = variable)) + geom_bar(position="dodge", stat = "identity")

## Be good to add lines to this for each stat, but will take a bit of learning of ggplot2


