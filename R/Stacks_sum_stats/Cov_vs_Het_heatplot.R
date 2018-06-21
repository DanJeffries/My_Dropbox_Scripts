library(calibrate)


setwd("/home/djeffrie/Data/Ribe_LM/")

par(mfrow = c(1,1), mar = c(4,4,4,4))

Lab.palette <- colorRampPalette(rev(heat.colors(100)), space = "Lab")

het_cov <- read.delim("Unfiltered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "Coverage", main = "Coverage Vs Het", ylim = c(0,200))
abline(v = 0.5, lwd = 3, lty = 3)
abline(lm(het_cov$Av_cov ~ het_cov$het ))

points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)


