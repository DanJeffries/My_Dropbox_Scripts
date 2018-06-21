library(calibrate)

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/populations_tests/")

cov <- read.delim("populations__p31_r07/Average_coverage.txt", header = F)
black_cov <- read.delim("populations__p31_r07/Average_Pure_cru_blacklisted_coverage.txt", header = F)

par(mfrow = c(1,2), mar = c(2,2,6,2))

hist(cov$V1, breaks = max(cov$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(cov$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(cov$V1)),cex = 1, pos = 4)

hist(black_cov$V1, breaks = max(black_cov$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(black_cov$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(black_cov$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M4_m4_N4/populations_tests/populations__p31_r07/")

cov <- read.delim("Average_coverage.txt", header = F)
black_cov <- read.delim("Average_Pure_cru_blacklisted_coverage.txt", header = F)

par(mfrow = c(1,2))

hist(cov$V1, breaks = max(cov$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage")
textxy(50,1100, paste("mean coverage =", round(mean(cov$V1), 3)), cex = 1, pos = 4)
textxy(50, 1000, paste("n Tags = ", length(cov$V1)),cex = 1, pos = 4)

hist(black_cov$V1, breaks = max(black_cov$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage")
textxy(50,1400, paste("mean coverage =", round(mean(black_cov$V1), 3)), cex = 1, pos = 4)
textxy(50, 1300, paste("n Tags = ", length(black_cov$V1)),cex = 1, pos = 4)


## Heterozygosity vs coverage for finding split loci! ----------------------------------------------------


par(mfrow = c(2,7), mar = c(4,4,4,1))

## Denovo ##### -------------------------------------------------------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M4_m4_N4/populations_tests/")

Lab.palette <- colorRampPalette(rev(heat.colors(100)), space = "Lab")

het_cov <- read.delim("populations__p25_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "Coverage", main = "p25")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p26_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p26")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p27_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p27")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p28_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p28")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p29_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p29")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p30_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p30")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p31_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, xlab = "Heterozygosity", ylab = "", main = "p31")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)



## Ref_aligned ##### -------------------------------------------------------------------------------------------

#par(mfrow = c(2,7), mar = c(4,4,4,1))

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/")

Lab.palette <- colorRampPalette(rev(heat.colors(100)), space = "Lab")

het_cov <- read.delim("populations__p25_r07/Pure_cru_Ho_filtered_cov_vs_het_ammended.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "Coverage", main = "p25")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p26_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p26")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p27_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p27")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p28_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p28")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p29_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p29")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p30_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p30")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)

het_cov <- read.delim("populations__p31_r07/Pure_cru_Ho_filtered_cov_vs_het.txt")
smoothScatter(het_cov$het, het_cov$Av_cov, colramp = Lab.palette, ylim = c(0,80), xlab = "Heterozygosity", ylab = "", main = "p31")
#points(het_cov$het, het_cov$Av_cov, pch = 16, cex = 0.1)


### Ho_filtering tag coverage ## ------------------------------------------------------------------------------

## De novo ---------------------
par(mfrow = c(4,2))

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED//All_samples_M4_m4_N4/populations_tests/")

#P25
unfiltered <- read.delim("populations__p25_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p25_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P26
unfiltered <- read.delim("populations__p26_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p26_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P27
unfiltered <- read.delim("populations__p27_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p27_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P28
unfiltered <- read.delim("populations__p28_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p28_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P29
unfiltered <- read.delim("populations__p29_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p29_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P30
unfiltered <- read.delim("populations__p30_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p30_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P31
unfiltered <- read.delim("populations__p31_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p31_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

## Ref aligned -----------------
par(mfrow = c(4,2))

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/")

#P25
unfiltered <- read.delim("populations__p25_r07/T unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p25_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P26
unfiltered <- read.delim("populations__p26_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p26_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P27
unfiltered <- read.delim("populations__p27_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p27_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P28
unfiltered <- read.delim("populations__p28_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p28_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P29
unfiltered <- read.delim("populations__p29_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p29_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P30
unfiltered <- read.delim("populations__p30_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p30_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P31
unfiltered <- read.delim("populations__p31_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p31_r07/Pure_cru_ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

## Pure Crucian only ------------------------------
par(mfrow = c(1,2))
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_tests/")

#P19
unfiltered <- read.delim("populations__p19_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p19_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P18
unfiltered <- read.delim("populations__p18_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p18_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P17
unfiltered <- read.delim("populations__p17_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p17_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 240, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P16
unfiltered <- read.delim("populations__p16_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p16_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P15
unfiltered <- read.delim("populations__p15_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p15_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)

#P14
unfiltered <- read.delim("populations__p14_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p14_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)


#P13
unfiltered <- read.delim("populations__p13_r07/unfiltered_av_tag_cov.txt", header = F)
Hofiltered <- read.delim("populations__p13_r07/Ho_filtered_av_tag_cov.txt", header = F)

hist(unfiltered$V1, breaks = max(unfiltered$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(unfiltered$V1), 2)), cex = 1, pos = 4)
textxy(45, 280, paste("n Tags = ", length(unfiltered$V1)),cex = 1, pos = 4)

hist(Hofiltered$V1, breaks = max(Hofiltered$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage", xlim = c(0,100))
textxy(45,300, paste("mean coverage =", round(mean(Hofiltered$V1), 3)), cex = 1, pos = 4)
textxy(45,280, paste("n Tags = ", length(Hofiltered$V1)),cex = 1, pos = 4)
title(main = "Reference aligned", outer = T, line = -1.5)





