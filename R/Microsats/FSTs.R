
library(hierfstat)


### Trying in Hierfstat -------------------------------------------------------

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DIY ABC/UK analyses - Final poolings/")


micros <- read.fstat("FSTAT_dataset_numbered.dat") ## DIYABC pooled dataset
micros_seppop <- read.fstat("FSTAT_seppop_noGF29_numbered.dat") ## Dataset with separated populations

fst.mat <- pairwise.WCfst(micros,diploid=TRUE) ## calculate fst matrix
fst_seppop.mat <- pairwise.WCfst(micros_seppop,diploid=TRUE) ## calculate fst matrix

bootstraps <- boot.ppfst(micros,nboot = 100) ## DIYABC pooled dataset
bootstraps_seppop <- boot.ppfst(micros_seppop, nboot = 100) ## Dataset with separated populations

write.fstat(fst.mat, fname = "FSTAT_DIYABC_pooled_noGF29_Fsts.txt")
write.fstat(fst_seppop.mat, fname = "FSTAT_seppop_noGF29_Fsts.txt")

write.fstat(bootstraps$ll, fname = "FSTAT_DIYABC_pooled_noGF29_Bootstrapped_CIs.txt")
write.fstat(bootstraps_seppop$ll, fname = "FSTAT_seppop_noGF29_Bootstrapped_CIs.txt")
bootstraps_seppop$ul

