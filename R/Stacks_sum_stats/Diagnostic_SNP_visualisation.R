
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_Pure_spp_only/Diagnostic_SNPs/")

spp.cols = c(heat.colors(8))

par(mfrow = c(1,3), mar=c(2,2,2,2))


## CRU_AU -----------------------------------------------------------------------

CRU_AU_all <- read.PLINK("./CRU_AU_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_AU_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(CRU_AU_all)))) ## Not too much missing data!


CRU_AU_diagnostic <-  read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/CRU_GOLD_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_AU_diagnostic, posi="topleft", yaxt = 'n', axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL, xlabels = CRU_AU_diagnostic$loc.names) ## Not too much missing data!


#CRU_GIB-------------------------------------------------------------------------

CRU_GI_all <- read.PLINK("./CRU_GI_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_GI_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(CRU_GI_all)))) ## Not too much missing data!

CRU_GI_diagnostic <- read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/CRU_GIB_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_GI_diagnostic, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL) ## Not too much missing data!


## CRU_COMM ---------------------------------------------------------------------

CRU_COMM_all <- read.PLINK("./CRU_COMM_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_COMM_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(CRU_COMM_all)))) ## Not too much missing data!

CRU_COMM_diagnostic <- read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/CRU_COMMON_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_COMM_diagnostic, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL) ## Not too much missing data!

## AU_COMM ---------------------------------------------------------------------
AU_COMM_all <- read.PLINK("./AU_COMM_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(AU_COMM_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(AU_COMM_all)))) ## Not too much missing data!


AU_COMM_diagnostic <- read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/GOLD_COMMON_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(AU_COMM_diagnostic, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL) ## Not too much missing data!

## GI_COMM ---------------------------------------------------------------------

GI_COMM_all <- read.PLINK("./GI_COMM_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(GI_COMM_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(GI_COMM_all)))) ## Not too much missing data!


GI_COMM_diagnostic <- read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/GIB_COMMON_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(GI_COMM_diagnostic, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL) ## Not too much missing data!

## AU_GI ---------------------------------------------------------------------

AU_GI_all <- read.PLINK("./AU_GI_populations/Ho_filtered_altered_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(AU_GI_all, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = seq(length(indNames(AU_GI_all)))) ## Not too much missing data!


AU_GI_diagnostic <- read.PLINK("~/Dropbox/PhD/Dans_PhD_Shared/Thesis//Chapter_5_Hybridisation_and_introgression//data/RAD/Diagnostic_snp_files/GOLD_GIB_fixed_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(AU_GI_diagnostic, posi="topleft", yaxt = 'n', xaxt = "n", axes = F, YaxisCEX = 0.4, XaxisCEX = 0.4, ylabels = NULL) ## Not too much missing data!


### Altered the gl plot slightly to incorporate axis names

myglPlot = function (x, col = NULL, legend = TRUE, posi = "bottomleft", 
                     bg = rgb(1, 1, 1, 0.5), XaxisCEX = 1, YaxisCEX = 1,  ylabels = NULL, xlabels = NULL, ...) 
{
  X <- t(as.matrix(x))
  X <- X[, ncol(X):1]
  ylabpos <- pretty(1:nInd(x), 5)
  if (is.null(col)) {
    myCol <- colorRampPalette(c("royalblue3", "firebrick1"))(max(X, 
                                                                 na.rm = TRUE) + 1)
  }
  else {
    myCol <- col
  }
  image(x = 1:nLoc(x), y = 1:nInd(x), z = X, xlab = "SNP index", ylab = "",
        yaxt = "n", col = myCol, ...)
  axis(side = 2, at = seq(nInd(x)), labels = rev(ylabels), cex.axis = YaxisCEX, las = 2)
  axis(side = 1, at = seq(nLoc(x)), labels = xlabels, cex.axis = XaxisCEX, las = 2)
  if (legend) {
    legend(posi, fill = myCol, legend = 0:max(X, na.rm = TRUE), 
           horiz = TRUE, bg = bg, title = "Number of 2nd allele")
  }
  return(invisible())
}


pop(CRU_AU_only)
nInd(CRU_AU_only)
