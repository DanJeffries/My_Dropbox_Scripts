
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_Pure_spp_only/Diagnostic_SNPs/")

spp.cols = c(heat.colors(8))

par(mar=c(5,7,2,2))


## CRU_AU
CRU_AU_all_snps <- read.PLINK("CRU_AU_populations/All_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(CRU_AU, posi="topleft", yaxt = 'n')

myglPlot(CRU_AU_all_snps, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = names$V1) ## Not too much missing data!

CRU_AU_all_diagnostic <- read.PLINK("CRU_AU_populations/All_diag_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
CRU_AU_sub<- CRU_AU_all_diagnostic[c(100:199)]

myglPlot(CRU_AU_all_diagnostic, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = indNames(CRU_AU_all_diagnostic)) ## Not too much missing data!

CRU_AU_super_diagnostic <-  read.PLINK("CRU_AU_populations/Super_diagnostic_tags_3_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_AU_super_diagnostic, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = indNames(CRU_AU_super_diagnostic)) ## Not too much missing data!


names <- read.delim("CRU_AU_populations/names.txt", header = F)
indNames(CRU_AU_super_diagnostic) <- names$V1


#CRU_GIB
CRU_GI_all_snps <- read.PLINK("CRU_GI_populations/plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

myglPlot(CRU_GI_all_snps, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = indNames(CRU_GI_all_snps)) ## Not too much missing data!

CRU_GI_all_diagnostic <- read.PLINK("CRU_GI_populations/All_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

myglPlot(CRU_GI_all_diagnostic, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = indNames(CRU_GI_all_diagnostic)) ## Not too much missing data!

CRU_GI_super_diagnostic3 <-  read.PLINK("CRU_GI_populations/Super_diagnostic_tags_3_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
myglPlot(CRU_GI_super_diagnostic3, posi="topleft", yaxt = 'n', axes = F, XaxisCEX = 1, YaxisCEX = 0.3, ylabels = indNames(CRU_GI_super_diagnostic3), xlabels = locNames(CRU_GI_super_diagnostic3)) ## Not too much missing data!


#CRU_DANUBE
CRU_NEU_DAN_all_snps <- read.PLINK("CRU_NEU_DANUBE_populations/All_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

myglPlot(CRU_NEU_DAN_all_snps, posi="topleft", yaxt = 'n', axes = F, axisCEX = 0.5, ylabels = indNames(CRU_NEU_DAN_all_snps)) ## Not too much missing data!

CRU_NEU_DAN_all_diagnostic <- read.PLINK("CRU_NEU_DANUBE_populations/All_diagnostic_snps_plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

myglPlot(CRU_NEU_DAN_all_diagnostic, posi="topleft", yaxt = 'n', axes = F, XaxisCEX = 1, YaxisCEX = 0.3, ylabels = indNames(CRU_NEU_DAN_all_diagnostic), xlabels = locNames(CRU_NEU_DAN_all_diagnostic)) ## Not too much missing data!



names <- read.delim("CRU_NEU_DANUBE_populations/names.txt", header = F)
indNames(CRU_NEU_DAN_all_snps) <- names$V1
indNames(CRU_NEU_DAN_all_diagnostic) <- names$V1


names <- read.delim("CRU_GI_populations/names_2.txt", header = F)
indNames(CRU_GI_all_snps) <- names$V1
indNames(CRU_GI_all_diagnostic) <- names$V1
indNames(CRU_GI_super_diagnostic3) <- names$V1

locNames(CRU_GI_super_diagnostic3)

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


