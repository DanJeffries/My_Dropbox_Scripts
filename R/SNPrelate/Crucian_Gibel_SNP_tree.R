library(gdsfmt)
library(SNPRelate)


## Read all Spp. vcf

setwd("//media/djeffrie/OS/Users/djeffrie/Data/Crucian_reRUN_2019/Populations_All_Spp_forSNPrelate/")

snpgdsVCF2GDS("batch_1.vcf", "batch_1.gds", method =  "biallelic.only")

snpgdsSummary("batch_1.gds")

genofile = snpgdsOpen("batch_1.gds", allow.duplicate = T)

set.seed(1000)


dissMatrix  =  snpgdsIBS(genofile,
                         sample.id = NULL,
                         snp.id = NULL,
                         num.thread=6, 
                         verbose=TRUE,
                         remove.monosnp = FALSE)

snpgdsClose(genofile)


snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)

snpgdsDrawTree(cutTree, main = "All Spp.",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular")
