?heatmap
alldat <- read.delim("/home/djeffrie/Data/RADseq/STOECK/Bviridis/All_samples_populations/ordered_plinks/plink_fam_sorted_for_heatmap2.raw", row.names = 1)

alldat_sorted <- alldat[order(row.names(alldat), decreasing = T), ]

alldatmat <- as.matrix(alldat_sorted )

heatmap(alldatmat, scale='none', col = c("royalblue", "orange", "darkblue"), cexRow = 1)

library(adegenet)

library(adegenet)
library(ggplot2)
source("~/Dropbox/My_Dropbox_Scripts/R/My_functions/Sex_sorter_PCA_plotting_function.txt")

## Reading in the data -------------------------------------------------------------------------------

setwd("/home/djeffrie/Data/RADseq/Pperezi/Sex_linked_markers/") 

## All SNPs (about 80,000)
XYdata <- read.PLINK('XY_linked_snps_freq_het_altered_adegenet.raw', chunkSize=1000, parallel = TRUE, n.cores=1, saveNbAlleles=T) 

## Sex information:
sexes <- read.delim("Sex_info_ID.txt", header = F)

## Running PCA ---------------------------------------------------------------------------------------

pca2 <- glPca(XYdata, parallel = TRUE, n.cores = 2, nf = 5) 

Sex_sorter(pca2,   ## PCA object
           sexes,  ## sex_info file
           c(1,2), ## which components to plot
           Title = "P perezi Sex linked only PCA", ## plot title
           xtitle = "PC1", ## X axis title
           ytitle = "PC2", ## Y axis title 
           Palette = c("violet", "black"),
           func_cex = 5) ## Size of points and text in plots) ## Colours to use

%%R -w 30 -h 30 -u cm

XY_heatmap_data <- read.delim("/home/djeffrie/Data/RADseq/Pperezi/Sex_linked_markers/XY_linked_snps_freq_het_altered_adegenet.raw", row.names = 1)
XY_heatmap_data_sorted <- XY_heatmap_data[order(row.names(XY_heatmap_data), decreasing = T), ]
XY_heatmap_data_datmat <- as.matrix(XY_heatmap_data_sorted)

heatmap(XY_heatmap_data_datmat, scale='none', col = c("royalblue", "orange", "darkblue"), cexRow = 1, main = "XY only family")
