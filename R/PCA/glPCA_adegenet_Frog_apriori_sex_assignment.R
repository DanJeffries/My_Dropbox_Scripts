#install.packages("adegenet") ## install package, only need to do this once, then you can comment out this line 

library(adegenet) ## load package into this session
library(ggplot2)

setwd("~/Data/RADseq/Rarvalis/Final_populations_outs/") ## wherever you put the input file.

## Reading in the data -------------------------------------------------------------------------------

## Rarvalis! Whole SNP dataset (something like 86,000 SNPs)

alldata <- read.PLINK('batch_1.plink.raw', chunkSize=1000, parallel = TRUE, n.cores=1, saveNbAlleles=T) ## If you are on linus you can change n.cores. But >1 core is not supported on windows.

## Sex information:

sexes <- read.delim("/home/djeffrie/Data/RADseq/Rarvalis/Final_populations_outs/sex_info.txt", header = F)

## A few checks on the data
summary(alldata)
pop(alldata)
indNames(alldata)
nInd(alldata)
ploidy(alldata)

## PCA analyses ---------------------------------------------------------------------------------------

??glPca

pca1 <- glPca(alldata, parallel = FALSE) ## again if on linux, you can set "parallel" to TRUE to use multi cores

## Use the function defined below to plot the data showing males and females separately

Sex_sorter(pca1,   ## PCA object
           sexes,  ## sex_info file
           c(1,3), ## which components to plot
           Title = "R. arvalis PCA", ## plot title
           xtitle = "PC1", ## X axis title
           ytitle = "PC3", ## Y axis title 
           Palette = c("violet", "black")) ## Colours to use



###  Now using just 150 sex linked SNPs --------------------------------------------------------------

setwd("~/Data/RADseq/Rarvalis/Final_populations_outs/") ## wherever you put the input file.

## Now just the sex linked markers -------------------------------------------------------------------

## Reading in the data -------------------------------------------------------------------------------

alldata <- read.PLINK('XY_linked_snps_plink.raw', chunkSize=1000, parallel = TRUE, n.cores=1, saveNbAlleles=T) ## If you are on linus you can change n.cores. But >1 core is not supported on windows.

## Sex information:

sexes <- read.delim("/home/djeffrie/Data/RADseq/Rarvalis/Final_populations_outs/sex_info.txt", header = F)

## A few checks on the data
summary(alldata)
pop(alldata)
indNames(alldata)
nInd(alldata)
ploidy(alldata)

## PCA analyses ---------------------------------------------------------------------------------------

pca1 <- glPca(alldata, parallel = 6) ## again if on linux, you can set "parallel" to TRUE to use multi cores

## Use the function defined below to plot the data showing males and females separately

Sex_sorter(pca1,   ## PCA object
           sexes,  ## sex_info file
           c(2,3), ## which components to plot
           Title = "R. arvalis Sex linked only PCA", ## plot title
           xtitle = "PC2", ## X axis title
           ytitle = "PC3", ## Y axis title 
           Palette = c("violet", "black")) ## Colours to use




## Sex_sorter function --------------------------------------------------------------------------------------------

## This function will take a PCA and plot it showing males and females (or whatever
## pops are specified in the sex_info file), in different colours for comparison. 

Sex_sorter <- function(pca, sex_info, PCs_to_plot = c(1,2), Title = "PCA",  xtitle = "X", ytitle = "Y", Palette = c('red','royalblue','green','violet','black'))
{
  ## Make sure sexes and PC data are in the right order
  
  pca_ordered = pca$scores[order(row.names(pca$scores)),]
  pca_df_ordered <- as.data.frame(pca_ordered) ## ggplots wants a dataframe, not a matrix
  
  sexes_ordered = sex_info[order(sexes$V1),]
  
  pca_ordered
  sexes_ordered
  
  ## Data for plotting
  
  X = pca_df_ordered[,PCs_to_plot[1]]
  Y = pca_df_ordered[,PCs_to_plot[2]]
  
  ## Labels
  Sample_names = sexes_ordered$V1
  Sex = sexes_ordered$V2
  
  ## Plotting
  
    a = ggplot(pca_df_ordered, aes(X, Y, colour=Sex)) + 
    geom_point(size = 10 ,alpha = 0.4) +
    geom_hline(aes(yintercept=0), color = "grey50") +
    geom_vline(aes(xintercept=0), color = "grey50") +
    scale_colour_manual(values=Palette) +
    annotate("text", x = X, y = Y, label = Sample_names, size = 2)+
    ggtitle(Title) + 
    xlab(xtitle) +
    ylab(ytitle) +
    theme(plot.title = element_text(size = 20, face="bold"), 
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey90"))
    
    ggsave(paste(Title, ".pdf", sep = ""), device = "pdf")
    print(a)
}


