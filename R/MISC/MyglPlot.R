library(adegenet)

setwd("/home/djeffrie/Data/RADseq/STOECK/Bviridis/All_samples_populations/ordered_plinks") 

## All SNPs (about 80,000)
alldata <- read.PLINK('plink_fam_sorted.raw', chunkSize=1000, parallel = TRUE, n.cores=1, saveNbAlleles=T) ## If you are on linus you can change n.cores. But >1 core is not supported on windows.

## Plot for visualisation
sample_names <- alldata@ind.names

sample_names_ordered <- sample_names[order(sample_names, decreasing = T)]

ticks = seq(1,length(sample_names))

MyglPlot(alldata, posi = 'topleft', cex = 0.5, col = c("blue", "green", "orange"), axes = FALSE, ann = FALSE)
axis(2, at = ticks, labels = ordered(sample_names_ordered), las = 2)

#axes=FALSE, frame.plot=TRUE, ann = FALSE, yaxt = 'no')

X <- t(as.matrix(alldata))
X <- X[, ncol(X):1]

X <- as.data.frame(X)
X <- as.matrix(X[ , order(names(X))])

glPlot

MyglPlot <- function(x, col = NULL, legend = TRUE, posi = "bottomleft", 
          bg = rgb(1, 1, 1, 0.5), ylabels = NULL , ...) 
{

  X <- t(as.matrix(alldata))
  X <- X[, ncol(X):1]
  
  X <- as.data.frame(X)
  X <- as.matrix(X[ , order(names(X), decreasing = T)])
  
  ylabpos <- pretty(1:nInd(x), 5)
  
  if (is.null(col)) {
    myCol <- colorRampPalette(c("royalblue3", "firebrick1"))(max(X, 
                                                                 na.rm = TRUE) + 1)
  }
  else {
    myCol <- col
  }
  image(x = 1:nLoc(x), y = 1:nInd(x), z = X, xlab = "SNP index", 
        ylab = "Individual index", col = myCol, ...)
  
  if (legend) {
    legend(posi, fill = myCol, legend = 0:max(X, na.rm = TRUE), 
           horiz = TRUE, bg = bg, title = "Number of 2nd allele")
  }
  return(invisible())
}

image.default
