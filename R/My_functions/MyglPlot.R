## Plots the adegenet glplot but allows lables for the samples. (Don't know why they didn't do this before!!)
## Sorts the samples by descending alphabetical order. Prepend sample names if order needed is different!

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
