library(maps)
#install.packages("mapdata")
library(mapdata)
#install.packages("mapplots")
library(mapplots)

setwd("/home/djeffrie/Data/RADseq/Rchensinensis/Sample_info")

par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(6,6)) ## set plot window size
mycol = rainbow(12, alpha = 0.6) # set colours

coords <- read.delim("Coordinates.txt", header = F)

## Zoomed out
map("worldHires", xlim=c(103, 107), ylim=c(30,34), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(coords$V3, coords$V2, pch=19, col="red", cex=0.9)



## Zoomed in
map("worldHires", xlim=c(90, 130), ylim=c(20,50), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)
points(coords$V3, coords$V2, pch=19, col="red", cex=0.9)



points(coords$V3, coords$V2, pch=19, col="red", cex=0.9)
points(coords$V3, coords$V2, col="black", cex=0.91)
text(coords$V3, coords$V2, coords$V1, cex = 0.7)