library(maps)

install.packages("mapdata")
library(mapdata)
install.packages("mapplots")
library(mapplots)

par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(6,6)) ## set plot window size
mycol = rainbow(12, alpha = 0.6) # set colours

map("worldHires", xlim=c(-100, -40), ylim=c(10,50), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)

map("worldHires", xlim=c(-86, -80), ylim=c(28,34), col="gray90", fill=TRUE)##plots the area of the map that I am interested in (just Europe, leave out x/ylim arguments for whole world view)

coords <- read.delim("/home/djeffrie/Data/RADseq/Lsphenosephalus/Sample_info/Sample_coords.txt", header = F)

points(coords$V3, coords$V2, pch=19, col="red", cex=0.9)
points(coords$V3, coords$V2, col="black", cex=0.91)
text(coords$V3, coords$V2, coords$V1, cex = 0.7)
