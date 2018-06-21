library(maps)
library(mapplots)
library(mapdata)
#install.packages("calibrate")
library(calibrate)
par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(7,7)) ## set plot window size

## Maps


setwd("~/../Dropbox/PhD/Dans_PhD_Shared/Triploid_gibel_hybrids/")

whole <- map("worldHires", xlim=c(5, 30), ylim=c(50,65), col="gray90", fill=TRUE)
cords <- read.csv("data//coordinates.txt", sep = "\t")

points(cords$lon, cords$lat, pch = 16, col = "black")

FINlay <- map("worldHires", xlim=c(22.5, 26), ylim=c(59.2,60.7), col="gray90", fill=TRUE)
