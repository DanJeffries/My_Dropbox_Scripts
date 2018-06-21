setwd("~/Desktop/WindeRmere/")

install.packages("RgoogleMaps") ## install google map package if its not already installed


library(RgoogleMaps)
library(calibrate)


## Define map property objects ##

lat <- c(54.23707, 54.446842 ) #define our map's ylim
lon <- c(-2.9, -2.85) #define our map's xlim
center = c(mean(lat), mean(lon))  #tell what point to center on
zoom <- 12  #zoom: 1 = furthest out (entire globe), larger numbers = closer in

### Get coordinates for all sites and each species separately ###

all_sites <- read.csv("./Coordinate_files/Sample_coords.csv")

perch <- read.csv("./Coordinate_files/perch_coords_only.csv")
pike <- read.csv("./Coordinate_files/pike_coords_only.csv")
roach <- read.csv("./Coordinate_files/roach_coords_only.csv")
charr <- read.csv("./Coordinate_files/charr_coords_only.csv")


### Plotting separate maps for each spp ## 

par(mfrow = c(1,4))

Wind_map <- GetMap(center= c(center), size = c(200, 490), zoom=zoom, maptype="terrain",SCALE = T, frame=T)

## Perch ##
PlotOnStaticMap(Wind_map, lat=all_sites$lat,lon=all_sites$lon, destfile="terrclose.png", pch = 1, cex = 1)
PlotOnStaticMap(Wind_map, lat=perch$lat,lon=perch$lon, destfile="terrclose.png", pch = 16, cex = 1, add = T)
textxy(-1, 225, "Perch", cex = 2)

## Pike ##
PlotOnStaticMap(Wind_map, lat=all_sites$lat,lon=all_sites$lon, destfile="terrclose.png", pch = 1, cex = 1)
PlotOnStaticMap(Wind_map, lat=pike$lat,lon=pike$lon, destfile="terrclose.png", pch = 16, cex = 1, add = T)
textxy(-1, 225, "Pike", cex = 2)

## Roach ##
PlotOnStaticMap(Wind_map, lat=all_sites$lat,lon=all_sites$lon, destfile="terrclose.png", pch = 1, cex = 1)
PlotOnStaticMap(Wind_map, lat=roach$lat,lon=roach$lon, destfile="terrclose.png", pch = 16, cex = 1, add = T)
textxy(-1, 225, "Roach", cex = 2)

## Charr ##
PlotOnStaticMap(Wind_map, lat=all_sites$lat,lon=all_sites$lon, destfile="terrclose.png", pch = 1, cex = 1)
PlotOnStaticMap(Wind_map, lat=charr$lat,lon=charr$lon, destfile="terrclose.png", pch = 16, cex = 1, add = T)
textxy(-1, 225, "Charr", cex = 2)



## Plotting Perch and Charr on same map ## 

par(mfrow = c(1,1), mar = c(0.1,0.1,0.1,0.1))

## All sites ##
PlotOnStaticMap(Wind_map, lat=all_sites$lat,lon=all_sites$lon, destfile="terrclose.png", pch = 1, cex = 1)
### Perch
PlotOnStaticMap(Wind_map, lat=perch$lat,lon=perch$lon, destfile="terrclose.png", pch = 16, cex = 1, add = T)
## Charr ##
PlotOnStaticMap(Wind_map, lat=charr$lat,lon=charr$lon, destfile="terrclose.png", pch = 16, cex = 0.7, col = "red", add = T)
title(main = "Perch and Charr presence absence")




