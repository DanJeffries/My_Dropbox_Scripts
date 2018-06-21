library(raster)
library(maptools)
library(rgdal)
library(calibrate)
library(sp)

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV_CLOUD_COMP/Climate_data/2.5_arc_minutes/")

par(mfrow = c(1,1))
par(mar = c(2,2,2,2))

## I downloaded these files maunally into the dir above, but there is a getData function which will automatically get data from worldclim:
# e.g. w = getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)


## Can convert the binary files (.bil) to GeoTiffs to read in using the rgdal and raster
#input_name=raster("tmax10.bil") #if you use 
#output_name=("tmax_10.tif")
#writeRaster(input_name, output_name,format="GTiff",datatype='INT1U',overwrite=TRUE)



## But actually you don't need to, it will read .bil files, despite the vingette not mentioning this

clip.extent <- as(extent(-10, 45,43,70), "SpatialPolygons")


## make raster objects from the files for each month and crop them to EU only
## the crop function is very useful, allows you to crop to the long and lat limits you want

EU.APR <- crop(raster("Temperature Data/tmean4.bil"), clip.extent)
EU.MAY <- crop(raster("Temperature Data/tmean5.bil"), clip.extent)
EU.JUN <- crop(raster("Temperature Data/tmean6.bil"), clip.extent)
EU.JULY <- crop(raster("Temperature Data/tmean7.bil"), clip.extent)
EU.AUG <- crop(raster("Temperature Data/tmean8.bil"), clip.extent)
EU.SEPT <- crop(raster("Temperature Data/tmean9.bil"), clip.extent)

## The temperature data comes in 12 files, one for each month, so I can average the data across these months and the result should be the average of the data at each cell....

EU.Apr_Sept <- overlay(EU.APR, EU.MAY, EU.JUN, EU.JULY, EU.AUG, EU.SEPT, fun = mean)

plot(EU.JULY, main = "Average European temperature for Apr") ## plots the average temp for these months


plot(EU.Apr_Sept, main = "Average European temperature for Apr - Sept") ## plots the average temp for these months

## Looks good

#### Now to Put my RAD points on there #######


RADpoints <- read.csv("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/RAD seq/Samples/RAD sample coordinates.csv", header = T)


## Subset samples into their project allocations.

Ipoints <- RADpoints[c(1:8),]
NBTpoints <- RADpoints[c(10:11,15:16,22),]
BTpoints <- RADpoints[c(9,12:14),]
NBBpoints <- RADpoints[c(19,21,24,26),]
BBpoints <- RADpoints[c(20,23,25,27),]
SEUpoints <- RADpoints[c(18,28),]
Goldpoints <-RADpoints[c(29,30),]
Gibpoints <-RADpoints[c(31,33),]
Commpoints <-RADpoints[14,]


## Plot Map


#map("worldHires", xlim=c(-10, 45), ylim=c(44,72), col="gray90", fill=TRUE)

##plot filled points

points(Ipoints$lon, Ipoints$lat, pch = 16, col = "black")
points(NBTpoints$lon, NBTpoints$lat, pch = 15, col = "blue")
points(BTpoints$lon, BTpoints$lat, pch = 17, col = "blue")
points(NBBpoints$lon, NBBpoints$lat, pch = 15, col = "red")
points(BBpoints$lon, BBpoints$lat, pch = 17, col = "red")
points(SEUpoints$lon, SEUpoints$lat, pch = 17, col = "darkgreen")
points(Goldpoints$lon, Goldpoints$lat, pch = 18, col = 'gold')
points(Gibpoints$lon, Gibpoints$lat, pch = 18, col = 'gray')
points(Commpoints$lon, Commpoints$lat, pch = 18, col = 'purple')

textxy(RADpoints$lon, RADpoints$lat, RADpoints$population, col = 'black', cex = 0.6, pos = 2) ## add pop names

#### Add Keys ####


points(-10, 65.01, pch = 16, col = "black")
textxy(-1, 65, "Introgression samples", pos = 2)

points(-10, 64, pch = 15, col = "blue")
textxy(1.2, 64, "Non-Bottlenecked, Temperate", pos = 2)

points(-10, 67, pch = 17, col = "blue")
textxy(-0.4, 67, "Bottlenecked Temperate", pos = 2)

points(-10, 66, pch = 15, col = "red")
textxy(-0.15, 65.95, "Non-Bottlenecked, Boreal", pos = 2)

points(-10, 68, pch = 17, col = "red")
textxy(-1.8, 68, "Bottlenecked Boreal", pos = 2)

points(-10, 63, pch = 17, col = "darkgreen")
textxy(-2.9, 63, "Southern Europe", pos = 2)

points(-10, 69, pch = 18, col = 'gold')
textxy(-4, 69, "Pure goldfish", pos = 2)

points(-10, 70, pch = 18, col = 'gray')
textxy(-5, 70, "Pure gibel", pos = 2)

points(-10, 71, pch = 18, col = 'purple')
textxy(-3.8, 71, "Pure common", pos = 2)

#Add outline to points

points(NBTpoints$lon, NBTpoints$lat, pch = 0, col = "black")
points(BTpoints$lon, BTpoints$lat, pch = 2, col = "black")
points(NBBpoints$lon, NBBpoints$lat, pch = 0, col = "black")
points(BBpoints$lon, BBpoints$lat, pch = 2, col = "black")


## Map sorted! ##

## So how do I get to the data? ##

# The easiest way was to do this...
Av.Temperature_data <- rasterToPoints(EU.Apr_Sept)

## Basically gives a matrix. 

## i put this into python to search for the cells at the right coordinates . . . 

pop_temps <- read.delim("./Temp_at_pop_coords_take_2.tsv", sep = "\t", header = F)
names(pop_temps) <- c("pop", "temp")
head(pop_temps)


pop_values <- split(pop_temps$temp, pop_temps$pop)
pop_means <- lapply(pop_values, mean)
pop_means
pop_means_temp<- as.matrix(pop_means)
pop_means_temp

write.csv(pop_means_temp, "./pop_mean_temps.csv")
means_tmep <- read.csv("./pop_mean_temps.csv", header = F)
names(means_tmep) <- c("pop", "temp")
means_tmep
means_tmep$temp <- means_tmep$temp/10

## Join it to RADcoords file

RADcoords_and_Av_temp <- cbind(RADcoordsall, means_tmep$temp)

## Quick plot to see what the relationship of lat and temp is in the samples...
plot(RADcoords_and_Av_temp$lat, RADcoords_and_Av_temp$means_tmep, pch = 16, col = "blue")
textxy(RADcoords_and_Av_temp$lat, RADcoords_and_Av_temp$means_tmep, RADcoords_and_Av_temp$population, cex = 0.7)


### Only Pure Cru RAD pops ####

Pure_sub_coords_and_temps <- RADcoords_and_Av_temp[c(1,3:5,15,17:21,23:30,32),]

plot(Pure_sub_coords_and_temps$lat, Pure_sub_coords_and_temps$means_tmep, pch = 16, col = "blue")
textxy(Pure_sub_coords_and_temps$lat, Pure_sub_coords_and_temps$means_tmep, Pure_sub_coords_and_temps$population, cex = 0.7)

## Highlight the bottlenecked and non_bottlenecked pops ##
row.names(Pure_sub_coords_and_temps) <- seq(1,19,1)
bottlenecked <- Pure_sub_coords_and_temps[c(14,3,6,16,5,2,11,1),]
non_bottlenecked <- Pure_sub_coords_and_temps[c(12,7,19,15,9,17,13,4,8),]

points(bottlenecked$lat, bottlenecked$"means_tmep$temp", col = "red", pch = 16)

### Quick regression ## 

abline(lm(non_bottlenecked$"means_tmep$temp"~non_bottlenecked$lat), col = "blue")
abline(lm(bottlenecked$"means_tmep$temp"~bottlenecked$lat), col = "red")

summary(lm(non_bottlenecked$"means_tmep$temp"~non_bottlenecked$lat))

#Call:
#  lm(formula = non_bottlenecked$"means_tmep$temp" ~ non_bottlenecked$lat)
#
#Residuals:
# Min       1Q   Median       3Q      Max 
#-0.66584 -0.08193  0.05862  0.16149  0.57394 
#
#Coefficients:
  #Estimate Std. Error t value Pr(>|t|)
#(Intercept)          34.35546    2.08428   16.48 7.38e-07
#non_bottlenecked$lat -0.37358    0.03624  -10.31 1.75e-05
#
#(Intercept)          ***
  #non_bottlenecked$lat ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.3866 on 7 degrees of freedom
#Multiple R-squared:  0.9382,  Adjusted R-squared:  0.9294 
#F-statistic: 106.2 on 1 and 7 DF,  p-value: 1.752e-05

summary(lm(bottlenecked$"means_tmep$temp"~bottlenecked$lat))

#Call:
#  lm(formula = bottlenecked$"means_tmep$temp" ~ bottlenecked$lat)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.9500 -0.3731 -0.1072  0.4713  1.0553 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      33.32675    2.43076  13.710 9.36e-06 ***
#  bottlenecked$lat -0.36835    0.04193  -8.784 0.000121 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.7111 on 6 degrees of freedom
#Multiple R-squared:  0.9279,  Adjusted R-squared:  0.9158 
#F-statistic: 77.16 on 1 and 6 DF,  p-value: 0.0001207

textxy(57,18, "Non-bottlenecked: adj. R^2 = 0.93, p < 0.001", col = "blue", cex = 1)
textxy(57,17.5, "Bottlenecked: adj. R^2 = 0.92, p < 0.001", col = "red", cex =1)
title(main = "Correlation between latitude and temperature in Pure RAD populations")
