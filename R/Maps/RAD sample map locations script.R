library(maps)
library(mapplots)
library(mapdata)
library(calibrate)
par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(7,7)) ## set plot window size


RADpoints <- read.csv("~/Dropbox/PhD/Dan's PhD (Shared)/RAD seq/Samples/RAD sample coordinates.csv", header = T)
RADpoints <- read.csv("~/Dropbox/PhD/Dans_PhD_Shared/RAD seq/Samples/RAD sample coordinates.csv", header = T)

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


map("worldHires", xlim=c(-10, 45), ylim=c(44,72), col="gray90", fill=TRUE)

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

