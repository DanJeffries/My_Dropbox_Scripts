library(raster)
library(maptools)
library(rgdal)
library(calibrate)
library(rgeos)

setwd("/home/dan/Desktop/River_systems/")

par(mfrow = c(1,1))
par(mar = c(2,2,1,1))

## Using readOGR function to read in the files downloaded from HydroSHEDS: http://hydrosheds.cr.usgs.gov/data.php 

BAS<- readOGR(".", layer = "eu_bas_15s_beta")

## Here the folder containing the files is the first arg, the second is the prefix of the files containing the basin data

#BASEU <- crop(BAS, extent(-10, 45,43,70)) 
clip.extent <- as(extent(-10, 45,43,70), "SpatialPolygons")

# This very quick method keeps whole countries
gI <- gIntersects(BAS , clip.extent , byid = TRUE )
Europe <- BAS[ which(gI) , ]
plot(Europe)


#If you want to crop the country boundaries, it's slightly more involved:
# This crops countries to your bounding box

out <- lapply( which(gI) , function(x){ gIntersection( BAS[x,] , clip.extent ) } )

# But let's look at what is returned
table( sapply( out , class ) )


# We want to keep only objects of class SpatialPolygons                 
keep <- sapply(out, class)
out <- out[keep == "SpatialPolygons"]


# Coerce list back to SpatialPolygons object
Europe <- SpatialPolygons( lapply( 1:length( out ) , function(i) { Pol <- slot(out[[i]], "polygons")[[1]];   slot(Pol, "ID") <- as.character(i)
                                                                   Pol
}))

plot( Europe )