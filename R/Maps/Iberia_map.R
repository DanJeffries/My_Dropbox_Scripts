##ggmap

#install_github(ggmap") ## load devtools first
library(ggmap)
#install.packages("ggplot2")
library(ggplot2)

setwd("/home/djeffrie/Data/RADseq/Riberica/New_RAD_data/Maps/")

all_samples <- read.delim("pops3.txt", stringsAsFactors = F)

map <- get_map(location = c(-10,39,0,44),  zoom = 6, maptype = "terrain-background")

?get_map

mapPoints <- ggmap(map) +
  geom_point(aes(x = LONG, y = LAT, size = N), data = all_samples, alpha = .5) +
  geom_text(aes(x = LONG, y = LAT, label = POP, size = 4), data = all_samples, hjust = 1.5)

svg("All_R_iberca_samples.svg")
mapPoints
dev.off()

