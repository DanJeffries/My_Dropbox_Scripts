detach("package:hierfstat")
library(adegenet)

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DIY ABC/UK analyses - Final poolings/")
##Read in with adegenet to get some information of populations
body <- read.fstat("FSTAT_seppop_named.dat") ## Load file - converts to GENIND object.

pop_names <- names(table(pop(body))) ## correct order - i checked

data("nancycats")

nancycats
body


summary(nancycats)
summary(body)

genpop_body <- genind2genpop(body)

genpop_body@tab

summ_genpop <- summary(genpop_body)

