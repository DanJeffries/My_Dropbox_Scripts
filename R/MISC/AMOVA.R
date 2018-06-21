
################## MICROSATS ##########################

M1 <- read.fstat("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/Purecruchecked.DAT")

MicroCoords <- read.csv("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv",header = TRUE, row.names = 1) ##table containing xy coordinates

############# All populations (M1) #####################

M1gen <- genind2genpop(M1)

M1popnames <- pop(M1) ## get names from full microsat
M1_Hap1_noSTECnames <- popnames[c(1:9,11,14:49)] ## get only M1 hap 1 and non Danubian names
M3names <- levels(pop(RADsub)) ## get names for RAD and micro only


##### AMOVA #######

M1Dgen <- dist.genpop(M1gen,method=2) ##Generates the genetic distance semi-matrix

M1amova <- amova(M1Dgen~M1popnames)
