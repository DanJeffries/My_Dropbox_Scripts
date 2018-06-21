library(calibrate)


##### Preliminary run, all populations ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/preliminary_run")

BAYESFACTS <- read.delim("bf_environ.pop_mean_stdzd_temps.env", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "Preliminary run, all pops")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY_prelim <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY_prelim, 100)

#plot(head(BAY_prelim$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY_prelim$V2, 100), head(BAY_prelim$V1, 100), cex = 1, offset = 0.3)





## non bottled rep1 ####

setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_1/")
BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7)

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)

##### Non-bottled rep_2 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_2/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7)

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)


##### Non-bottled rep_3 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_3/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7)

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)


##### Non-bottled rep_4 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_4/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7)

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)


##### Non-bottled rep_5_1000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_5_100000/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7)

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)

### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY5 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY5)

#plot(head(BAY5$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY5$V2, 100), head(BAY5$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY5, 100),"./top100_BFs.csv")


##### Non-bottled rep_6_1000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_6_100000/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 6")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)

### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY6 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY6)

#plot(head(BAY6$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY6$V2, 100), head(BAY6$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY6, 100),"./top100_BFs.csv")


##### Non-bottled rep_7 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_7_100000/")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 7_100000")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY7 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY7, 100)

#plot(head(BAY7$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY7$V2, 100), head(BAY7$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY7, 100),"./top100_BFs.csv")



##### Non-bottled rep_8_500000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_8_500000//")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 8_500000")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY8 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY8, 100)

#plot(head(BAY7$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY7$V2, 100), head(BAY7$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY8, 100),"./top100_BFs.csv")


##### Non-bottled rep_9_500000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_9_500000//")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 9_500000")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY9 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY8, 100)

#plot(head(BAY9$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY7$V2, 100), head(BAY7$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY9, 100),"./top100_BFs.csv")



##### Non-bottled rep_10_500000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_10_500000//")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 10_500000")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY10 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY10, 100)

#plot(head(BAY9$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY7$V2, 100), head(BAY7$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY10, 100),"./top100_BFs.csv")



##### Non-bottled rep_11_500000 ########


setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/non_bottlenecked/rep_11_500000//")

BAYESFACTS <- read.delim("bf_environ.non_bottled_temp_stdzd", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2, cex = 0.7, main = "rep 11_500000")

abline(h = ((mean(BAYESFACTS$V2))+2*(sd(BAYESFACTS$V2))), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.5, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)



### ordering BFs

ord <- order(BAYESFACTS[,2]) ## creates a vector of line numbers (row.names) in order from smalles value of V2 to highest

BAY11 <- BAYESFACTS[rev(ord),] ## "subsets" (but doesn't lose any data) in the reverse order of "ord" (to get largest to smallest)
head(BAY11, 100)

#plot(head(BAY11$V2, 100), pch = 16, cex = 1)
#textxy(seq(1,100,1), head(BAY11$V2, 100), head(BAY11$V1, 100), cex = 1, offset = 0.3)

write.csv(head(BAY11, 100),"./top100_BFs.csv")




#### Bottlenecked ####

setwd("/home/dan/RAD_programs/BAYENV/bayenv_release/bottlenecked/")


BAYESFACTS <- read.delim("bf_environ.bottlenecked_stdized_mean_temps.", header = F)

head(BAYESFACTS)
BAYESFACTS$V3 <- NULL

plot(BAYESFACTS$V2[1:5000], cex = 0.7)

abline(h = (((max(BAYESFACTS$V2))/100)*90), col = "blue", lty = 3)

textxy(row.names(BAYESFACTS), BAYESFACTS$V2, BAYESFACTS$V1, offset = -0.6, cex = 1)

mean(BAYESFACTS$V2)
sd(BAYESFACTS$V2)
((max(BAYESFACTS$V2))/100)*90
