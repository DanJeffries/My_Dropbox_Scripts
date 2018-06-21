library(ggplot2)
library(grid)
library(diveRsity) ## contains the multiplot function for ggplot2

setwd("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/structure")

str_outs <- read.delim("Str_assignments_formatted.out")
str_samples <- read.delim("../NAMES_final.txt", header = F)
str_pops <- read.delim("../pops.lab", header = F)

## Annotate data properly
str_outs$pop <- str_pops$V1
str_outs$sample <- str_samples$V1
head(str_outs)
tail(str_outs)

str_melted1 <- melt(str_outs[c(1:98),])
str_melted2 <- melt(str_outs[c(99:207),])
str_melted3 <- melt(str_outs[c(208:303),])
str_melted4 <- melt(str_outs[c(304:439),])
str_melted5 <- melt(str_outs[c(440:544),])
str_melted6 <- melt(str_outs[c(545:653),])
str_melted7 <- melt(str_outs[c(661:694),])
str_melted8 <- melt(str_outs[c(695:799),])
str_melted9 <- melt(str_outs[c(800:903),])
str_melted10 <- melt(str_outs[c(904:1002),])
str_melted11 <- melt(str_outs[c(1003:1183),])


p1 <- ggplot(str_melted1, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

p2 <- ggplot(str_melted2, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p3 <- ggplot(str_melted3, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p4 <- ggplot(str_melted4, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p5 <- ggplot(str_melted5, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p6 <- ggplot(str_melted6, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p7 <- ggplot(str_melted7, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p8 <- ggplot(str_melted8, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p9 <- ggplot(str_melted9, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p10 <- ggplot(str_melted10, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p11 <- ggplot(str_melted11, aes(x = sample, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('sample',values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))



multiplot(p1, p2, p3, p4, p5, cols=1)

multiplot(p6, p7, p8, p9, p10, p11, cols=1)
