#install.packages("polysat")
library(polysat)
#install.packages("combinat")
library(combinat)
par(mar = c(4,4,4,4))
library(reshape2)
library(ggplot2)



?meandistance.matrix
?read.Structure

## read in from structure file - need to specify ploidy- ie number of lines in data file per indiv.
## if mixed ploidy in file- keep same number of lines for each indiv, and fill spare lines in lower ploidy samples with -9

data <- read.Structure("/Users/Dan/Dropbox/PhD/Dans_PhD_Shared/Triploid_gibel_hybrids/data/gibel_finland_NP.txt", ploidy = 4)
summary(data)

fil <- read.delim("/Users/Dan/Dropbox/PhD/Dans_PhD_Shared/Triploid_gibel_hybrids/data/gibel_finland_NP.txt", header = T)


## then after reading in, run the estimate ploidy function and then use summary to check/edit
data <- estimatePloidy(data)
summary(data)

# And if you want to use the "Bruvo" distance index, which apparently uses mutation models then
# you need to add the repeat length for each locus in as below . . . I imagine doing this would be more accurate
Usatnts(data) <- c(2,2,2,2,2,2)

PopNames(data)

## Then get genetic distances between samples

testmat <- meandistance.matrix(data)

## can now do a principle coordinate analysis on this
pca2 <- cmdscale(testmat, eig = TRUE, k = 5)
palette = c("black", "blue", "darkblue", "green", "turquoise", "magenta", "green", "red", "yellow", "purple") ## add as many cols a populations if thats what you want to colour by

cols = c(rep("pink", 26), rep("red",7), "yellow", rep("red", 8), "yellow", rep("red",2), "yellow", 
rep("red",33), "green", rep("blue", 8), rep("yellow", 31),"darkblue", "darkblue", "darkblue", rep("green", 9), rep("lightblue", 13), "black", "black", rep("purple", 4), "red")

colour_lables <- data.frame(cols, row.names(pca2$points))

row.names(pca2$points) <- colour_lables$cols

plot(pca2$points[,1], pca2$points[,2], pch = 21, col="black", bg = palette[colour_lables$cols], main = "PCA with Bruvo distance", cex = 2)
text(pca2$points[,1], pca2$points[,2], colour_lables$row.names.pca2., cex = 0.7)

## Note the colours aren't quite right - the GF population won't colour properly. Fixed this in Inkscape.




## I am not sure if the colours for population etc are correct. Make sure of this.
## Same for individual names



head(data)

length(pca2[,1])



# Can also get allele freqs in order to generate Fst matrices. . . but
#"For mixed ploidy systems, the simpleFreq function is available, but will be biased toward underestimating common allele
#frequencies and overestimating rare allele frequencies, which will cause an underestimation of FST" this isn'too much of a probelm for 
#us - just want to find diagnostic alleles.

freqs <- simpleFreq(data)
head(freqs)

## find alleles for which crucian have fixed diagnostic alleles. . . .

cru_sep_data <- read.Structure("/Users/Dan/Dropbox/PhD/Dans_PhD_Shared/Triploid_gibel_hybrids/data/gibel_crucian_separate_pop.txt", ploidy = 4)

cru_sep_data <- estimatePloidy(cru_sep_data)
summary(cru_sep_data)


freqs <- simpleFreq(cru_sep_data)
tfrqs <- t(freqs[,c(2:100)])

melted <- melt(tfrqs)

cols= c("grey", "grey", "grey", "grey", "grey","grey", "green", "red", "grey","grey")

ggplot(melted, aes(Var1,value,fill=Var2, width =2))+
  geom_bar(position="dodge",stat="identity")+
  scale_fill_manual(values= cols)+
  theme(axis.text.x = element_text(angle = 90))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))


barplot(freqs)

cru_diag_loci <- freqs[which(freqs$)] ## pop 99 is my crucian




fsts <- calcFst(freqs)


## Note that apparently we should be dealing with allotetraploid genomes separately, by coding the loci for which allele came from which genome . . 
# this should be relatively easy in our data . . . but right now i don't have time . . . 


