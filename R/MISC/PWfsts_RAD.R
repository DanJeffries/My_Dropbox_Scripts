install.packages("diveRsity")
library(diveRsity)
vignette("diveRsity")
library(parallel)
install.packages("doParallel")
install.packages("foreach")
install.packages("iterators")
install.packages("StAMPP")
library(StAMPP)
library(adegenet)
help("StAMPP")
library(hierfstat)
setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Complete dataset outputs/")

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17_Ho_filetered")

### STaMPP #####

## read in plink using adegenet ##
?pp.fst
RADdata <- read.PLINK('plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
RADdata_df <- read.fstat('batch_1_FSTAT.dat')
RADdata <- df2genind(RADdata_df)

pops <- read.delim("../populations.txt", header = F)

pop(RADdata) <- pops$V1

RADdata_sep <- seppop(RADdata)
names(RADdata_sep)

levels(pop(RADdata))

allele_SUMs <- data.frame(BF = glSum(RADdata_sep$BF), BOR = glSum(RADdata_sep$BOR), CAKE = glSum(RADdata_sep$CAKE), CALK = glSum(RADdata_sep$CALK), COP = glSum(RADdata_sep$COP), MOAT = glSum(RADdata_sep$MOAT), OBY = glSum(RADdata_sep$OBY), OU = glSum(RADdata_sep$OU), PED = glSum(RADdata_sep$PED), POLEN = glSum(RADdata_sep$POLEN), PRO = glSum(RADdata_sep$PRO), SD = glSum(RADdata_sep$SD), SK = glSum(RADdata_sep$SK), STEC = glSum(RADdata_sep$STEC), STYV = glSum(RADdata_sep$STYV), TROM = glSum(RADdata_sep$TROM), TU = glSum(RADdata_sep$TU), V = glSum(RADdata_sep$V), WEN = glSum(RADdata_sep$WEN))
sum(allele_SUM$Allele_counts)
Allele_counts <- data.frame(Pop_allele_counts = colSums(allele_SUMs))
mean(Allele_counts$Pop_allele_counts)


RADgenp <- genind2genpop(RADdata, missing ="chi2")
RADgenp$pop.names

RADdata$loc.names

RADhier <- genind2hierfstat(RADdata)

RADdata$tab[,c(1:10)]

RADagain <- df2genind(RADdata$tab[,c(1:10)])


RADdata_meanNAs <- na.replace(RADdata$tab, "chi2")

?na.replace

head(RADhier[c(1:20)])

PWfsts <- pp.fst(RADhier)

stamppRAD <- stamppConvert(RADdata, type = "genlight") ## convert to stampp format

PWfsts <- stamppFst(stamppRAD, nboots = 1, percent = 95, nclusters = 1) ## calculate PWfsts
?stamppFst


####### diveRsity ########

data(Test_data)
system.time({
  pwDiff <- fastDivPart(Test_data, pairwise = TRUE,
                        bs_pairwise = TRUE, fst = TRUE,
                        boots = 1000, para = TRUE)
}) 

pwDiff <- fastDivPart(Test_data, pairwise = TRUE,
                      bs_pairwise = TRUE, fst = TRUE,
                      boots = 1000, para = TRUE)





install.packages("hierfstat")
library(hierfstat)

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Complete dataset outputs/")



All_micros <- read.fstat("PurecruplusManc_stripped_and_checked_numbered.DAT")
All_RAD_micros <- All_micros[c(14:40,281:334,342:378,537:557,583:642,654:672,683:726,744:825),]

pop_names <- c("CAKE", "TU", "POLEN", "MOAT", "BF", "SD", "STEC", "STYV", "SK", "CALK", "OU", "PRO", "COP", "OBY","PED","TROM", "WEN")
length(pop_names)
row.names(All_RAD_MicroFSTs) <- pop_names
names(All_RAD_MicroFSTs) <- pop_names
All_RAD_MicroFSTs ## OK named - check that they are in right order though!
RAD_micros_vect <- as.vector(t(All_RAD_MicroFSTs))
RAD_micros_vect

header_1 <- rep(pop_names, 16)
header_2 <- paste(c(rep("CAKE", 17), rep("TU", 17), rep("POLEN", 17), rep("MOAT", 17), rep("BF", 17), rep("SD", 17), rep("STEC", 17), rep("STYV", 17), rep("SK", 17), rep("CALK", 17), rep("OU", 17), rep("PRO", 17), rep("COP", 17), rep("OBY", 17), rep("PED", 17), rep("TROM", 17), rep("WEN", 17)))

final_micro_header <- paste(header_1, header_2, sep = "_")
final_micro_header

Micro_RAD_fsts_vect_labelled <- data.frame(final_micro_header, RAD_micros_vect)
Micro_RAD_fsts_vect_labelled

length(All_RAD_micros$Pop)
table(All_micros$Pop)

PWfsts <- pp.fst(All_RAD_micros)

write.table(PWfsts$fst.pp, "PW_fsts_npop49", sep = "\t")

###### RAD ############

RAD_Fsts <- read.table("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17_Ho_filetered/batch_1.fst_summary.tsv", sep = "\t", header = T, row.names = 1)
?as.vector

RAD_fst_vect <- as.vector(t(RAD_Fsts))
RAD_fst_vect

### make pop headers to go with ##

RAD_pop_names <- names(RAD_Fsts)

top_header <- rep(RAD_pop_names,18)
sub_header <- paste(c(rep("BF", 19),rep("BOR", 19),rep("CAKE", 19),rep("CA.LK", 19),rep("COP", 19),rep("MOAT", 19),rep("OBY", 19),rep("OU", 19),rep("PED", 19),rep("POLEN", 19),rep("PRO", 19),rep("SD", 19),rep("SK", 19),rep("STEC", 19),rep("STYV", 19),rep("TROM", 19),rep("TU", 19),rep("V", 19)))

Final_header <- paste(top_header, sub_header, sep = "_")

RAD_fst_vect_labelled <- data.frame(Final_header, RAD_fst_vect) ## this is the final vecter version of the matrix
RAD_fst_vect_labelled

### Comparing vectors ####

RAD_ordered <- RAD_fst_vect_labelled[order(Final_header),]

Micro_ordered <- Micro_RAD_fsts_vect_labelled[order(final_micro_header),]

length(RAD_ordered$Final_header)
length(Micro_ordered$final_micro_header)


#### write these to files so they can be edited in Bash to remove the pops in the RAD dataset that aren't in the Microsats

write.csv(RAD_ordered, "RAD_PWfsts_vectorised", sep = "\t")
write.csv(Micro_ordered, "Micro_PWfsts_vectorised", sep = "\t")

RAD_formatted <- read.delim("RAD_PWfsts_formatted_final.tsv")
Micro_formatted <- read.delim("Micro_PWfsts_formatted_no_NAs_sorted.tsv")
head(Micro_formatted, 20)
head(RAD_formatted, 20)
length(RAD_formatted$Pop)
length(Micro_formatted$final_micro_header)

RAD_Micro_PWfsts <- data.frame(Pop = RAD_formatted$Pop, RAD = RAD_formatted$PW_fst, Micro = Micro_formatted$RAD_micros_vect)
write.table(RAD_Micro_PWfsts, "~/Dropbox/PhD/Dans_PhD_Shared/Papers/Phylogeography paper/Final_data_files/RAD_M3_PW_fsts.tsv", sep = "\t")

RAD_Micro_PWfsts.melted <- melt(RAD_Micro_PWfsts)
library(ggplot2)
library(reshape2)
gplot <- ggplot(RAD_Micro_PWfsts.melted, aes(Pop, value, fill = variable)) + geom_bar(position="dodge", stat = "identity")
gplot + theme(axis.text.x=element_text(angle=-90))

Test <- cor.test(RAD_Micro_PWfsts$RAD, RAD_Micro_PWfsts$Micro)

Test

### Write to file ##

All_RAD_MicroFSTs <- read.table("PW_fsts_npop_19", sep = "\t")

All_RAD_MicroFSTs

###########################################
par(mfrow=c(1,1))
par(mar = c(4,4,4,4))

RAD_Micro_PWfsts
RAD_Micro_PWfsts$Micro_minus_RAD <- (RAD_Micro_PWfsts$Micro + RAD_Micro_PWfsts$RAD)
RAD_Micro_PWfsts_ordered <- RAD_Micro_PWfsts[order(RAD_Micro_PWfsts$Micro_minus_RAD),]

FSTbar <- barplot(RAD_Micro_PWfsts_ordered$Micro_minus_RAD)
axis(1, RAD_Micro_PWfsts_ordered$Pop, at = FSTbar, las = 2, cex.axis = 0.5)
?axis
hets <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Data/Micro_RAD_comp/Micro_vs_RAD_Hobs.txt")

hets$Micro_minus_RAD <- hets$Micro - hets$RAD
bar <- barplot(hets$Micro_minus_RAD)
axis(1, hets$Pop, at = bar)
title(main = "Microsats Ho - RAD Ho")
