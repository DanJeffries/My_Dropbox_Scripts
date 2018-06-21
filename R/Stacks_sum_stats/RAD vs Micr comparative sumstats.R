
#### RAD summary #####

setwd("/media/dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_only/populations_r07_p17_m8/")

## made the PLINK input files using vcftools from the batch_1.vcf file outputted by stacks

RAD <- read.genepop("batch_1_single_snp.gen")
RAD

RADsummary <- summary(RAD)

mean(RADsummary$pop.eff)
sd(RADsummary$pop.eff)
mean(RADsummary$loc.nall)
sd(RADsummary$loc.nall)
mean(RADsummary$pop.nall)
sd(RADsummary$pop.nall)

length(RADsummary$loc.nall)
length(RADsummary$)

###### Full Microsat dataset ##########

microsats<- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/PurecruplusManc_stripped_and_checked.DAT") ## Load file - converts to GENIND object.

alleles(microsats)

microsats_summary <- summary(microsats)
microsats_summary$pop.nall


###### Summary stats fig ##########
par(mfrow=c(2,2))
allele_bar <- barplot(microsats_summary$pop.nall,main="Number of Alleles per population",ylab="Number of genotypes",las=3, xaxt = 'n')
axis(1, levels(pop(microsats)), at = allele_bar, las = 2, cex.axis = 0.8)
abline(h = mean(microsats_summary$pop.nall), col = "red", lty = 2)
loc_bar <- barplot(microsats_summary$loc.nall,ylab="Number of alleles", main="Number of alleles per locus", xaxt = 'n')
axis(1, microsats$loc.names, at = loc_bar, cex.axis = 0.8, las = 2)
abline(h = mean(microsats_summary$loc.nall), col = "red", lty = 2)
barplot(microsats_summary$Hexp-microsats_summary$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
abline(h = mean(microsats_summary$Hexp-microsats_summary$Hobs), col = "red", lty = 2)
barplot(microsats_summary$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)
abline(h = mean(microsats_summary$pop.eff), col = "red", lty = 2)


mean(microsats_summary$pop.eff)
sd(microsats_summary$pop.eff)
mean(microsats_summary$loc.nall)
sd(microsats_summary$loc.nall)
mean(microsats_summary$pop.nall)
sd(microsats_summary$pop.nall)

########## RAD indivs only ###########

RAD_only_microsats <- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_genotypes_1.dat") ## Load file - converts to GENIND object.
### this file contains only the RAD individuals! This should match up exactly with the RAD data file - ie: 170 individuals, all from the same pops. 
## note the amount of missing data is quite bad in some pops. 

RAD_only_microsats_pops <- read.delim("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/RAD/RADsample_microsat_POPS.txt", header = F)
RAD_only_microsats_pops
pop(RAD_only_microsats) 

pop(RAD_only_microsats) <- RAD_only_microsats_pops$V1

obj<- seppop(RAD_only_microsats)
obj

RAD_only_microsats_subnoSTEC <- repool(obj$BF, obj$CAKE , obj$TU , obj$POLEN , obj$MOAT, obj$SD, obj$STYV, obj$SK, obj$'CA-LK', obj$OU, obj$PRO, obj$COP, obj$OBY, obj$PED, obj$TROM, obj$WEN)

pop(RAD_only_microsats_subnoSTEC)

RAD_only_summary <- summary(RAD_only_microsats_subnoSTEC)



###### Summary stats fig ##########
par(mfrow=c(2,2))
allele_bar <- barplot(RAD_only_summary$pop.nall,main="Number of Alleles per population",ylab="Number of genotypes",las=3, xaxt = 'n')
axis(1, levels(pop(RAD_only_microsats_subnoSTEC)), at = allele_bar, las = 2, cex.axis = 0.8)
abline(h = mean(RAD_only_summary$pop.nall), col = "red", lty = 2)
loc_bar <- barplot(RAD_only_summary$loc.nall,ylab="Number of alleles", main="Number of alleles per locus", xaxt = 'n')
axis(1, RAD_only_microsats_subnoSTEC$loc.names, at = loc_bar, cex.axis = 0.8, las = 2)
abline(h = mean(RAD_only_summary$loc.nall), col = "red", lty = 2)
barplot(RAD_only_summary$Hexp-RAD_only_summary$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
abline(h = mean(RAD_only_summary$Hexp-RAD_only_summary$Hobs), col = "red", lty = 2)
barplot(RAD_only_summary$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)
abline(h = mean(RAD_only_summary$pop.eff), col = "red", lty = 2)

mean(RAD_only_summary$pop.eff)
sd(RAD_only_summary$pop.eff)
mean(RAD_only_summary$loc.nall)
sd(RAD_only_summary$loc.nall)
mean(RAD_only_summary$pop.nall)
sd(RAD_only_summary$pop.nall)



##### All indivs in RAD pops ###############
microsats<- read.fstat("/home/dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/PurecruplusManc_stripped_and_checked.DAT") ## Load file - converts to GENIND object.
### this file contains all the microsat data!

### Subset! ##### (NO STEC)
obj2<- seppop(microsats)

MicrsatFULL_RAD_POPSub <- repool(obj2$BF, obj2$CAKE , obj2$TU , obj2$POLEN , obj2$MOAT, obj2$SD, obj2$STYV, obj2$SK, obj2$CALK, obj2$OU, obj2$PRO, obj2$COP, obj2$OBY, obj2$PED, obj2$TROM, obj2$WEN)

pop(MicrsatFULL_RAD_POPSub)
MicrsatFULL_RAD_summary<- summary(MicrsatFULL_RAD_POPSub)
MicrsatFULL_RAD_summary
###### Summary stats fig ##########
par(mfrow=c(2,2))
allele_bar <- barplot(MicrsatFULL_RAD_summary$pop.nall,main="Number of Alleles per population",ylab="Number of genotypes",las=3, xaxt = 'n')
axis(1, levels(pop(MicrsatFULL_RAD_POPSub)), at = allele_bar, las = 2, cex.axis = 0.8)
abline(h = mean(MicrsatFULL_RAD_summary$pop.nall), col = "red", lty = 2)
loc_bar <- barplot(MicrsatFULL_RAD_summary$loc.nall,ylab="Number of alleles", main="Number of alleles per locus", xaxt = 'n')
axis(1, MicrsatFULL_RAD_POPSub$loc.names, at = loc_bar, cex.axis = 0.8, las = 2)
abline(h = mean(MicrsatFULL_RAD_summary$loc.nall), col = "red", lty = 2)
barplot(MicrsatFULL_RAD_summary$Hexp-MicrsatFULL_RAD_summary$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
abline(h = mean(MicrsatFULL_RAD_summary$Hexp-MicrsatFULL_RAD_summary$Hobs), col = "red", lty = 2)
barplot(MicrsatFULL_RAD_summary$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)
abline(h = mean(MicrsatFULL_RAD_summary$pop.eff), col = "red", lty = 2)

mean(MicrsatFULL_RAD_summary$pop.eff)
sd(MicrsatFULL_RAD_summary$pop.eff)
mean(MicrsatFULL_RAD_summary$loc.nall)
sd(MicrsatFULL_RAD_summary$loc.nall)
mean(MicrsatFULL_RAD_summary$pop.nall)
sd(MicrsatFULL_RAD_summary$pop.nall)

########### Comparative sumstats fig ###########

comp <- read.csv("~/Dropbox/PhD/Dan's PhD (Shared)/Data/Micro_RAD_comp/comparative_sum_stats_table.csv", header = T)
comp
library(ggplot2)
library(reshape)

comp_melted <- melt(comp)
comp_melted

ggplot(comp_melted, aes(Dataset, value, fill = variable))+geom_bar(position = "dodge", stat = "identity")

comp
comp_Micro <- comp[2:4,c(1,3:6)]
comp_Micro
comp_micro_melted <- melt(comp_Micro)
ggplot(comp_micro_melted, aes(Dataset, value, fill = variable))+geom_bar(position = "dodge", stat = "identity", color = "white")

##### stats ########

## Paired t-tests ##

# 1. Samples/pop 

M2_vs_M3_samples_per_pop <- t.test(RAD_only_summary$pop.eff, MicrsatFULL_RAD_summary$pop.eff, paired = T)
M2_vs_M3_samples_per_pop

# 2. alleles/pop

M2_vs_M3_alleles_per_pop <- t.test(RAD_only_summary$pop.nall, MicrsatFULL_RAD_summary$pop.nall, paired = T)
M2_vs_M3_alleles_per_pop

#3. alleles/loc

M2_vs_M3_alleles_per_loc <- t.test(RAD_only_summary$loc.nall, MicrsatFULL_RAD_summary$loc.nall, paired = T)
M2_vs_M3_alleles_per_loc


####### looking at the differences between UK and Denmark across M2 and M3 ####################

## pooling UK and Danish samples ##


## M2 subset ##
names(obj) ## RAD indivs only
UK_DM_M2 <- repool(obj$PED,  obj$COP, obj$SK, obj$BF, obj$MOAT, obj$CAKE)
pop(UK_DM_M2)

## M3 subset ##
names(obj2) ## all indivs in RAD pops
UK_DM_M3 <- repool(obj2$PED,  obj2$COP, obj2$SK, obj2$BF, obj2$MOAT, obj2$CAKE)
pop(UK_DM_M3)

sum(UK_DM_M2$loc.nall)
sum(UK_DM_M3$loc.nall)

UK_DM_M2_summ <- summary(UK_DM_M2)
UK_DM_M3_summ <- summary(UK_DM_M3)

UK_DM_M2_summ$pop.nall
UK_DM_M3_summ$pop.nall






######### Allele frequency differences between populations #############

#M1
microsats_genpop <- genind2genpop(microsats)
microsats_freqs <- makefreq(microsats_genpop, missing = NA, truenames = TRUE)
microsats_freqs

#M2 
RAD_only_microsats_subnoSTEC
RAD_only_microsats_subnoSTEC_genpop <- genind2genpop(RAD_only_microsats_subnoSTEC)
RAD_only_microsats_subnoSTEC_freqs <- makefreq(RAD_only_microsats_subnoSTEC_genpop, missing = NA, truenames = T)
write.csv(RAD_only_microsats_subnoSTEC_freqs$tab, "~/Dropbox/PhD/Dan's PhD (Shared)/Data/Micro_RAD_comp/M2_allele_freqs.csv")


#M3
MicrsatFULL_RAD_POPSub_genpop <- genind2genpop(MicrsatFULL_RAD_POPSub)
MicrsatFULL_RAD_POPSub_freqs <- makefreq(MicrsatFULL_RAD_POPSub_genpop, missing = NA, truenames = TRUE)
write.csv(MicrsatFULL_RAD_POPSub_freqs$tab, "~/Dropbox/PhD/Dan's PhD (Shared)/Data/Micro_RAD_comp/M3_allele_freqs.csv")


alleles(microsats)
alleles(RAD_only_microsats_subnoSTEC)$L08
alleles(MicrsatFULL_RAD_POPSub)$L08

