## Looking at the SNP representation in the data

################# batch 1 only ##################
## checking the difference

## However, I don't think I used any r or p constraints on this populations run!

setwd('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_v.1/Cru_only')

pdf('./SNPrepresentation_pure_cru_batch1only.pdf')
hist(SNPcov$INFO, breaks = 76, main = "SNP representation in sample batch 1 pure crucian", xlab = "#Individuals", ylab = '#Number of Loci')
dev.off()

####### ALL DATA (BATCH1 and BATCH2)  ###############

setwd('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/populations/crucian_only_M8_N8/Pure_cruc_N_EU_r_0.8_P15')

SNPcov <- read.delim('SNP_coverage_only.txt', header = T) ## file prepared with a quick grep, cut bash one-liner and sed, from VCF see below
## this vcf was from a populations run using -r 0.8 and -p 17. No m threshold

# grep -v '^##' batch_1.vcf | cut -f3,8 |  cut -d';' -f1 | > SNP_coverage.txt                                                                                                                              [ 2:02PM]
# sed 's/NS=//g' SNP_coverage.txt > SNP_coverage_only.txt 

head(SNPcov) ## check
summary(SNPcov)

## Make histogram

pdf('./SNPrepresentation_pure_cru.pdf') 
hist(SNPcov$INFO, breaks = 76, main = "SNP representation in pure crucian", xlab = "#Individuals", ylab = '#Number of Loci')
dev.off()
### There are three worringly distinct peaks in this histogram - meaning that locus dropout is certainly not random! 
## So there are some populations that are dropping out in a lot of loci.


####### Looking at the distribution of the missing data accross my populations ###### 

## This data comes from my little 'dist of missing data in VCF.py' program ##

misSNPs <- read.delim('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_r0.8_p_17_m8/missing_loci_dist.txt', header = T)
misSNPsnames <- read.delim('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/populations/crucian_only_M8_N8/r_0.7_P_17_no_m_thresh/batch_1_missing_data_samplenames.txt', header = T)
summary(misSNPsnames)
misSNPs
bp <- barplot(misSNPs$X._missing_loci, axes = F)
axis(2)
axis(1,at = bp, labels = misSNPs$sample, las = 2, cex.axis = 0.3)
