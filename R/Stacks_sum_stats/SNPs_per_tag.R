## made input file from the batch_X.vcf file using the below bash oneliners:
# dan@BioPC4[NEU_r_0.7_p_16_m8_all_SNPs] grep -v '^#' batch_1.vcf |cut -f3 | uniq -c > SNPcount_per_tag.txt ## cut the column we want
# sed 's/  //g' SNPcount_per_tag.txt > SNPcount_per_tag_sripped.txt ## strip multiple spaces so R can read in nicely

setwd('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_r0.7_p_16_m8/')

body <- read.delim('SNPcount_per_tag.txt', sep = ' ', header = F)

names(body) = c('count', 'Tag_ID')

head(body)

summary(body$count)

library(calibrate)

pdf('./SNP_number_per_tag')
hist(body$count, breaks = 21, main = 'Number of SNPs per Tag', xlab = 'Number of SNPs', ylab = 'Tag count', xaxt = 'n')
axis(1, at = seq(1,21,1))

textxy(6,5000,"Dataset contains only first SNP of each locus. Whole dataset contained 26142 SNPs
- 13344 SNPs produced using the --write-single-SNP flag
- Mean number of SNPs per tag = 1.96.", pos = 4, cex = 1)

dev.off()



