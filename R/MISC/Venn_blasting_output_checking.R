setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Blast_tag_sharing_analysis/fully_shared_tags")
par(mfrow = c(1,2))
library(calibrate)

####### CRUCIAN vs GIBEL #########

crucian_gibel_blast <- read.delim("0_1_blast_results_tab.txt", header = F)

head(crucian_gibel_blast)
names(crucian_gibel_blast)  <- c("seqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
summary(crucian_gibel_blast$pident)

hist(crucian_gibel_blast$pident, breaks = 19, main = "% identity")
hist(crucian_gibel_blast$gapopen, main= "Number of gaps (indels)")

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Blast_tag_sharing_analysis")


cru_hyb_VS_gib <- read.delim("1_2_vs_0_blast_results_tab_max1.txt", header = F)
cru_hyb_VS_gib_hitcounts <- read.delim("1_2_vs_0_blast_results_N_Hits.txt", header = F, sep = " ")

names(cru_hyb_VS_gib)  <- c("seqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(cru_hyb_VS_gib)
head(cru_hyb_VS_gib_hitcounts)
summary(cru_hyb_VS_gib_hitcounts$V1)
summary(cru_hyb_VS_gib$gapopen)

par(mfrow = c(2,2), mar = c(4,4,1,1))

## Alignment identity histogram ###

hist(cru_hyb_VS_gib$pident, breaks = 19, main = "% identity")

## Number of INDELS ##

hist(cru_hyb_VS_gib$gapopen, main= "Number of gaps (indels)", breaks = 9)
table(cru_hyb_VS_gib$gapopen)
## Number of hits per query ###

myhist <- hist(cru_hyb_VS_gib_hitcounts$V1, main= "Number of hits per query", breaks = 140, ylim = c(0,1000), xaxt = 'n')
axis(1,  seq(0,140,5), at = seq(0,140,5),  las = 2, cex = 0.5)
textxy(5,800, "<----- 1 HIT", col = "red", cex = 1)

## length of each alignment ##

summary(cru_hyb_VS_gib$length)
hist(cru_hyb_VS_gib$length, breaks = 70, xlim = c(0,100), xaxt = 'n')
axis(1, labels = seq(0,100,5), at = seq(0,100,5), cex.axis = 0.7, las = 2)
textxy(65, 2500, "92bp ---->", col = "red", cex = 1)


##### Gib + Hyb v crucian #########



gib_hyb_VS_cru <- read.delim("0_2_Vs_1_blast_results_tab_max1.txt", header = F)

names(gib_hyb_VS_cru)  <- c("seqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(gib_hyb_VS_cru)

summary(gib_hyb_VS_cru$pident)
summary(gib_hyb_VS_cru$gapopen)

par(mfrow = c(2,2), mar = c(4,4,1,1))

## Alignment identity histogram ###

hist(gib_hyb_VS_cru$pident, breaks = 19, main = "% identity")

## Number of INDELS ##

hist(gib_hyb_VS_cru$gapopen, main= "Number of gaps (indels)", breaks =20)
table(gib_hyb_VS_cru$gapopen)

## length of each alignment ##

summary(gib_hyb_VS_cru$length)
hist(gib_hyb_VS_cru$length, breaks = 70, xlim = c(0,100), xaxt = 'n')
axis(1, labels = seq(0,100,5), at = seq(0,100,5), cex.axis = 0.7, las = 2)
textxy(65, 500, "92bp ---->", col = "red", cex = 1)
?lapply

