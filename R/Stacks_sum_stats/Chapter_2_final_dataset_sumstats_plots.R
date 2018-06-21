library(ggplot2)
#install.packages("gridBase")
library(gridBase)
library(grid)
#install.packages("gridExtra")
library(gridExtra)

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/Final_data_files/")


#### Incremental USTACKS tests ############# ----------------------------------------------------
#-----------------------------------------------------------------------------------------------

setwd("./Incremental_U_tests/")

## M tests -----------------------------------------------------------------

#### Change in tag plot ##
all <- read.delim("Ustacks_M_mismatch_test_tag_numbers.tsv", header = T) ## make sure data is in the right format

tags_M <- ggplot(all) + geom_line(aes(value, tag_numb, group = sample, color =  ID), size = 0.5) +
                        geom_point(aes(value, tag_numb, group = sample, color =  ID), size = 2)

## Change in coverage plot ## 

all_cov_M <- read.delim("All_samples_Mtests_av_coverage_data.txt", header = F) ## make sure data is in the right format
names(all_cov_M) <- c("sample", "M", "coverage")
all_cov_M$ID <- all$ID  ## AS the order of samples is the same, I can just add these spp ids here
cov_M <- ggplot(all_cov_M) + geom_line(aes(M, coverage, group = sample, color =  ID), size = 0.5) +
                             geom_point(aes(M, coverage, group = sample, color =  ID), size = 2)


grid.arrange(tags_M, cov_M, ncol=2, main="Incremental M tests", clip = T)


## m tests -----------------------------------------------------------------

#### Change in tag plot ##
all_tags_m <- read.delim("Ustacks_m_cov_test_tag_numbers.tsv", header = T) ## make sure data is in the right format

tags_m <- ggplot(all_tags_m) + geom_line(aes(value, tag_numb, group = sample, color =  ID), size = 0.5) +
                  geom_point(aes(value, tag_numb, group = sample, color =  ID), size = 2)


## Change in coverage plot ## 

all_cov_m <- read.delim("All_samples_m_cov_tests_av_coverage_data.txt", header = F) ## make sure data is in the right format
names(all_cov_m) <- c("sample", "m", "coverage")
all_cov_m$ID <- all$ID  ## AS the order of samples is the same, I can just add these spp ids here
cov_m <- ggplot(all_cov_m) + geom_line(aes(m, coverage, group = sample, color =  ID), size = 0.5) +
                             geom_point(aes(m, coverage, group = sample, color =  ID), size = 2)

grid.arrange(tags_m, cov_m, ncol=2, main="Incremental m tests", clip = T)




# INCREMENTAL_C ----------------------------------------------------------------

## N tests

## Tag_share plot horizontal bargraphs

par(mfrow = c(1,1))
tags <- read.delim("./Incremental_C_example_tag_share_numbers.txt", row.names = 1)
ttags <- as.table(t(tags))

barplot(ttags, beside = T, horiz = T)





## WHOLE DATASETS - POPULATIONS TESTS--------------------------------------------------------
#--------------------------------------------------------------------------------------------

### Treemix all pure species dataset ###


Treemix_data <- read.PLINK("Treemix_data/plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
names <- read.delim("Treemix_data/indiv_names.txt", header = F)

myglPlot(Treemix_data, posi="topleft", yaxt = 'n', axes = F, XaxisCEX = 0.1, YaxisCEX = 0.3, ylabels = names$V1) ## Not too much missing data!

cov_data <- read.delim("Treemix_data/Average_coverage.txt", header = F)

hist(cov_data$V1, xlim = c(0, 100), breaks = max(cov_data$V1))



## Populations tests ## ---------------------------------------------------------

setwd("Populations_tests/")


pop_tests <- read.delim("treemix_data_populations_tests_snp_numbers.txt")

p_tests <- pop_tests[c(1,3:5,8, 20),]
r_tests <- pop_tests[c(2,6:8,19),]
r_tests <- r_tests[order(r_tests$r),]
m_tests <- pop_tests[c(8:18),]
m_tests <- m_tests[order(m_tests$m),]


mat <- matrix(c(1,4,2,4,3,4), 2,3) ## Nifty way of making custom plot layout
layout(mat)

## -p tests

plot((26-p_tests$p), p_tests$tags, type = "c", main = "N SNPs with change in p", xlab = "Minimum N pops per tag", ylab = "N tags")
points((26-p_tests$p), p_tests$tags)

## -r tests

plot(r_tests$r, r_tests$tags, type = "c", main = "N SNPs with change in r", xlab = "Minimum N samples/Pop/tag", ylab = "N tags")
points(r_tests$r, r_tests$tags)

## -m tests
plot(m_tests$m, m_tests$tags, type = "c", main = "N SNPs with change in m", xlab = "Minimum read depth per tag", ylab = "N tags")
points(m_tests$m, m_tests$tags)



## Tag dropout plots

tag_dropouts <- read.delim("Populations_tests/treemix_data_pop_tests_dropouts.txt")

p_tests <- tag_dropouts[which(tag_dropouts$r == "r07" & tag_dropouts$m == "0"),]
p_tests_ordered <- p_tests[order(p_tests$p),]

p_tests_plot <- table(data.frame(p_tests_ordered$p, p_tests_ordered$pop, p_tests_ordered$dropouts))


par(mfrow = c(2,2))

## This bit allows me to add a ggplot to a panel that already has r plots on it (dont know how though!)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(.5,.5,.5,.5)) 

cols <- grey.colors(6)

p <- ggplot(p_tests_ordered, aes(x=factor(pop), y=dropouts)) +
  geom_bar(position="dodge", stat="identity",aes(fill=factor(p)))+
  scale_fill_manual("p", values = cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, colour = "black"))+
  theme(axis.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.line = element_line(colour = 'grey'))+
  ggtitle("SNP dropout across populations") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))

print(p,vp = vp1)


## All samples denovo dataset ## ------------------------------------------------------------------------

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/Final_data_files/Populations_tests/")

#svg(file = "~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/MS Figures/Test_dev.svg")

pop_tests1 <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/Final_data_files/Populations_tests/denovo_all_samples_populations_tests_snp_numbers.txt")
pop_tests <- pop_tests1[c(1,2,4:21),] ## get rid of duplocated test

p_tests <- pop_tests[which(pop_tests$r == 0.7 & pop_tests$m == 0),]
p_tests <- p_tests[order(p_tests$p),]
r_tests <- pop_tests[which(pop_tests$p == 27 & pop_tests$m == 0),]
r_tests <- r_tests[order(r_tests$r),]
m_tests <- pop_tests[which(pop_tests$p == 27 & pop_tests$r == 0.7),]
m_tests <- m_tests[order(m_tests$m),]


mat <- matrix(c(1,4,7,2,5,7,3,6,7), 3,3) ## Nifty way of making custom plot layout
layout(mat)

## -p tests

plot(p_tests$p, p_tests$SNPs, type = "c", main = "-p", xlab = "", ylab = "N tags")
points(p_tests$p, p_tests$SNPs)

## -r tests

plot(r_tests$r, r_tests$SNPs, type = "c", main = "-r", xlab = "", ylab = "")
points(r_tests$r, r_tests$SNPs)

## -m tests
plot(m_tests$m, m_tests$SNPs, type = "c", main = "-m", xlab = "", ylab = "")
points(m_tests$m, m_tests$SNPs)



## Avererage sample coverage (from VCF) for changing -p ###---------------------------------------------------------------------


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M4_m4_N4/populations_tests/")

#p31
av_covs <- read.delim("populations__p31_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)

#p30
av_covs <- read.delim("populations__p30_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p29
av_covs <- read.delim("populations__p29_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p28
av_covs <- read.delim("populations__p28_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p27
av_covs <- read.delim("populations__p27_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p26
av_covs <- read.delim("populations__p26_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p25
av_covs <- read.delim("populations__p25_r07/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


## data for change in average coverage for -p

pure_means_df <- data.frame(p=rev(seq(25,31)), pure_means)
hyb_means_df <- data.frame(p=rev(seq(25,31)), hyb_means)

plot(hyb_means_df$p, hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", ylab = "Average coverage per sample", lwd = 1.5)
lines(pure_means_df$p, pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)



## coverage for changing -r -------------------------------------------------------

#r .5
av_covs <- read.delim("populations_p27_r00.5/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)

# r .6
av_covs <- read.delim("populations_p27_r00.6/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

# r .7
av_covs <- read.delim("populations_p27_r00.7/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

# r .8
av_covs <- read.delim("populations_p27_r00.8/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


# r .9
av_covs <- read.delim("populations_p27_r00.9/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


# r 1
av_covs <- read.delim("populations_p27_r01.0/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


pure_means_df <- data.frame(p=rev(seq(0.5, 1, 0.1)), pure_means)
hyb_means_df <- data.frame(p=rev(seq(0.5, 1, 0.1)), hyb_means)

plot(rev(hyb_means_df$p), hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", lwd = 1.5, ylim = c(20,32), ylab = "")
lines(rev(pure_means_df$p), pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)


## coverage for changing -m --------------------------------------------------

#m1
av_covs <- read.delim("populations_p27_r07_m1/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)

#m2
av_covs <- read.delim("populations_p27_r07_m2/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#m3
av_covs <- read.delim("populations_p27_r07_m3/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


#m4
av_covs <- read.delim("populations_p27_r07_m4/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


#m5
av_covs <- read.delim("populations_p27_r07_m5/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


#m6
av_covs <- read.delim("populations_p27_r07_m6/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))



#m7
av_covs <- read.delim("populations_p27_r07_m7/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


#m8
av_covs <- read.delim("populations_p27_r07_m7/sample_average_coverages_ammended.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


pure_means_df <- data.frame(p=rev(seq(1,8)), pure_means)
hyb_means_df <- data.frame(p=rev(seq(1,8)), hyb_means)

plot(rev(hyb_means_df$p), hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", lwd = 1.5, ylim = c(20,32), ylab = "")
lines(rev(pure_means_df$p), pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)



## Tag dropout plots

tag_dropouts <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/Final_data_files/Populations_tests/All_samples_denovo_populations_tests_snp_dropout_code_replacde.txt")

p_tests <- tag_dropouts[which(tag_dropouts$r == "r07" & tag_dropouts$m == "0"),]
p_tests_ordered <- p_tests[order(p_tests$p),]


## This bit allows me to add a ggplot to a panel that already has r plots on it (dont know how though!)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(.5,.5,.5,.5)) 

cols <- grey.colors(7)

p <- ggplot(p_tests_ordered, aes(x=factor(pop), y=dropouts)) +
  geom_bar(position="dodge", stat="identity",aes(fill=factor(p)))+
  scale_fill_manual("p", values = cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, colour = "black"))+
  theme(axis.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.line = element_line(colour = 'grey'))+
  ggtitle("SNP dropout across populations") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))

print(p,vp = vp1)


## All_samples_ref_aligned dataset ## ------------------------------------------------------------------------

pop_tests1 <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/populations_tests_snp_numbers.txt")


pop_tests <- pop_tests1[c(1,2,4:21),]

p_tests <- pop_tests[which(pop_tests$r == 0.7 & pop_tests$m == 0),]
p_tests <- p_tests[order(p_tests$p),]
r_tests <- pop_tests[which(pop_tests$p == 27 & pop_tests$m == 0),]
r_tests <- r_tests[order(r_tests$r),]
m_tests <- pop_tests[which(pop_tests$p == 27 & pop_tests$r == 0.7),]
m_tests <- m_tests[order(m_tests$m),]


mat <- matrix(c(1,4,7,2,5,7,3,6,7), 3,3) ## Nifty way of making custom plot layout
layout(mat)

## -p tests

plot(p_tests$p,p_tests$tags,  type = "c", main = "-p", xlab = "", ylab = "N tags")
points(p_tests$p,p_tests$tags)

## -r tests

plot(r_tests$r, r_tests$tags, type = "c", main = "-r", xlab = "", ylab = "")
points(r_tests$r, r_tests$tags)

## -m tests
plot(m_tests$m, m_tests$tags, type = "c", main = "-m", xlab = "", ylab = "")
points(m_tests$m, m_tests$tags)


## coverage plots ------------------

## Avererage sample coverage (from VCF) for changing -p ###---------------------------------------------------------------------


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/")

#p31
av_covs <- read.delim("populations__p31_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)

#p30
av_covs <- read.delim("populations__p30_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p29
av_covs <- read.delim("populations__p29_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p28
av_covs <- read.delim("populations__p28_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p27
av_covs <- read.delim("populations__p27_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p26
av_covs <- read.delim("populations__p26_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

#p25
av_covs <- read.delim("populations__p25_r07/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]

row.names(av_covs_sorted) <- seq(length(av_covs_sorted$V1))

av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:238),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238),]
pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


## data for change in average coverage for -p

pure_means_df <- data.frame(p=rev(seq(25,31)), pure_means)
hyb_means_df <- data.frame(p=rev(seq(25,31)), hyb_means)

plot(hyb_means_df$p, hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", ylab = "Average coverage per sample", lwd = 1.5)
lines(pure_means_df$p, pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)
## change in cov for -r --------------------------------

setwd("/media/dan//34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/")

av_covs <- read.delim("populations_p27_r00.5/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)

av_covs <- read.delim("populations_p27_r00.6/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r00.7/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r00.8/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r00.9/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r01.0//sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,233:237),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

pure_means_df <- data.frame(r=rev(seq(0.5, 1, 0.1)), pure_means)
hyb_means_df <- data.frame(r=rev(seq(0.5, 1, 0.1)), hyb_means)

plot(rev(hyb_means_df$r), hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", lwd = 1.5, ylim = c(27,31), ylab = "")
lines(rev(pure_means_df$r), pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)


## -- Change coverage change for m ------------------------------------------

av_covs <- read.delim("populations_p27_r07_m1/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

hyb_means <- mean(av_covs_hybs$V2)
pure_means <- mean(av_covs_pure$V2)


av_covs <- read.delim("populations_p27_r07_m2/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r07_m3/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r07_m4/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r07_m5/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

av_covs <- read.delim("populations_p27_r07_m6/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


av_covs <- read.delim("populations_p27_r07_m7/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))


av_covs <- read.delim("populations_p27_r07_m8/sample_average_coverages.txt", header = F)
av_covs_sorted <- av_covs[order(av_covs$V1),]
av_covs_hybs <- av_covs_sorted[c(49,51,75:89,232:236),]
av_covs_pure <- av_covs_sorted[c(-49,-51,-75:-89,-233:-238, -204:-207),]

pure_means <- append(pure_means, mean(av_covs_pure$V2))
hyb_means <- append(hyb_means, mean(av_covs_hybs$V2))

pure_means_df <- data.frame(m=rev(seq(1,8)), pure_means)
hyb_means_df <- data.frame(m=rev(seq(1,8)), hyb_means)

plot(rev(hyb_means_df$m), hyb_means_df$hyb_means, type = "b", col = "firebrick1", xlab = "", lwd = 1.5, ylim = c(27,33), ylab = "")
lines(rev(pure_means_df$m), pure_means_df$pure_means, type = "b", col = "blue", lwd = 1.5)


## Tag dropout plot

tag_dropouts <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_m_4/populations_tests/All_samples_ref_aligned_snp_dropout_code_replacde.txt")

p_tests <- tag_dropouts[which(tag_dropouts$r == "r07" & tag_dropouts$m == "0"),]
p_tests_ordered <- p_tests[order(p_tests$p),]



## This bit allows me to add a ggplot to a panel that already has r plots on it (dont know how though!)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(.5,.5,.5,.5)) 

cols <- grey.colors(7)

p <- ggplot(p_tests_ordered, aes(x=factor(pop), y=dropouts)) +
  geom_bar(position="dodge", stat="identity",aes(fill=factor(p)))+
  scale_fill_manual("p", values = cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, colour = "black"))+
  theme(axis.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.line = element_line(colour = 'grey'))+
  ggtitle("SNP dropout across populations") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))

print(p,vp = vp1)


## Pure_cru denovo dataset ## ------------------------------------------------------------------------


pop_tests1 <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_tests/populations_tests_snp_numbers.txt")
pop_tests <- pop_tests1[c(1,2,4:21),]

p_tests <- pop_tests[which(pop_tests$r == 0.7 & pop_tests$m == 0),]
p_tests <- p_tests[order(p_tests$p),]
r_tests <- pop_tests1[which(pop_tests$p == 27 & pop_tests$m == 0),]
r_tests <- r_tests[order(r_tests$r),]
m_tests <- pop_tests[which(pop_tests$p == 27 & pop_tests$r == 0.7),]
m_tests <- m_tests[order(m_tests$m),]


mat <- matrix(c(1,4,7,2,5,7,3,6,7), 3,3) ## Nifty way of making custom plot layout
layout(mat)

## -p tests

plot((19-p_tests$p), p_tests$snps, type = "c", main = "-p", xlab = "", ylab = "N tags")
points((19-p_tests$p), p_tests$snps)

## -r tests

plot(r_tests$r, r_tests$snps, type = "c", main = "-r", xlab = "", ylab = "")
points(r_tests$r, r_tests$snps)

## -m tests
plot(m_tests$m, m_tests$snps, type = "c", main = "-m", xlab = "", ylab = "")
points(m_tests$m, m_tests$snps)


# -p coverage plots -------------------------------------------
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_tests/")

av_covs <- read.delim("populations__p13_r07/sample_average_coverages.txt", header = F)
cov_means <- mean(av_covs$V2)
av_covs <- read.delim("populations__p14_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations__p15_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations__p16_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations__p17_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations__p18_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations__p19_r07/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))

means_df <- data.frame(r=rev(seq(13, 19)), cov_means)

plot(rev(means_df$r), means_df$cov_means, type = "b", col = "blue", xlab = "", lwd = 1.5, ylab = "")

# -r coverage plots ------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_tests/")

av_covs <- read.delim("populations_p27_r00.5/sample_average_coverages.txt", header = F)
cov_means <- mean(av_covs$V2)
av_covs <- read.delim("populations_p27_r00.6/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r00.7/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r00.8/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r00.9/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r01.0/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))

means_df <- data.frame(r=rev(seq(0.5, 1, 0.1)), cov_means)

plot(rev(means_df$r), means_df$cov_means, type = "b", col = "blue", xlab = "", lwd = 1.5, ylab = "")



## -- Change coverage change for m ------------------------------------------

av_covs <- read.delim("populations_p27_r07_m1/sample_average_coverages.txt", header = F)
cov_means <- mean(av_covs$V2)
av_covs <- read.delim("populations_p27_r07_m2/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m3/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m4/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m5/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m6/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m7/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))
av_covs <- read.delim("populations_p27_r07_m8/sample_average_coverages.txt", header = F)
cov_means <- append(cov_means  , mean(av_covs$V2))

means_df <- data.frame(r=rev(seq(1,8)), cov_means)

plot(rev(means_df$r), means_df$cov_means, type = "b", col = "firebrick1", xlab = "", lwd = 1.5, ylab = "")
lines(rev(means_df$r), means_df$cov_means, type = "b", col = "blue", lwd = 1.5)



## Tag dropout plots

tag_dropouts <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_2_RAD_methods_chapter/Final_data_files/Populations_tests/Pure_CRU_populations_tests_all_data_dropout_code_replacde.txt")

p_tests <- tag_dropouts[which(tag_dropouts$r == "r07" & tag_dropouts$m == "0"),]
p_tests_ordered <- p_tests[order(p_tests$p),]



## This bit allows me to add a ggplot to a panel that already has r plots on it (dont know how though!)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(.5,.5,.5,.5)) 

cols <- grey.colors(7)

p <- ggplot(p_tests_ordered, aes(x=factor(pop), y=dropouts)) +
  geom_bar(position="dodge", stat="identity",aes(fill=factor(p)))+
  scale_fill_manual("p", values = cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, colour = "black"))+
  theme(axis.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.line = element_line(colour = 'grey'))+
  ggtitle("SNP dropout across populations") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))

print(p,vp = vp1)



## Function for plotting the glPlot from adegenet, modified to add axis labels

myglPlot = function (x, col = NULL, legend = TRUE, posi = "bottomleft", 
                     bg = rgb(1, 1, 1, 0.5), XaxisCEX = 1, YaxisCEX = 1,  ylabels = NULL, xlabels = NULL, ...) 
{
  X <- t(as.matrix(x))
  X <- X[, ncol(X):1]
  ylabpos <- pretty(1:nInd(x), 5)
  if (is.null(col)) {
    myCol <- colorRampPalette(c("royalblue3", "firebrick1"))(max(X, 
                                                                 na.rm = TRUE) + 1)
  }
  else {
    myCol <- col
  }
  image(x = 1:nLoc(x), y = 1:nInd(x), z = X, xlab = "SNP index", ylab = "",
        yaxt = "n", col = myCol, ...)
  axis(side = 2, at = seq(nInd(x)), labels = rev(ylabels), cex.axis = YaxisCEX, las = 2)
  axis(side = 1, at = seq(nLoc(x)), labels = xlabels, cex.axis = XaxisCEX, las = 2)
  if (legend) {
    legend(posi, fill = myCol, legend = 0:max(X, na.rm = TRUE), 
           horiz = TRUE, bg = bg, title = "Number of 2nd allele")
  }
  return(invisible())
}

