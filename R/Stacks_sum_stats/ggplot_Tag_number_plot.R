library(ggplot2) 

######################################### BATCH 1 ###############################################

## Cru_gold_hyb_1 ###

setwd('/media//dan//34D5D1CE642D7E36//2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Incremental_tests/Batch_1/cru_gold_hyb_1/')
#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()

### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()


### Cru_gold_hyb_2 ###

setwd('../cru_gold_hyb_2/')
#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()


### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()


### Cru_gib_hyb_1 ###

setwd('../cru_gib_hyb_1/')

#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()


### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()



### Cru_gib_hyb_2 ###

setwd('../cru_gib_hyb_2/')

#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()


### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()


################################### BATCH 2 ###########################################
### all_cru ###

setwd('../../Batch_2/all_cru/')

#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()


### cru_gold_putHyb ##########

setwd('../cru_gold_putHyb/')

### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()

#### Change in tag plot ##
CH <- read.delim("all_tag_changes.txt", header = T) ## make sure data is in the right format
p <- qplot(data = CH,m_value,No_Tags, colour= Sample, 
           geom = c("point", "line"), 
           scale_x_continuous(breaks = seq(0,14,by=2)), 
           main = "Change in Tag Number with read-depth(m) incrementation") + scale_x_continuous(breaks=seq(0,14,by=2))
p
## save to pdf ##
pdf("./ggplot_m_tests_change_in_Tags_all.pdf")
p
dev.off()


### Change in coverage plot ###

allcov <- read.delim("all_m_test_coverage_data.txt", header = F)
names(allcov) = c("Sample", "m_value", "Tag_coverage")
allcov


q <- qplot(data = allcov,m_value, Tag_coverage, colour= Sample, geom = c("point", "line"), scale_x_continuous(breaks = seq(0,14,by=2)), main = "Change in coverage with read-depth(m) incrementation")
q

pdf("./ggplot_m_test_change_in_cov_all.pdf")
q
dev.off()


## MS tests tag number plots










