setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/PopPool_analysis_r0.5_p_all4/")

data <- read.delim("single_snp_missing_data_per_sample.txt", header = T)

head(data)

library(ggplot2)

bar <- qplot(sample, data = data, weight = X._missing_loci, geom = "bar")

par(mfrow = c(1,1), mar = c(10,2,2,2))
mybar <- barplot(data$X._missing_loci, space = 3, main = "Missing data across the reference guided stacks analysis")
axis(1, labels = data$sample, at = mybar, cex.axis = 0.5, las = 2)


######### Hits per scaffold ##########

scaffs <- read.delim("Number_SNPs_per_scaffold.txt", header = F, sep = "\t")
head(scaffs)

scaffBar <- barplot(scaffs$V1, space = 3, main = "Number of SNPs per scaffold")
axis(1, labels = scaffs$V2, at = scaffBar, cex.axis = 0.1, las = 2)



###### FLash PCA ###########


pcs <- read.delim("PCS_names_pops_form.txt", header = F, sep = ",")
names(pcs) <- c("pop", "samples", "PC1",  "PC2" , "PC3" , "PC4" , "PC5" , "PC6" , "PC7" , "PC8" , "PC9" , "PC10")

head(pcs)
library(ggplot2)

####### First component #######

## PC1 contains crucian -> goldfish -> gibel variation
## PC2 contains crucian -> common carp variation.

pca <- qplot(PC1,PC2, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC1, y = pcs$PC2, label = pcs$samples, cex = 4)
?annotate

pca <- qplot(PC1,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC1, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate



pca <- qplot(PC2,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC2, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate
