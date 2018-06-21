setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples/populations_r07_p27_m8")

pcs <- read.delim("pcs_names_and_pops_formatted.txt", header = F, sep = ",")
names(pcs) <- c("samples", "pop", "PC1",  "PC2" , "PC3" , "PC4" , "PC5" , "PC6" , "PC7" , "PC8" , "PC9" , "PC10")

head(pcs)
library(ggplot2)

pca <- qplot(PC1,PC2, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC1, y = pcs$PC2, label = pcs$samples, cex = 4)
?annotate

pca <- qplot(PC1,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC1, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate



pca <- qplot(PC2,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC2, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate
