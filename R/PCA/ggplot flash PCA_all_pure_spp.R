
setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_Pure_spp_only/treemix_unpooled_cru/populations_p_24_r07_blacklisted/FlashPCA")



pcs <- read.delim("PCS_names_formatted.csv", header = F, sep = ",")
names(pcs) <- c("samples", "pop","pools", "PC1",  "PC2" , "PC3" , "PC4" , "PC5" , "PC6" , "PC7" , "PC8" , "PC9" , "PC10")

head(pcs)

head(pcs)
library(ggplot2)
library(directlabel)

pca <- qplot(PC1,PC2, data = pcs, color = pop, cex = 2)
pca + annotate("text", x = pcs$PC1, y = pcs$PC2, label = pcs$samples, cex = 4)

pcs[,c(1:4)]
crucian <- pcs[c(1:43,49:54,63:69,73:82,88:134,137:181,183:207,209:217),]

crucian
cru_pca <- qplot(PC1,PC2, data = crucian, color = pop, cex = 4)
cru_pca + annotate("text", x = crucian$PC1, y = crucian$PC2, label = crucian$samples, cex = 4)



pca <- qplot(PC1,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC1, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate



pca <- qplot(PC2,PC3, data = pcs, color = pop, cex = 2)

pca + annotate("text", x = pcs$PC2, y = pcs$PC3, label = pcs$samples, cex = 4)
?annotate
