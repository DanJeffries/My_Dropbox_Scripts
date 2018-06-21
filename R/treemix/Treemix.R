setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_Pure_spp_only/treemix_unpooled_cru/")

source("~/RAD_programs/TreeMix/treemix-1.12/src/plotting_funcs.R")


## Looking at datasets ---------------------------------------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_Pure_spp_only/treemix_spp_level/populations_p_5_r07_blacklisted")

pooled_spp <- read.PLINK('plink.raw' , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(pooled_spp, posi="topleft", yaxt = 'n') ## Not too much missing data!
pooled_spp$n.loc
p24 <- read.PLINK('populations_no_SD_p24_r07/plink.raw' , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(p24, posi="topleft", yaxt = 'n') ## Not too much missing data!
p24$n.loc

spp_codes <- read.delim("Spp_codes.txt", header = F)

pop(p24) <- spp_codes$V1


seploc

table(pop(p24))

data.frame(spp_assignments)

pop_clust$V3
pop(p24)


p25 <- read.PLINK('populations_no_SD_p25_r07/plink.raw' , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(p25, posi="topleft", yaxt = 'n') ## Not too much missing data!
p25$n.loc

p26 <- read.PLINK('populations_no_SD_p26_r07/plink.raw' , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(p26, posi="topleft", yaxt = 'n') ## Not too much missing data!
p26$n.loc


plot_tree("populations_no_SD_p25_r07/Treemix_mig_0/out0")
fs <- get_f("populations_no_SD_p25_r07/Treemix_mig_0/out0")

plot_tree("populations_no_SD_p25_r07/Treemix_mig_1/out1")
fs <- append(fs, get_f("populations_no_SD_p25_r07/Treemix_mig_1/out1"))

plot_tree("populations_no_SD_p25_r07/Treemix_mig_2/out2")
fs <- append(fs, get_f("populations_no_SD_p25_r07/Treemix_mig_2/out2"))

plot_tree("populations_no_SD_p25_r07/Treemix_mig_3/out3")
fs <- append(fs, get_f("populations_no_SD_p25_r07/Treemix_mig_3/out3"))

plot_tree("populations_no_SD_p25_r07/Treemix_mig_4/out4")
fs <- append(fs, get_f("populations_no_SD_p25_r07/Treemix_mig_4/out4"))

plot_tree("populations_no_SD_p25_r07/Treemix_mig_5/out5")
fs <- append(fs, get_f("populations_no_SD_p25_r07/Treemix_mig_5/out5"))



plot(seq(4),fs[c(1:4)], pch = 16)

get_f = function(stem){
  d = paste(stem, ".cov.gz", sep = "")
  d2 = paste(stem, ".modelcov.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  d2 = read.table(gzfile(d2), as.is = T, comment.char = "", quote = "")
  d = d[order(names(d)), order(names(d))]
  d2 = d2[order(names(d2)), order(names(d2))]
  tmpcf = vector()
  tmpmcf = vector()
  for (j in 1:nrow(d)){
    for (k in (j+1):nrow(d)){
      tmpcf = append(tmpcf, d[j,k])
      tmpmcf = append(tmpmcf, d[j,k] - d2[j,k])
    }
  }
  tmpv = var(tmpmcf)/var(tmpcf)
  return(1-tmpv)
  
}
