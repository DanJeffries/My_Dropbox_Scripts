### ------------------BAYENV ---------------------###

####---- MASTER -------------------------------------------------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV/Sub-samples/MASTER_npops_17")

bf_1 <- read.delim("bf_environ.STNDRDZED_pop_temps_1.env", header = F, row.names = 1)
bf_1$V3 <- NULL

bf_2 <- read.delim("bf_environ.STNDRDZED_pop_temps_2.env", header = F, row.names = 1)
bf_2$V3 <- NULL

bf_3 <- read.delim("bf_environ.STNDRDZED_pop_temps_3.env", header = F, row.names = 1)
bf_3$V3 <- NULL

bf_4 <- read.delim("bf_environ.STNDRDZED_pop_temps_4.env", header = F, row.names = 1)
bf_4$V3 <- NULL

bf_5 <- read.delim("bf_environ.STNDRDZED_pop_temps_5.env", header = F, row.names = 1)
bf_5$V3 <- NULL



# ---- Spearmans correlations ------------------------------------------------------------------------------------


cor.test(bf_1$V2, bf_2$V2, method = "spearman") # ***
cor.test(bf_1$V2, bf_3$V2, method = "spearman") # ***
cor.test(bf_1$V2, bf_4$V2, method = "spearman") # ***
cor.test(bf_1$V2, bf_5$V2, method = "spearman") # ***
cor.test(bf_2$V2, bf_3$V2, method = "spearman") # ***
cor.test(bf_2$V2, bf_4$V2, method = "spearman") # ***
cor.test(bf_2$V2, bf_5$V2, method = "spearman") # ***
cor.test(bf_3$V2, bf_4$V2, method = "spearman") # ***
cor.test(bf_3$V2, bf_5$V2, method = "spearman") # ***
cor.test(bf_4$V2, bf_5$V2, method = "spearman") # ***

## -Plotting -----------------------------------------------------------------------------------------------------------------

#hist(log(bf_1$V2), breaks = length(bf_1$V2))
#cutoff <- quantile(bf_1$V2, 0.95)
#abline(lty = 3, v = cutoff, col = "red")


bf_1_top <- plot_Bayenv2(perc_cutoff = 0.99, id = T, bayenv_file = "bf_environ.STNDRDZED_pop_temps_1.env", map_file = "Loc_keys.txt", ylab = "BFs")
length(bf_1_top$BAYENV_ID)

bf_2_top <- plot_Bayenv2(perc_cutoff = 0.99, id = T, bayenv_file = "bf_environ.STNDRDZED_pop_temps_2.env", map_file = "Loc_keys.txt", ylab = "BFs")
length(bf_2_top$BAYENV_ID)

bf_3_top <- plot_Bayenv2(perc_cutoff = 0.99, id = T, bayenv_file = "bf_environ.STNDRDZED_pop_temps_3.env", map_file = "Loc_keys.txt", ylab = "BFs")
length(bf_3_top$BAYENV_ID)

bf_4_top <- plot_Bayenv2(perc_cutoff = 0.99, id = T, bayenv_file = "bf_environ.STNDRDZED_pop_temps_4.env", map_file = "Loc_keys.txt", ylab = "BFs")
length(bf_4_top$BAYENV_ID)

bf_5_top <- plot_Bayenv2(perc_cutoff = 0.99, id = T, bayenv_file = "bf_environ.STNDRDZED_pop_temps_5.env", map_file = "Loc_keys.txt", ylab = "BFs")
length(bf_5_top$BAYENV_ID)


## Combining top 1% of each run 

All_tops <- data.frame(run_1 = bf_1_top$Stacks_ID,run_2 = bf_2_top$Stacks_ID ,run_3 = bf_3_top$Stacks_ID ,run_4 = bf_4_top$Stacks_ID ,run_5 = bf_5_top$Stacks_ID)

All_tops

All_combined <- data.frame(all = ordered(c(All_tops$run_1, All_tops$run_2, All_tops$run_3, All_tops$run_4, All_tops$run_5)))

Top_loc_counts <- as.data.frame(table(All_combined$all))
names(Top_loc_counts) <- c("Stacks_ID", "Run_count")
Top_loc_counts_ordered <- Top_loc_counts[order(Top_loc_counts$Run_count),]
Top_loc_counts_ordered

length(Top_loc_counts_ordered$Stacks_ID[which(Top_loc_counts_ordered, Top_loc_counts_ordered$Run_count > 5)])

b <- barplot(Top_loc_counts_ordered$Run_count, main = "Loci with BFs in top 1% of each run", ylab = "# Runs", col = "darkblue")
axis(1, labels = Top_loc_counts_ordered$Stacks_ID, at = b, las = 2, cex.axis = 0.7)



## Plotting Function-----------------------------------------------------------------------------------------------------------


#### define function

plot_Bayenv2 <- function(bayenv_file,map_file,perc_cutoff=0,abs_cutoff=0,id=FALSE,xlab="",ylab="",log="") {
  data <- read.delim(bayenv_file, header=F)
  map <- read.delim(map_file, header=F)
  x_all = data
  x_all[,3] <- map$V1[1:length(x_all$V2)] ## change
  x_all[,4] <- seq(1,length(x_all$V2))
  #x_all
  
  cutoff = 0
  if (perc_cutoff > 0){
    cutoff <- quantile(x_all$V2, perc_cutoff, names=F)
  }
  if (abs_cutoff > 0){
    cutoff = abs_cutoff
  }
  
  sub <- x_all[ which(x_all$V2 > cutoff),]
  
  #  sub_main <- x_all[ which(x_all$V2 <= perc),]
  #sub_top_perc
  
  plot(x_all$V4,x_all$V2, cex=0.4, pch=19, col="grey", ylab=ylab, xlab=xlab, log=log)
  #points(sub_main$V4,sub_main$V2, cex=0.5, pch=19, col="grey")
  if (cutoff > 0){
    abline(cutoff,0,lty=3,untf = TRUE)
  }
  points(sub$V4,sub$V2, cex=0.5, pch=19)
  if (id == TRUE){
    text(sub$V4, sub$V2, labels = sub$V3, cex=0.6, pos=4)
  }
  names(sub) <- c("BAYENV_ID", "BF", "Stacks_ID", "x_coord")
  return(sub)
}


