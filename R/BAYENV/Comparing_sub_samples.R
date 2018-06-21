
## Comparing BAYENV results

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV/Sub-samples")

Bayenv_results <- read.delim("Numb_shared_loc_per_Subsample.tsv", sep = '\t', header = F)

names(Bayenv_results) = c("Numb_loci", "Numb_bottled_pops", "subsample")
length(Bayenv_results$V1)

bottled_stats <- read.delim("Sub_sample_stats/Mean_bottleneck_stats_per_subsample.tsv")
bottled_stats

bottled_stats_winners <- bottled_stats[c(4,7,8,9,10,12,14,17,21,22,27,28,30,32,38,39), ]
GW_winners <-  bottled_stats[c(4,7,8,9,10,12,14,17,21,22,27,28,30,32,38,39), 5]
GW_winners
Bayenv_results$G_W <- GW_winners
Bayenv_results

plot(Bayenv_results$G_W, Bayenv_results$Numb_loci, xlab = "Garza-Williamson index", ylab = "# loci shared with master", pch = 16)
text(Bayenv_results$G_W,Bayenv_results$Numb_loci, gsub("_bottlenecked_pops/rep_", "/",paste(Bayenv_results$Numb_bottled_pops, Bayenv_results$subsample, sep = "/")), pos = 2)




subs_temp_ranges <- read.delim("Sub_sample_stats/Sub_samp_temp_ranges.txt", header = F)
subs_temp_ranges
axis(1, )

winning_subs_temp_ranges <- subs_temp_ranges[c(4,7,8,9,10,12,14,17,21,22,27,28,30,32,38,39), ]
winning_subs_temp_ranges
Bayenv_results$Temp_range <- winning_subs_temp_ranges$V2

b <- barplot(subs_temp_ranges$V2)
axis(1, gsub("_bottlenecked_pops/rep_", "/", subs_temp_ranges$V1) , at = b, las = 3, cex.axis = 1)

barplot(Bayenv_results$Temp_range)


plot(Bayenv_results$Temp_range, Bayenv_results$Numb_loci, pch = 16)
text(Bayenv_results$Temp_range,Bayenv_results$Numb_loci, gsub("_bottlenecked_pops/rep_", "/",paste(Bayenv_results$Numb_bottled_pops, Bayenv_results$subsample, sep = "/")), pos = 2)
length(bottled_stats_winners$X.Pop)

plot(Bayenv_results$G_W, Bayenv_results$Temp_range, ylab = "Temp range", xlab = "GW")
text(Bayenv_results$G_W, Bayenv_results$Temp_range, gsub("_bottlenecked_pops/rep_", "/",paste(Bayenv_results$Numb_bottled_pops, Bayenv_results$subsample, sep = "/")), pos = 2)
Bayenv_results
cor.test(Bayenv_results$G_W[c(1:3,6:9,11:12,14:16)], Bayenv_results$Temp_range[c(1:3,6:9,11:12,14:16)], method = "spearman" )

plot(Bayenv_results$G_W[c(1:3,6:9,11:12,14:16)], Bayenv_results$Temp_range[c(1:3,6:9,11:12,14:16)], ylab = "Temp range", xlab = "GW")
text(Bayenv_results$G_W[c(1:3,6:9,11:12,14:16)], Bayenv_results$Temp_range[c(1:3,6:9,11:12,14:16)], gsub("_bottlenecked_pops/rep_", "/",paste(Bayenv_results$Numb_bottled_pops[c(1:3,6:9,11:13,14:16)], Bayenv_results$subsample[c(1:3,6:9,11:13,14:16)], sep = "/")), pos = 2)

# --- With zeroes ------------------------------------------------------------------------------------------------------------

counts <- read.delim("~/Desktop/BAYENV_counts", header = F)

plot(bottled_stats$GW_M_Ratio[c(3:55)], counts$V1, ylab = "Nloc", xlab = "GW", pch = 16)
text(bottled_stats$GW_M_Ratio[c(3:55)],counts$V1, gsub("_bottlenecked_pops/rep_", "/",paste(Bayenv_results$Numb_bottled_pops, Bayenv_results$subsample, sep = "/")), pos = 2)

plot(Bayenv_results$G_W, Bayenv_results$Temp_range, ylab = "Temp range", xlab = "GW")


## Locus counts ----------------------------------------------------------------------------------------------------------------

loc_counts <- read.delim("Shared_locus_counts.tsv", header = T, sep = " ")

plot(loc_counts$X1, pch = 16, col = "blue", xlab = "SNP", ylab = "Number of sub-samples")
text(seq(35), loc_counts$X1, loc_counts$Shared, pos = 3, cex = 0.7)





### 5% cut-off analysis ## --------------------------------------------------------------------------

par(mfrow = c(2,2), mar = c(4,4,4,4))

shared_counts <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV/Sub-samples/Shared_winners_5_percent_counts.tsv", header = F)
subs_stats <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYENV/Sub-samples/Sub_sample_stats/Mean_bottleneck_stats_per_subsample.tsv")

names <- paste(shared_counts$V2, shared_counts$V3)
names2 <- gsub("_bottlenecked_pops rep_","/", names)

## Numb identified loci Vs Numb bottlenecked pops.

plot(shared_counts$V1, xaxt = "n", pch = 16, col = "blue", main = "#Loci vs #Bottlenecked pops", xlab = "replicate", ylab = "Number of Co-ID'd loci")
axis(1, names2, at = seq(1, length(names2)), cex.axis = 0.5, las = 2)


## Number IDd loci Vs GW

plot(subs_stats$GW_M_Ratio[3:55], shared_counts$V1, pch = 16, col = "blue", main = "#Loci vs GW index", xlab = "GW index", ylab = "Number of Co-ID'd loci" )
abline(lm(shared_counts$V1 ~subs_stats$GW_M_Ratio[3:55] ), lty = 2, col = "blue")


## Numb loci vs Temp
temps <- read.delim("Sub_sample_stats/Sub_samp_temp_ranges.txt", header = F)
names3 <- gsub("_bottlenecked_pops/rep_","/", temps$V1)

plot(temps$V2[3:55], shared_counts$V1, pch = 16, col = "blue", las = 2, main = "#Loci vs temperature range in subsample", xlab = "Temp range (oC)", ylab = ("Number of Co-ID'd loci"))
#axis(1, names3, at = seq(1, length(names3)), cex = 0.7, las = 2)
abline(lm(shared_counts$V1 ~ temps$V2[3:55]), lty = 2, col = "blue")


av_temps <- read.delim("Ave_temp_per_bottlenecked_class.tsv")
plot(av_temps, pch = 16, col = "blue", main = "Averege temp range in bottleneck class", xlab = "N bottlenecked pops per subsample", ylab = "Temp range (oC)")
abline(lm(av_temps$Av_temp_range ~ av_temps$N_bottled), lty = 2, col = "blue")


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
