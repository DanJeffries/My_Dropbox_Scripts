setwd("/home/djeffrie/Data/Pperezi/Stacks_outs")

nsamps <- read.delim("n_samples.temp", header = F)
head(nsamps)

hist(nsamps$V1, breaks = 43)

## Data comes from the bash one-liner in the Pnig evernote

tags_per_samp <- read.delim("tags_per_samp.txt", header = F, sep = " ")
reads_per_samp <- read.delim("../demultiplexed/raw_reads_per_sample.txt", header = F)
par(mar = c(10,2,2,2), mfrow = c(2,1))

barplot(tags_per_samp$V2, names = tags_per_samp$V1, las = 2)
barplot(reads_per_samp$V1, names = reads_per_samp$V2, las = 2)

cd?ordered

setwd("/home/djeffrie/Data/Pnig_RAD/Stacks_outs/populations_1pop_kept_samples_outs")

all_freqs <- read.delim("batch_1.vcf.All_frequency_ratios.tsv")
head(all_freqs)
hist(all_freqs$female.male_freq, breaks = 100)
abline(v = mean(all_freqs$female.male_freq), col = "red", lty = 2)


