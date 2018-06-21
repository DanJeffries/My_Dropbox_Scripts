library(ggplot2)
library(diveRsity)
library(scales)

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Thesis_Incremental_runs")

tag_share <- read.delim("Total_and_shared_catalog_tags.txt")

ref_tag_share <- read.delim("reference_aligned/Total_and_shared_catalog_tags.txt")


Colorsv <- c("lightblue", "darkblue", "darkred", "firebrick1", "orange", "darkgreen", "green", "lightgreen")

p <- ggplot(tag_share, aes(N)) + 
  geom_line(aes(y = Total_tags, colour = Catalog), size = 1, linetype = "twodash") + 
  geom_line(aes(y = Shared_tags, colour = Catalog), size = 1) +
  scale_colour_manual(name="Catalog",values=Colorsv)
  
tag_share_stacked <- read.delim("Total_and_shared_catalog_tags_Rformatted_stacked_barplot.txt")

t <- ggplot(tag_share_stacked, aes(N, tags, fill = tot_shar)) +
  geom_bar(position = "identity", stat = "identity",  colour ="black", aes(fill = tot_shar ))+
  facet_wrap(~ Catalog)+
  scale_fill_manual( values = grey.colors(2))+
  scale_y_continuous(labels = comma)+
  scale_x_continuous(labels = comma,breaks = seq(0,10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))
t

ggplot(ref_tag_share, aes(x= Catalog, y = tags, fill = tot_shar))+
  geom_bar(stat = "identity", position = "identity", colour = "black")+
  scale_fill_manual(values = grey.colors(2))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


q <- ggplot(tag_share_stacked, aes(Catalog, tags, fill = tot_shar)) +
  geom_bar(position="fill", stat = "identity")+
  ylab("Percent")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


r <- ggplot(ref_tag_share, aes(Catalog, tags, fill = tot_shar)) +
  geom_bar(position="fill", stat = "identity")+
  ylab("Percent")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

multiplot(q,r, cols = 2)







