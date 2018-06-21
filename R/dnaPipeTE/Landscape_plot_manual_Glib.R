library(ggplot2)

### Using only one illumina shotgun read file (Spain_0_)

setwd("~/Data/Transposable_elements/Nanorana/0.2/dnaPipeTE_outs/")

#read the table with reads and identity selected for landscape analysis
file1 = "reads_landscape"
file2 = "factors_and_colors"

land = read.table(file1)


#read.table("reads_landscape")->land
names(land)=c("id", "annot", "fam1", "fam")
land$div=100-land$id

#read the corresponding factor order and color table
fac_col = read.table(file2)

#order factors and colors
land$fam1=factor(land$fam1, levels=as.character(fac_col$V1))

#plot the landscape graph
ggplot(land, aes(div, fill=fam1)) +
      geom_histogram(binwidth=1.1) +
      labs(list(x="div", y="reads count")) + 
      scale_x_continuous(limits = c(0, 35)) + 
      scale_fill_manual(values=as.character(fac_col$V3)) + 
      guides(fill=guide_legend(ncol=3)) + 
      theme(legend.direction ="vertical", legend.position = c(0.75, 0.55))
      
ggsave("landscape.pdf", height=12, width=15)


