## Sex_sorter function --------------------------------------------------------------------------------------------

## This function will take a PCA and plot it showing males and females (or whatever
## pops are specified in the sex_info file), in different colours for comparison. 
## The Sex_sorter_plus function takes a third collumn in the Sex_info file, which specifies a
## binary character (e.g. 1 or 0) with which to plot black lines around points. 

Sex_sorter_plus <- function(pca, sex_info, PCs_to_plot = c(1,2), Title = "PCA",  xtitle = "X", ytitle = "Y", leg_title = "Legend", Palette = c('red','royalblue','green','violet','black'), func_cex = 2)
{
  ## Make sure sexes and PC data are in the right order
  
  pca_ordered = pca$scores[order(row.names(pca$scores)),]
  pca_df_ordered <- as.data.frame(pca_ordered) ## ggplots wants a dataframe, not a matrix
  
  sexes_ordered = sex_info[order(sex_info$V1),]
  
  ## Data for plotting
  
  X = pca_df_ordered[,PCs_to_plot[1]]
  Y = pca_df_ordered[,PCs_to_plot[2]]
  
  Data = data.frame(Sample_names = sexes_ordered$V1, X = X, Y = Y, Sex = sexes_ordered$V2, Clones = sexes_ordered$V3)

  ## Labels
  Sample_names = sexes_ordered$V1
  
  a = ggplot(Data, aes(x = X,
                       y = Y,
                       fill = Sex,
                       color = factor(Clones))) + 
    
    geom_point(shape = 21, size = func_cex*5 ,alpha = 0.4, stroke = 1) +

    scale_fill_manual(name = "Sex",
                      values=Palette) +
    scale_color_manual(name = "Clones",
                       values=c("white", "black")) +
    
    annotate("text", x = X, y = Y, label = Sample_names, size = func_cex) +
    
    geom_hline(aes(yintercept=0), color = "grey50") +
    geom_vline(aes(xintercept=0), color = "grey50") +
    ggtitle(Title) + 
    xlab(xtitle) +
    ylab(ytitle) +
    theme(plot.title = element_text(size = 20, face="bold"), 
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey90")) + 
    labs(colour=leg_title)
  
  ggsave(paste(Title, ".pdf", sep = ""), device = "pdf")
  print(a)
  
}

