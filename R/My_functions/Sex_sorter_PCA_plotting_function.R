
## Sex_sorter function --------------------------------------------------------------------------------------------

## This function will take a PCA and plot it showing males and females (or whatever
## pops are specified in the sex_info file), in different colours for comparison. 

Sex_sorter <- function(pca, sex_info, PCs_to_plot = c(1,2), Title = "PCA",  xtitle = "X", ytitle = "Y", leg_title = "Legend", Palette = c('red','royalblue','green','violet','black'), func_cex = 2)
{
  ## Make sure sexes and PC data are in the right order
  
  pca_ordered = pca$scores[order(row.names(pca$scores)),]
  pca_df_ordered <- as.data.frame(pca_ordered) ## ggplots wants a dataframe, not a matrix
  
  sexes_ordered = sex_info[order(sex_info$V1),]

  pca_ordered
  sexes_ordered
  
  ## Data for plotting
  
  X = pca_df_ordered[,PCs_to_plot[1]]
  Y = pca_df_ordered[,PCs_to_plot[2]]
  
  ## Labels
  Sample_names = sexes_ordered$V1
  Sex = sexes_ordered$V2
  
  ## Plotting
  
  a = ggplot(pca_df_ordered, aes(X, Y, colour=Sex)) + 
    geom_point(size = func_cex*5 ,alpha = 0.4) +
    geom_hline(aes(yintercept=0), color = "grey50") +
    geom_vline(aes(xintercept=0), color = "grey50") +
    scale_colour_manual(values=Palette) +
    annotate("text", x = X, y = Y, label = Sample_names, size = func_cex)+
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


