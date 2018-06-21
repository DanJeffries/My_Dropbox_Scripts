#install.packages("ape")
#install.packages("phytools")
library(ape)
library(phytools)

setwd("~/Dropbox/My_Dropbox_Scripts/R/Anc_states_Rana/")

## make the tree (species only)

tre<-read.tree("BDNF_3_tree_ann_v1_NEWICK")
plot.phylo(tre)

print(tre$tip.label)  ## get names of taxa

## read in the state probability matrix
Five_states_matrix<-as.matrix(read.table("Tips_prob_matrix_5_states.txt", h=T, row.names = 1))
rowSums(Five_states_matrix)  ## check sum of rows (should add up to 1)

## plot tip pies (i.e. input state matrix)

tiplabels(pie=Five_states_matrix,cex=0.3)

## Stochastic mapping - 5 states

SM_5_states_1000sim <- make.simmap(tre, Five_states_matrix, model = "ARD", nsim =1000, tips = T)
sum_SM_5_states_1000sim <- summary(SM_5_states_1000sim, plot = F)
plot(sum_SM_5_states_1000sim) ## Plots piecharts at nodes summarising all simulations. 

## Stochastic mapping - 10 states

# read in the 10 state probability matrix
Ten_states_matrix<-as.matrix(read.table("Tips_prob_matrix_10_states.csv", h=T, row.names = 1))
rowSums(Five_states_matrix)  ## check sum of rows (mine add up to 1.002 in some cases but still works)
tiplabels(pie=Ten_states_matrix,cex=0.3)

SM_10states_1000sim <- make.simmap(tre, Ten_states_matrix, model = "ARD", nsim =1000, tips = T)
sum_SM_10states_1000sim <- summary(SM_10states_1000sim, plot = F)
plot(sum_SM_10states_1000sim)


############################
## PLOTTING ##
############################

## Output edge data for plotting in python

tree_counter = 0
for (tree in SM_100_10states) {
  tree_counter <- tree_counter+1
  edge_counter = 0
  for (edge in SM_100_10states[[tree_counter]]$maps) {
    edge_counter <- edge_counter + 1
    state_counter = 0
    write(paste0("tree ", tree_counter), "Ten_states_node_transitions.txt")
    write(paste0("edge ", edge_counter), "Ten_states_node_transitions.txt")
    for (state in SM_100_10states[[tree_counter]]$maps[[edge_counter]]){
      state_counter <- state_counter + 1
      write(paste0(names(SM_100_10states[[tree_counter]]$maps[[edge_counter]][state_counter])," ", SM_100_10states[[tree_counter]]$maps[[edge_counter]][[state_counter]]), "Ten_states_node_transitions.txt")
    }
  }
}
