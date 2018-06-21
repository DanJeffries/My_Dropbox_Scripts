setwd("~/Data/Rana_phylo/Rana_subset/")


library(ape)
library(phytools)

###################### UPDATED FUNCTION THAT ACCEPTS POLYTOMIES #######################
rerootingMethod<-function(tree,x,model=c("ER","SYM"),...){
  if(!is.binary.tree(tree)){
    mt<-tree
    tree<-multi2di(tree)
    multif<-TRUE
  } else multif<-FALSE
  if(hasArg(tips)) tips<-list(...)$tips
  else tips<-NULL
  if(!is.matrix(model)) model<-model[1]
  n<-length(tree$tip.label)
  # if vector convert to binary matrix
  if(!is.matrix(x)){ 
    yy<-to.matrix(x,sort(unique(x)))
    if(is.null(tips)) tips<-FALSE
  } else { 
    if(is.null(tips)) tips<-TRUE
    yy<-x
  }
  yy<-yy[tree$tip.label,]
  yy<-yy/rowSums(yy)
  YY<-apeAce(tree,yy,model=model)
  Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],length(YY$states),length(YY$states),
            dimnames=list(YY$states,YY$states))
  diag(Q)<--colSums(Q,na.rm=TRUE)
  nn<-if(tips) c(1:n,2:tree$Nnode+n) else 2:tree$Nnode+n
  ff<-function(nn){
    tt<-reroot(tree,node.number=nn,position=tree$edge.length[which(tree$edge[,2]==nn)])
    apeAce(tt,yy,model=model,fixedQ=Q)$lik.anc[1,]
  }
  XX<-t(sapply(nn,ff))
  if(tips) XX<-rbind(XX[1:n,],YY$lik.anc[1,],XX[(n+1):nrow(XX),])
  else XX<-rbind(yy,YY$lik.anc[1,],XX)
  rownames(XX)<-1:(tree$Nnode+n)
  if(tips) rownames(XX)[1:n]<-tree$tip.label
  XX<-if(tips) XX else XX[1:tree$Nnode+n,]
  if(multif){
    M<-matchNodes(mt,tree)
    ii<-sapply(M[,2],function(x,y) which(y==x),y=as.numeric(rownames(XX)))
    XX<-XX[ii,]
    rownames(XX)<-M[,1]		
  }
  return(list(loglik=YY$loglik,Q=Q,marginal.anc=XX))
}

#######################################################################################



##Species only

tre<-read.tree("BEAST/100M_MCMC/Rep2/BDNF_3_tree_ann_v1_NEWICK")
plot.phylo(tre)

print(tre$tip.label)  ## number of taxa

chr_sp<-c(8,1,1,2,"?","?",5,"?",1,"?","?",5,5,1,"?","?",1,2,2,"?","?",3,"?",3,3)
mat_chr2<-read.table("chr_matrix.txt", h=T)
mat_chr<-as.matrix(mat_chr2[,3:7])
row.names(mat_chr)<-tre$tip.label
colnames(mat_chr)<-c("1","2","3","5","8")
names(chr_sp)<-tre$tip.label
#write.table(chr_sp, "chr.txt", row.names=T, col.names=F)

anc_sp <- rerootingMethod(tre, mat_chr, model = "ER")

tiplabels(pie = mat_chr, cex=0.7)#to.matrix(mat_chr, sort(unique(mat_chr))), cex = 1)
nodelabels(node = as.numeric(rownames(anc_sp$marginal.anc)), pie = anc_sp$marginal.anc, cex = 0.6, col = couleurs)
nodelabels(names(anc_sp$marginal.anc))
couleurs<-c("#189e77", "#d86013", "#7470b2", "#e52d89", "#e6ab01")
for(i in 1:length(anc_sp$marginal.anc)) {
  pie(anc_sp$marginal.anc[i,], labels="", main=i, col=couleurs)
}
##

mat_chr_with_undiff2<-read.table("chr_matrix_with_undifferetiated.txt", h=T)
mat_chr_with_undiff<-as.matrix(mat_chr_with_undiff2[,3:8])
row.names(mat_chr_with_undiff)<-tre$tip.label
colnames(mat_chr_with_undiff)<-c("1","2","3","5","8","undiff")

##







## Separated -- you can forget about it, I didn't use it in the end

tree<-read.tree("rerooted_rana_poly_separated.newick.txt")
plot.phylo(tree, show.node.label=T, type="phylogram", root.edge=T)

het<-c("M","F","M","?","M","M","M","M","?","?","?","M","?","M","M","?","M","?","?","M","M","M","?","?","M","M","M","?","M","M","?","M","M","?","M","M","M","?")
names(het)<-tree$tip.label

chr<-as.character(c(8,8,1,"?",1,"1_2",1,"2_3","?","?","?",5,"?",1,3,"?",1,"?","?",5,5,1,"?","?",1,1,2,"?",2,"?","?",3,"?",3,"?",3,5,"?"))
names(chr)<-tree$tip.label