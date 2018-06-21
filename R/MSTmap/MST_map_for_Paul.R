#install.packages("ASMap")
library(ASMap)

setwd("/home/djeffrie/Data/RADseq/R_temp_fams/Populations_ST01/MSTmap/")

## FAMILY ST01 ============================================================================ ##
## ======================================================================================== ##
## Female map 
## ======================================================================================== ##

## pval 1e-6
Femdat <- read.delim("./Female_map/MST_female_map_input.txt", stringsAsFactors = F) ## load input made by MST_input_maker_cline.py

MstFem_1e06 <- mstmap.data.frame(Femdat, pop.type = "DH", dist.fun = "kosambi",  ## this is the function that makes the map. See help for all the options
                  objective.fun = "COUNT", p.value = 1e-6, noMap.dist = 30, noMap.size = 1,
                  miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                  return.imputed = TRUE, trace = FALSE)
MstFem_1e06$geno

## play with p-value to get optimal loci clustering. Looking for 13 LGs with good numbers of markers and 
## a big difference between them and the "leftovers". 

#=======================================================================================================
## Plot the linkage map (using the function at the bottom of this script (run those lines first to define function))
## see my function description of what it does
Kept_female_1e06 <- LG_sorter(MstFem_1e06, 4) ## note, plots can take a while to appear - they are big

#=======================================================================================================

## If you need to you can write the marker names from a chosen linkage group to a file for further analyses.
write(names(MstFem$geno$L18$map), file = "./Female_map/Fem_LG18_tag_IDs.txt") ## just change the name of the linkage group

## ==================================================================================== ##
## Male map ---------------------------------------------------------------------------
## ==================================================================================== ##

# (same process as above)

Maldat <- read.delim("./Male_map/MST_male_map_input.txt", stringsAsFactors = F)

MstMale <- mstmap(Maldat, pop.type = "DH", dist.fun = "kosambi",
                 objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 30, noMap.size = 1,
                 miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                 return.imputed = TRUE, trace = FALSE)

Kept_male <- LG_sorter(MstMale, 15)


##### ======= Other functions ========== ####

## Can take just a subset of the map to plot

kept_LM <- subset(MstMale, chr = 'L4')
heatMap(kept_LM, lmax = 50)

## Again, output marker names from interesting linkage groups. 
write(names(MstMale$geno$L4$map), file = "Male_LG4_tag_IDs.txt")

## Write marker positions to a file (function at bottom of this script )
LG_writer(Kept_MstFem_EcoRI_1e09, "FemaleMap_marker_positions.txt") 


##########################################################################################
#########################   FUNCTIONS  ###################################################
##########################################################################################


## Dan's LG sorter function ==============================================================
#=========================================================================================

### This function is a modification of the ASmap plotting function. It takes an MST object created
### and a threshold number of SNPs provided by the user. Linkage groups are kept and plotted only if
### they contain more SNPs than the threshold number. 

LG_sorter <- function(MST_obj,N_snp_thresh)
{
  ## A couple of cool functions to get odd/even = TRUE/FALSE
  is.even <- function(x) x %% 2 == 0
  is.odd <- function(x) x %% 2 != 0
  
  for(i in names(MST_obj$imputed.geno))
  {
    i_int <- is.even(strtoi(strsplit(i, "L")[[1]][2]))
    N_snps <- length(names(MST_obj$geno[[i]]$map))
    if(N_snps >= N_snp_thresh)
    { 
      if(i_int == T)
      {
        if(exists("Kept_LGs") == FALSE)
        {
          Kept_LGs = i
        }
        else
        {
          Kept_LGs = c(Kept_LGs, i)
        }
      }
    }
  }
  print(c("Number of kept linkage groups:",length(Kept_LGs)))
  
  kept_LM <- subset(MST_obj, chr = Kept_LGs)
  heatMap(kept_LM, lmax = 50)
  rm("Kept_LGs")
  
  return(kept_LM)
}

## FUNCTION TO WRITE MARKER POSITIONS TO A FILE ========================================================
#=======================================================================================================

LG_writer <- function(MST_obj, file_name)
{
  df = data.frame()
  for (i in names(MST_obj$geno)) {
    print(i)
    df = rbind(df, data.frame(marker = names(MST_obj$geno[[i]]$map), LG = rep(i, length(MST_obj$geno[[i]]$map)), pos = unname(MST_obj$geno[[i]]$map)))
  }
  write.table(df, file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}





