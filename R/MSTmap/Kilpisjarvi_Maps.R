#install.packages("ASMap")
library(ASMap)

setwd("/home/djeffrie/Data/RADseq/Kilpisjarvi_maps/MSTmap/KAC1")

## KAC1 Family ============================================================================ ##

## Female map 

# pval 1e-6
Femdat_KAC1 <- read.delim("./Female_map/MST_female_map_input.txt", stringsAsFactors = F) ## load input made by MST_input_maker_cline.py

MstFem_KAC1 <- mstmap.data.frame(Femdat, pop.type = "DH", dist.fun = "kosambi",  ## this is the function that makes the map. See help for all the options
                  objective.fun = "COUNT", p.value = 2e-4, noMap.dist = 30, noMap.size = 5,
                  miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                  return.imputed = TRUE, trace = FALSE)

## Plot and keep largest linkage blocks

Kept_female_KAC1 <- LG_sorter(MstFem_KAC1, 20) ## note, plots can take a while to appear - they are big

# OUPUT marker positions and names ===============================================================================================

LG_writer(Kept_female_KAC1, "./Female_map/KAC1_F_2e4_positions.txt")


## ==================================================================================== ##
## Male map ---------------------------------------------------------------------------
## ==================================================================================== ##

# (same process as above)

Maldat_KAC1 <- read.delim("./Male_map/MST_male_map_input.txt", stringsAsFactors = F)

MstMale_KAC1 <- mstmap(Maldat, pop.type = "DH", dist.fun = "kosambi",
                 objective.fun = "COUNT", p.value = 1e-4, noMap.dist = 30, noMap.size = 1,
                 miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                 return.imputed = TRUE, trace = FALSE)

Kept_male_KAC1 <- LG_sorter(MstMale_KAC1, 25)

Kept_male_KAC1$geno$L4$map

# OUPUT marker positions and names ===============================================================================================

LG_writer(Kept_male_KAC1, "KAC1_M_1e4_positions.txt")


## KAC2 Family ============================================================================ ##
## ======================================================================================== ##

setwd("/home/djeffrie/Data/RADseq/Kilpisjarvi_maps/MSTmap/KAC2/MST_map_haplotypes/")

## Female map 
## ======================================================================================== ##

## pval 1e-6
Femdat_KAC2 <- read.delim("./Female_map/MST_female_map_input.txt", stringsAsFactors = F) ## load input made by MST_input_maker_cline.py

MstFem_KAC2 <- mstmap.data.frame(Femdat, pop.type = "DH", dist.fun = "kosambi",  ## this is the function that makes the map. See help for all the options
                                 objective.fun = "COUNT", p.value = 6e-5, noMap.dist = 30, noMap.size = 5,
                                 miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                                 return.imputed = TRUE, trace = FALSE)



Kept_female_KAC2 <- LG_sorter(MstFem_KAC2, 25) ## note, plots can take a while to appear - they are big

# OUPUT marker positions and names ===============================================================================================

LG_writer(Kept_female_KAC2 , "KAC2_F_6e5_positions.txt")

## ==================================================================================== ##
## Male map ---------------------------------------------------------------------------
## ==================================================================================== ##

# (same process as above)

Maldat <- read.delim("./Male_map/MST_male_map_input.txt", stringsAsFactors = F)

MstMale <- mstmap(Maldat, pop.type = "DH", dist.fun = "kosambi",
                  objective.fun = "COUNT", p.value = 3e-6, noMap.dist = 30, noMap.size = 1,
                  miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                  return.imputed = TRUE, trace = FALSE)

Kept_male <- LG_sorter(MstMale, 25)

Kept_male$geno$L2$map

LG_writer(Kept_male, "KAC2_M_3e6_positions.txt")

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





