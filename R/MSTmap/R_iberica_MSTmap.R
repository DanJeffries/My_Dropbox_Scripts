install.packages("ASMap")

library(ASMap)
setwd("/home/djeffrie/Data/Ribe_LM/MSTmap_haplotypes/V2_genome_mapping/")
setwd("/path/to/MST_maker_for_glib/")

## Female map -------------------------------------------------------------------------

Femdat <- read.delim("./Female_map/MST_female_map_input.txt", stringsAsFactors = F) ## load input made by MST_input_maker_cline.py

## Make the linkage map
MstFem <- mstmap.data.frame(Femdat, pop.type = "DH", dist.fun = "kosambi",
                  objective.fun = "COUNT", p.value = 5e-05, noMap.dist = 30, noMap.size = 1,
                  miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                  return.imputed = TRUE, trace = FALSE)

## Plot the linkage map
Kept_female <- LG_sorter(MstFem, 4) ## see my function at bottom of this script for description


## Make the linkage map without compli loci

Femdat2 <- read.delim("./Female_map/MST_female_map_input_wo_compli.txt", stringsAsFactors = F) ## load input made by MST_input_maker_cline.py


MstFem2 <- mstmap.data.frame(Femdat, pop.type = "DH", dist.fun = "kosambi",
                            objective.fun = "COUNT", p.value = 5e-05, noMap.dist = 30, noMap.size = 1,
                            miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                            return.imputed = TRUE, trace = FALSE)



## Plot the linkage map
Kept_female2 <- LG_sorter(MstFem2, 4) ## see my function at bottom of this script for description








## Write the marker names from a chosen linkage group to a file for further analyses.
write(names(MstFem$geno$L18$map), file = "./Female_map/Fem_LG18_tag_IDs.txt")


## Male map ---------------------------------------------------------------------------
# (same process as above)

Maldat <- read.delim("./Male_map/MST_male_map_input.txt", stringsAsFactors = F)

MstMale <- mstmap(Maldat, pop.type = "DH", dist.fun = "kosambi",
                 objective.fun = "COUNT", p.value = 5e-06, noMap.dist = 30, noMap.size = 1,
                 miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE, as.cross = TRUE,
                 return.imputed = TRUE, trace = FALSE)

Kept_male <- LG_sorter(MstMale, 4)

heatMap(Kept_male, lmax = 50)

## Can take just a subset of the map to plot

kept_LM <- subset(MstMale, chr = 'L4')
heatMap(kept_LM, lmax = 50)

## Again, output marker names from interesting linkage groups. 
write(names(MstMale$geno$L4$map), file = "Male_LG4_tag_IDs.txt")

mylabs<-row.names(Maldat)

str(kept_LM)

MST_combined <- combineMap(Kept_male, Kept_female)

heatMap(MST_combined, lmax = 50)

Kept_combined<- LG_sorter(MST_combined, 4)





mstmap.cross(Kept_combined)


## Dan's LG sorter function -------------------------

### This function is a modification of the ASmap plotting function. It takes an MST object created
### and a threshold number of SNPs provided by the user. If linkage groups are kept and plotted only if
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



