## Installation of SimRAD dependencies was tricky. Not done through CRAN
library(SimRAD)

#### Estimating the number of ddRAD sites in our species ### ---------------------------------------------------

setwd("~/Data/Genomes/Rtemp/V2/")

rfsq <- ref.DNAseq("RtempV2.0.fa" , subselect.contigs = TRUE, prop.contigs = 0.1) ## make a reference database from the fasta of the assembly.
## With these arguments, the function subsets the contigs (0.1 in this point). It essentially then just appends all the fragments together.

## Restriction enzyme cut sites

## SBF1
cs_5p1 <- "CCTGCA" ## 5 prime section of restriction site
cs_3p1 <- "GG" ## 3 prime section of restriction site

#MSe1
cs_5p2 <- "T"
cs_3p2 <- "TAA"

ref.digest <- insilico.digest(rfsq,  cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE) ## Do the insilico digestion, returns some summary numbers:
## For the Rtemp current assembly version:

# Number of restriction sites for the first enzyme: 12695
# Number of restriction sites for the second enzyme: 3059417
# Number of type AB and BA fragments:24486
# Number of type AA fragments:452
# Number of type BB fragments:3047173

ref.selected <- adapt.select(ref.digest, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2) ## then select only fragments with the AB+BA flanking pattern
## These are fragments that are flanked by both restriction enzyme sites - used for ddRAD.

size_sel_ref.dig <- size.select(ref.selected, min.size = 400, max.size = 500, graph=TRUE, verbose=TRUE) ## Then do the size selection
## In our RAD libraries we use 400-500
## Using this range I got 1407. Multiplying by 10 (because we are only using 0.1 of the contigs in the assembly), gives about 14000 tags.
## It is difficult to compare between libraries as the DNA quality and other things change a lot. But this isn't that close to what we normally find

## But assuming the size selection on the gel isn't very accurate, it could be that we select for closer to 350 - 550
size_sel_ref.dig <- size.select(ref.selected, min.size = 350, max.size = 550, graph=TRUE, verbose=TRUE) 
## In which case the number of tags comes out as about 30,000. Which is much closer to our data. 

## I think with this insilico approach, we are more likely to overestimate the nnumber of tags though, as low coverage positions will drop out of the dataset.
## So all in all, I am not that convinced by this. It might be the quality of the genome that is the problem. Can try with Xenopus but it has a much smaller genome
## and would likely not be that accurate. 



