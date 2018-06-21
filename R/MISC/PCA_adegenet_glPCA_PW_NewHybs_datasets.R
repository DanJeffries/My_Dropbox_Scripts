library(adegenet)

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M4_m4_N4/NewHybrids_PW_datafiles/")

# --------DATASETS-------------------------------------------------------------------------------------------

Cru_comm <- read.PLINK('Cru_comm_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_comm_names <- read.delim("Cru_comm_pop.names", header = F)
pop(Cru_comm) <- Cru_comm_names$V1

Cru_gold <- read.PLINK('Cru_gold_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gold_names <- read.delim("Cru_gold_pop.names", header = F)
pop(Cru_gold) <- Cru_gold_names$V1

Cru_gib <- read.PLINK('Cru_gib_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gib_names <- read.delim("Cru_gib_pop.names", header = F)
pop(Cru_gib) <- Cru_gib_names$V1

## -------- PCAs -----------------------------------------------------------------------------------------------

