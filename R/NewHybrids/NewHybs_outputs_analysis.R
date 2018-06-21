
library(diveRsity)
library(reshape2)
library(ggplot2)


#### RAD NewHybs outputs -----------------------------------------------------------------------------------------------

## CRU_GOLD ----------------------------------------

RAD_Cru_AU <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/NewHybrids/Ref_aligned_analysis/Cru_gold_withref_diag_formatted.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/NewHybrids/Ref_aligned_analysis/Cru_gold_withref_diag_formatted.dat_Newhybs_outputs/Indiv_names_code_replacde.txt", header = F) ## Add sample names

head(RAD_Cru_AU)
RAD_Cru_AU$X500000_sweeps <- indiv.names$V1
RAD_Cru_AU$IndivName <- NULL ## remove column

RAD_Cru_AU_melted <- melt(RAD_Cru_AU, id.vars = "X500000_sweeps") ## reshape the data


Cru_AU_DAN_plot <- ggplot(RAD_Cru_AU_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6))
Cru_AU_DAN_plot


## CRU_GIB ----------------------------------------

RAD_Cru_GIB <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/RAD/CRU_GIB_no_Danube/Newhybs_outs/aa-PofZ.txt")
indiv.names <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/RAD/CRU_GIB_no_Danube/Indiv_names_code_replacde.txt", header = F) ## Add sample names

head(RAD_Cru_GIB)
RAD_Cru_GIB$X301_sweeps <- indiv.names$V1
RAD_Cru_GIB$IndivName <- NULL ## remove column

RAD_Cru_GIB_melted <- melt(RAD_Cru_GIB, id.vars = "X301_sweeps") ## reshape the data


RAD_Cru_GIB_plot <- ggplot(RAD_Cru_GIB_melted, aes(x = X301_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X301_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6))
RAD_Cru_GIB_plot


## CRU_CYP ----------------------------------------

RAD_Cru_Comm <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/NewHybrids/Ref_aligned_analysis/Cru_comm_withref_diag_formatted.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/NewHybrids/Ref_aligned_analysis/Cru_comm_withref_diag_formatted.dat_Newhybs_outputs/Indiv_names_code_replacde.txt", header = F) ## Add sample names

head(RAD_Cru_Comm)
RAD_Cru_Comm$X500000_sweeps <- indiv.names$V1
RAD_Cru_Comm$IndivName <- NULL ## remove column

RAD_Cru_Comm_melted <- melt(RAD_Cru_Comm, id.vars = "X500000_sweeps") ## reshape the data


RAD_Cru_Comm_plot <- ggplot(RAD_Cru_Comm_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6))
RAD_Cru_Comm_plot


### MICROSAT Newhybrids -----------------------------------------------------------------------------------

setwd("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_AU_DANUBE_ONLY/")

## Danube gib only --------------------------

cru_danube <- read.delim("CRU_GIB_DANUBE_ONLY_NEWHYBS.dat_Newhybs_outputs/aa-PofZ.txt")
names <- read.delim("CRU_GIB_DANUBE_ONLY_NEWHYBS.dat_Newhybs_outputs/Indiv_names.txt", header = F)

names(cru_danube) <- c("X500000_sweeps", "IndivName", "X1.000.0.000.0.000.0.000", "X0.000.0.000.0.000.1.000", "X0.250.0.250.0.250.0.250", "X0.000.0.500.0.500.0.000",  "X0.500.0.250.0.250.0.000", "X0.000.0.250.0.250.0.500")

cru_danube$X500000_sweeps <- names$V1
cru_danube$IndivName <- NULL ## remove column

cru_danube_melted <- melt(cru_danube,id.vars = "X500000_sweeps")

q <- ggplot(cru_danube_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity')+
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2",  "yellow","purple", "green", "orange"))+
  theme_bw() 
q + theme(axis.text.x = element_text(angle = 90))



### Crucian Common --------------------------------------------------------------------------------------------------------

Cru_comm_outputs <- read.delim("")
indiv.names <- read.delim("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_CYP_NO_DANUBE/Indiv_names.txt", header = F) ## Add sample names
Cru_comm_outputs$X500000_sweeps <- indiv.names$V1
Cru_comm_outputs$IndivName <- NULL ## remove column

Cru_comm_melted <- melt(Cru_comm_outputs,id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
q <- ggplot(Cru_comm_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity')
q + theme(axis.text.x = element_text(angle = 90))




### Crucian Common --------------------------------------------------------------------------------------------------------

Cru_comm_outputs <- read.delim("Cru_comm_NewHybs_Diag_SNPs.dat_Newhybs_outputs/new_run/aa-PofZ.txt")
indiv.names <- read.delim("Cru_comm_NewHybs_Diag_SNPs.dat_Newhybs_outputs/Indiv.names.formatted", header = F) ## Add sample names
Cru_comm_outputs$X500000_sweeps <- indiv.names$V1
Cru_comm_outputs$IndivName <- NULL ## remove column

Cru_comm_melted <- melt(Cru_comm_outputs,id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
q <- ggplot(Cru_comm_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity')
q + theme(axis.text.x = element_text(angle = 90))

## Crucian goldfish --------------------------------------------------------------------------------------------------------

Cru_gold_outputs <- read.delim("Cru_gold_NewHybs_Diag_SNPs.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("Cru_gold_NewHybs_Diag_SNPs.dat_Newhybs_outputs/Indiv.names.formatted", header = F) ## Add sample names
Cru_gold_outputs$X500000_sweeps <- indiv.names$V1
Cru_gold_outputs$IndivName <- NULL ## remove column

Cru_gold_melted <- melt(Cru_gold_outputs,id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
p <- ggplot(Cru_gold_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity')
p + theme(axis.text.x = element_text(angle = 90))


## Crucian gibel --------------------------------------------------------------------------------------------------------

Cru_gib_outputs <- read.delim("Cru_gib_NewHybs_Diag_SNPs.dat_Newhybs_outputs/new_run/aa-PofZ.txt")
indiv.names <- read.delim("Cru_gib_NewHybs_Diag_SNPs.dat_Newhybs_outputs/Indiv.names.formatted", header = F) ## Add sample names
Cru_gib_outputs$X500000_sweeps <- indiv.names$V1
Cru_gib_outputs$IndivName <- NULL ## remove column

Cru_gib_melted <- melt(Cru_gib_outputs,id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
r <- ggplot(Cru_gib_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity') +
  scale_fill_manual('X500000_sweeps', values = c("red", "black", "dodgerblue2", "yellow", "green", "orange"))+
  theme_bw() 
r + theme(axis.text.x = element_text(angle = 90))



## -- MICROSATS  ------------------------------------------------------------------------------------------


#### REPLOT THIS WHEN IT HAS BEEN RE-RUN!!

### Crucian Common --------------------------------------------------------------------------------------------------------

Cru_comm_outputs <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_CYP_NO_DANUBE/Cru_CYP_NEWHYBS_4Loc_only.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_CYP_NO_DANUBE/Indiv_names_code_replacde.txt", header = F) ## Add sample names
Cru_comm_outputs$X500000_sweeps <- indiv.names$V1
Cru_comm_outputs$IndivName <- NULL ## remove column

head(Cru_comm_outputs)

Cru_comm_melted1 <- melt(Cru_comm_outputs[c(1:236),],id.vars = "X500000_sweeps") ## reshape the data
Cru_comm_melted2 <- melt(Cru_comm_outputs[c(237:553),],id.vars = "X500000_sweeps") ## reshape the data
Cru_comm_melted3 <- melt(Cru_comm_outputs[c(554:756),],id.vars = "X500000_sweeps") ## reshape the data
Cru_comm_melted4 <- melt(Cru_comm_outputs[c(757:1021),],id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
Cru_comm_1 <- ggplot(Cru_comm_melted1, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange", "black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 4))

Cru_comm_2 <- ggplot(Cru_comm_melted2, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange", "black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 4), legend.position="none")

Cru_comm_3 <- ggplot(Cru_comm_melted3, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange", "black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position="none")

Cru_comm_4 <- ggplot(Cru_comm_melted4, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange", "black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position="none")

multiplot(Cru_comm_1, Cru_comm_2, Cru_comm_3, Cru_comm_4, cols = 1)



### Crucian Gold DANUBE ONLY  --------------------------------------------------------------------------------------------------------


Cru_AU_DAN_only <- read.delim("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_AU_DANUBE_ONLY/CRU_AU_DAN_ONLY_Loc_removed/CRU_GIB_DANUBE_ONLY_NEWHYBS_no_J7_or_GF29.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_AU_DANUBE_ONLY/CRU_GIB_DANUBE_ONLY_NEWHYBS.dat_Newhybs_outputs/Indiv_names_code_replacde.txt", header = F) ## Add sample names

head(Cru_AU_DAN_only)
Cru_AU_DAN_only$X500000_sweeps <- indiv.names$V1
Cru_AU_DAN_only$IndivName <- NULL ## remove column


Cru_AU_DAN_only_melted <- melt(Cru_AU_DAN_only, id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
Cru_AU_DAN_plot <- ggplot(Cru_AU_DAN_only_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("dodgerblue2", "red", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6))
Cru_AU_DAN_plot



### Crucian Gold NO DANUBE  --------------------------------------------------------------------------------------------------------

Cru_AU_noDAN <- read.delim("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_AU_NO_DANUBE/CRU_AU_NO_DANUBE_NEWHYBS.txt_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/CRU_AU_NO_DANUBE/Indiv_names_code_replacde.txt", header = F) ## Add sample names

head(Cru_AU_noDAN)
Cru_AU_noDAN$X500000_sweeps <- indiv.names$V1
Cru_AU_noDAN$IndivName <- NULL ## remove column

Cru_AU_noDAN_melted <- melt(Cru_AU_noDAN, id.vars = "X500000_sweeps") ## reshape the data


Cru_AU_noDAN_melted1 <- melt(Cru_AU_noDAN[c(1:312),],id.vars = "X500000_sweeps") ## reshape the data
Cru_AU_noDAN_melted2 <- melt(Cru_AU_noDAN[c(313:599),],id.vars = "X500000_sweeps") ## reshape the data
Cru_AU_noDAN_melted3 <- melt(Cru_AU_noDAN[c(600:867),],id.vars = "X500000_sweeps") ## reshape the data
Cru_AU_noDAN_melted4 <- melt(Cru_AU_noDAN[c(868:1216),],id.vars = "X500000_sweeps") ## reshape the data




Cru_noDAN_1 <- ggplot(Cru_AU_noDAN_melted1, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5))

Cru_noDAN_2 <- ggplot(Cru_AU_noDAN_melted2, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position="none")

Cru_noDAN_3 <- ggplot(Cru_AU_noDAN_melted3, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position="none")

Cru_noDAN_4 <- ggplot(Cru_AU_noDAN_melted4, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position="none")

multiplot(Cru_noDAN_1, Cru_noDAN_2, Cru_noDAN_3, Cru_noDAN_4, cols = 1)




## Crucian goldfish --------------------------------------------------------------------------------------------------------

Cru_gold_outputs <- read.delim("Cru_gold_withref_diag_formatted.dat_Newhybs_outputs/aa-PofZ.txt")
indiv.names <- read.delim("Cru_gold_withref_diag_formatted.dat_Newhybs_outputs/Indiv.names", header = F) ## Add sample names
Cru_gold_outputs$X500000_sweeps <- indiv.names$V1
Cru_gold_outputs$IndivName <- NULL ## remove column

Cru_gold_melted <- melt(Cru_gold_outputs,id.vars = "X500000_sweeps") ## reshape the data

## Plot as stacked bar chart
p <- ggplot(Cru_gold_melted, aes(x = X500000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual('X500000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange")) +
  theme_bw() 
p + theme(axis.text.x = element_text(angle = 90))



### Crucian AU_spp --------------------------------------------------------------------------------------------------------


Cru_AU_outputs <- read.delim("FINAL_DATSET_NEWHYBS/CRU_AU/aa-PofZ.txt")
indiv.names <- read.delim("FINAL_DATSET_NEWHYBS/CRU_AU/Cru_Aur_names.txt", header = F)
Cru_AU_outputs$X10000_sweeps <- indiv.names$V1
Cru_AU_outputs$IndivName <- NULL ## remove column

head(Cru_AU_outputs)

length(Cru_AU_outputs[,1])

Cru_AU_melted1 <- melt(Cru_AU_outputs[c(1:312),],id.vars = "X10000_sweeps") ## reshape the data
Cru_AU_melted2 <- melt(Cru_AU_outputs[c(313:599),],id.vars = "X10000_sweeps") ## reshape the data
Cru_AU_melted3 <- melt(Cru_AU_outputs[c(600:867),],id.vars = "X10000_sweeps") ## reshape the data
Cru_AU_melted4 <- melt(Cru_AU_outputs[c(868:1114),],id.vars = "X10000_sweeps") ## reshape the data

## Plot as stacked bar chart
Cru_AU_1 <- ggplot(Cru_AU_melted1, aes(x = X10000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X10000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 2))

Cru_AU_2 <- ggplot(Cru_AU_melted2, aes(x = X10000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X10000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

Cru_AU_3 <- ggplot(Cru_AU_melted3, aes(x = X10000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X10000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

Cru_AU_4 <- ggplot(Cru_AU_melted4, aes(x = X10000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X10000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

multiplot(Cru_AU_1, Cru_AU_2, Cru_AU_3, Cru_AU_4, cols = 1)


### DANUBE_AU ONLY --------------------------------------------------------------------------------------------------------

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats/ALL_DATA_INDV_PRIORS_ANALYSIS/")


Cru_DAN_outputs <- read.delim("aa-PofZ.txt")
indiv.names <- read.delim("CRU_AU_INC_DANUBE_names.txt", header = F)
Cru_DAN_outputs$X100000_sweeps <- indiv.names$V1
Cru_DAN_outputs$IndivName <- NULL ## remove column

tail(Cru_DAN_outputs)

length(Cru_DAN_outputs[,1])

Cru_DAN_melted1 <- melt(Cru_DAN_outputs[c(1:300),],id.vars = "X100000_sweeps") ## reshape the data
Cru_DAN_melted2 <- melt(Cru_DAN_outputs[c(301:600),],id.vars = "X100000_sweeps") ## reshape the data
Cru_DAN_melted3 <- melt(Cru_DAN_outputs[c(601:900),],id.vars = "X100000_sweeps") ## reshape the data
Cru_DAN_melted4 <- melt(Cru_DAN_outputs[c(901:1250, 1259:1299),],id.vars = "X100000_sweeps") ## reshape the data


## Plot as stacked bar chart
Cru_AU_1 <- ggplot(Cru_DAN_melted1, aes(x = X100000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X100000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3))

Cru_AU_2 <- ggplot(Cru_DAN_melted2, aes(x = X100000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X100000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

Cru_AU_3 <- ggplot(Cru_DAN_melted3, aes(x = X100000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X100000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

Cru_AU_4 <- ggplot(Cru_DAN_melted4, aes(x = X100000_sweeps, y = value,fill=variable)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual('X100000_sweeps', values = c("red", "dodgerblue2", "purple", "yellow", "green", "orange"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 3), legend.position="none")

multiplot(Cru_AU_1, Cru_AU_2, Cru_AU_3, Cru_AU_4, cols = 1)
