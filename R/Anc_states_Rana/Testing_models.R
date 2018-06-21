## testing the transition model liklihoods
?ace  # <- Some model details in here. 

anc_sp_ER <- rerootingMethod(tre, Five_states_matrix, type = "discrete", model = "ER", method = "ML", tips = TRUE) ## model options are "ER", "SYM" or "ARD" - see ?ace for details
anc_sp_SYM <- rerootingMethod(tre, Five_states_matrix, type = "discrete", model = "SYM", tips = TRUE) ## model options are "ER", "SYM" or "ARD" - see ?ace for details
anc_sp_ARD <- rerootingMethod(tre, Five_states_matrix, type = "discrete", model = "ARD", tips = TRUE) ## model options are "ER", "SYM" or "ARD" - see ?ace for details

## Get the number of parameters (k) simply count number of different transition rates

anc_sp_ER$Q
anc_sp_SYM$Q
anc_sp_ARD$Q

nPars_ER = 1
nPars_SYM = 15
nPars_ARD = 30

## Calculate AIC for each model

AIC_ER = 2*1-(2*anc_sp_ER$loglik)
AIC_SYM = 2*15-(2*anc_sp_SYM$loglik)
AIC_ARD = 2*30-(2*anc_sp_ARD$loglik)

AIC_ER
AIC_SYM
AIC_ARD