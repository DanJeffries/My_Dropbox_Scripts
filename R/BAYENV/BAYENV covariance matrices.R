######### Preliminary run #########
?image
setwd("/home/dan/RAD_programs/BAYENV/bayenv_release")
par(mar = c(4,4,2,2))

## first subset (head -n2000) - ie first 1000 SNPs
test_mat <- read.delim("preliminary_run/old_test_runs/test_BAYENV_sub1_30000iter_matrix.out", header = F)
popnames <- c("BF","BOR", "CAKE", "CALK", "COP", "MOAT", "OBY", "OU", "PED", "POLEN" ,"PRO" ,"SD", "SK", "STEC", "STYV", "TROM", "TU", "V", "WEN")

## looking at the cov matrix
test_mat

## get rid of last column of NAs
test_mat$V20 <- NULL

#convert to matrix
test_mat_mat <- as.matrix(test_mat)


image(test_mat_mat, xaxt = 'n', yaxt = 'n', col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for 1st 1000 snp loci out of 19077") 

## 2nd subset (tail 2000) - ie last 1000 SNPs
test_mat_2 <- read.delim("test_BAYENV_sub2_30000iter_matrix.out", header = F)

test_mat_2

## get rid of last column of NAs
test_mat_2$V20 <- NULL

#convert to matrix
test_mat_mat_2 <- as.matrix(test_mat_2)


image(test_mat_mat_2, xaxt = 'n', yaxt = 'n', col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for last 1000 snp loci out of 19077") 

## These are slightly different.

## Sub 3 - middle 1000 SNPs

test_mat_3 <- read.delim("test_BAYENV_sub3_30000iter_matrix.out", header = F)

test_mat_3

## get rid of last column of NAs
test_mat_3$V20 <- NULL

#convert to matrix
test_mat_mat_3 <- as.matrix(test_mat_3)


image(test_mat_mat_3, xaxt = 'n', yaxt = 'n',col = heat.colors(12)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for middle 1000 snp loci out of 19077") 


### 5000 Randomly sampled SNPs from vcf 30000 iterations ###

test_mat_5000 <- read.delim("BAYENV_format_5000snps_matrix_30000iters.out", header = F)

test_mat_5000

## get rid of last column of NAs
test_mat_5000$V20 <- NULL

#convert to matrix
test_mat_mat_5000 <- as.matrix(test_mat_5000)


image(test_mat_mat_5000, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077") 



### 5000 Randomly sampled SNPs from vcf 25000 iterations - to check for convergence ###

test_mat_5000_25000its <- read.delim("BAYENV_format_5000snps_matrix_25000iters.out", header = F)

test_mat_5000_25000its

## get rid of last column of NAs
test_mat_5000_25000its$V20 <- NULL

#convert to matrix
test_mat_mat_5000_25000its <- as.matrix(test_mat_5000_25000its)


image(test_mat_mat_5000_25000its, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077, 25000its") 

### 5000 Randomly sampled SNPs from vcf 25000 iterations - to check for convergence ###

test_mat_5000_20000its <- read.delim("BAYENV_format_5000snps_matrix_20000iters.out", header = F)

test_mat_5000_20000its

## get rid of last column of NAs
test_mat_5000_20000its$V20 <- NULL

#convert to matrix
test_mat_mat_5000_20000its <- as.matrix(test_mat_5000_20000its)


image(test_mat_mat_5000_20000its, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077, 25000its") 

### Still not happy that its converging ###

### 5000 Randomly sampled SNPs from vcf 95000 iterations - to check for convergence ###

mat_95000its <- read.delim("95000its_matrix.out", header = F)

mat_95000its

## get rid of last column of NAs
mat_95000its$V20 <- NULL

#convert to matrix
mat_mat_95000its <- as.matrix(mat_95000its)


image(mat_mat_95000its, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077, 95000its") 


### 5000 Randomly sampled SNPs from vcf 100000 iterations - to check for convergence ###

mat_100000its <- read.delim("100000its_matrix.out", header = F)

mat_100000its

## get rid of last column of NAs
mat_100000its$V20 <- NULL

#convert to matrix
mat_mat_100000its <- as.matrix(mat_100000its)


image(mat_mat_100000its, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames, at = seq(0,1,(1/18)), las = 2) ## add pop names to axes
axis(2, popnames, at = seq(0,1,(1/18)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077, 100000its") 


########### non-Bottlenecked ##############


### 5000 Randomly sampled SNPs from vcf 100000 iterations - to check for convergence ###

mat_100000its <- read.delim("mat_100000.out", header = F)

mat_100000its

## get rid of last column of NAs
mat_100000its$V10 <- NULL

#convert to matrix
mat_mat_100000its <- as.matrix(mat_100000its)
popnames = read.delim("non_bottlenecked/non_bottlenecked_pop_order.txt", header = F)
popnames
image(mat_mat_100000its, xaxt = 'n', yaxt = 'n',col = heat.colors(24)) ## image function for looking at matrices
axis(1, popnames$V1, at = seq(0,1,(1/8)), las = 2) ## add pop names to axes
axis(2, popnames$V1, at = seq(0,1,(1/8)), las = 2)
title(main = "Cov Matrix for random 5000 snp loci out of 19077, 100000its") 