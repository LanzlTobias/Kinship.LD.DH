############################################################
############################################################
## This script contains functions for calculating different kinship and distance matrices for DH lines
## Input is always a individual x marker matrix coded as 0 and 2
############################################################
############################################################


## Testing matrix with 100 individuals x 1000 marker
# geno = matrix(sample(c(0, 2), 100000, replace = TRUE), ncol = 1000)

## Recodes genotypes to 0, 1, and 2 and sorts them after the minor allele frequency
recode_geno = function(geno){
  geno_centered = geno - min(geno)
  geno_scaled = geno_centered * (2 / max(geno_centered))
  
  
  vapply(1:ncol(geno_scaled), function(x){
    if(mean(geno_scaled[, x]) <= 1){
      return(geno_scaled[, x])
    }  else {
      (geno_scaled[, x]-2)*(-1)}
  },
  FUN.VALUE = numeric(nrow(geno)))
}

## Calculates the simple matching coefficient between an lines (Sneath & Sokal 1973)
SM_mtx = function(W) {
  M = ncol(W)
  
  N = nrow(W)
  
  J_NxM = matrix(rep(1, N * M), ncol = M)
  
  J_NxN = matrix(rep(1, N * N), ncol = N)
  
  return(((W - J_NxM) %*% t(W - J_NxM) + M * J_NxN) / (2 * M))
}

## Calculates the realized genomic kinship matrix based on the simple matching coefficient (Hayes & Goddard 2008)
SM_kinship = function(W) {
  M = ncol(W)
  
  N = nrow(W)
  
  J_NxM = matrix(rep(1, N * M), ncol = M)
  
  J_NxN = matrix(rep(1, N * N), ncol = N)
  
  SM_matrix = ((W - J_NxM) %*% t(W - J_NxM) + M * J_NxN) / (2 * M)
  
  SM_min = min(SM_matrix)
  
  return((SM_matrix - SM_min) / (1 - SM_min))
}


## Calculates the VanRaden kinship matrix for fully homozygous inbred lines 
VanRaden_kinship = function(W) {
  
  # number of markers
  M = ncol(W)
  
  # number of individuals
  N = nrow(W)
  
  # minor allele frequencies * 2
  maf = colMeans(W)
  
  P = matrix(rep(maf, N),
             ncol = M,
             byrow = TRUE)
  
  # calculate Z
  Z <- W - P
  
  # calculate kin
  return(2*(Z %*% t(Z)) / sum(maf * (2 - maf)))
}

## Calculates the matrix of pairwise modified rogers distances (Reif et al. 2005, Wright 1978, Goodman & Stuber 1983)
mrd_mtx = function(W){
  M = ncol(W)
  
  euclidean_dist = as.matrix(dist(W, method = 'euclidean'))
  return(euclidean_dist/sqrt(M * 4)) 
}


