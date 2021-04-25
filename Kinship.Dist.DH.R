############################################################
############################################################
## This script contains functions for calculating different kinship and distance matrices for DH lines
## Input is always a individual x marker matrix coded as 0 and 2
############################################################
############################################################


## Testing matrix with 100 individuals x 1000 marker
# geno = matrix(sample(c(0, 2), 100000, replace = TRUE), ncol = 1000)

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


## Calculates the VanRaden kinship matrix for fully homozygous inbred lines (Albrecht et al. 2011, VanRaden 2008, Habier et al. 2007)
VanRaden_kinship = function(W) {
  
  M = ncol(W)
  
  N = nrow(W)
  
  freq = colMeans(W / 2)
  
  p_m = 2 * ifelse(freq <= 0.5, freq, (1 - freq))
  
  P = matrix(rep(p_m, N),
             ncol = M,
             byrow = TRUE)
  
  return(((W - P) %*% t(W - P)) / (8 * sum(p_m * (1 - p_m))))
}

## Calculates the matrix of pairwise modified rogers distances (Reif et al. 2005, Wright 1978, Goodman & Stuber 1983)
mrd_mtx = function(W){
  M = ncol(W)
  
  euclidean_dist = as.matrix(dist(W, method = 'euclidean'))
  return(euclidean_dist/sqrt(M * 4)) 
}


