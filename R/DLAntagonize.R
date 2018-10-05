DLAntagonize = function(M, antprob){
  V = M * 0 # create V matrix with the same size of mat
  P = matrix(runif(dim(M)[1]*dim(M)[2], min = 0, max = 1),  # create matrix of probabilities
             ncol = dim(M)[1], nrow = dim(M)[2])
 
  M[antprob >= P] = 0 # if the probability is lower or equal than antprob,the link in M in zero
  V[antprob >= P] = 1 # and the link in V is 1
  diag(V) = 0 # diagonal's matrices must be zero

  ind = length(which(V == 1)) # quantos antagonismos tem
  dex = sample(which(V == 1), ind * 0.15)
  M[dex] = 1
  
  mats = list(M, V) # create list with the antagonistic and mutualistic matrices
  return(mats) # return a list with the antagonistic and mutualistic matrices
  
}

