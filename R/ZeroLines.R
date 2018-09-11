ZeroLines = function(M, V, n_sp, antprob){
  # Check if there are any lines full of zero`s in adjacency matrix
  # Args:
  #  M = Mutualistic adjacency matrix
  #  V = Antagonistic adjacency matrix
  #
  #  Return:
  #    In case of lines full of zeros, generate new M and V matrices
  #
  # This function works with the Antagonize.R function to create adjacency matrices
  source("Antagonize.R") # load Antagonize functiom
  
  m = rowSums(M) # m and v are the number of 1's in the lines of M and V
  v = rowSums(V)
  
  if(any(m == 0)){ # if a certain mutualism line is filled with zeros
    M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
    diag(M) = 0 # no intraespecific interactions
    matrices = Antagonize(M, antprob)
    return(matrices)
  }
  
  if(any(v == 0)){

    
    mats = list(M, V)
    return(mats)
  }
}