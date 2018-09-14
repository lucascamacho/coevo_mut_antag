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
  source("EndInteraction.R") # load EndInteraction function
  
  m = rowSums(M) # m and v are the number of 1's in the lines of M and V
  v = rowSums(V)
  
  if(any(m == 0) | any(v == 0)){ # if a certain mutualism line is filled with zeros
    warning("M or V has a line filled with zeros. Run the code again") # show a warning message
    
    M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
    diag(M) = 0 # no intraespecific interactions
    antagonize = Antagonize(M, antprob) # generate V matrix
    M = antagonize[[1]]
    V = antagonize[[2]]
    
    end = EndInteraction(M, V, "antagonism")
    M = end[[1]]
    V = end[[2]]
    
    matrices = list(M, V)
    return(matrices)
  }
}