ZeroLines = function(M, V, n_sp, antprob){
  # Check if there are any lines full of zero's in M adjacency matrix
  # Args:
  #  M = Mutualistic adjacency matrix
  #  V = Antagonistic adjacency matrix
  #  n_sp = number os species
  #  antprob = probability of a link become antagonist
  #
  #  Return:
  #    In case of lines full of zeros, make new adjacency matrices
  #
  # Obs: this function works with the Antagonize.R and EndInteraction.R functions.
  #      Make sure to have then in your work directory

  # load packages for new matrices if necessary
  source("Antagonize.R")
  source("EndInteraction.R")
  
  # try several types get a M matrix without lines of zero's
  for(z in 1:1000000){
    m = rowSums(M) # check for zero lines
    
    # make new matrices if there's a line full of zero's in M
    if(any(m == 0)){ 
      M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
      diag(M) = 0 # no intraespecific interactions
      
      # Antagonize M
      antagonize = Antagonize(M, antprob)
      M = antagonize[[1]]
      V = antagonize[[2]]
      
      # End antagonism AA
      end = EndInteraction(M, V, "antagonism")
      M = end[[1]]
      V = end[[2]]     
    }
    
    # if there are no zero lines, stop
    else{
      break
    }
  }
  
  # return the current M and V matrices
  mats = list(M, V)
  return(mats)
}