ZeroLines = function(M, V, n_sp, antprob){
  # Check if there are any lines full of zero`s in adjacency matrix
  # Args:
  #  M = Mutualistic adjacency matrix
  #  V = Antagonistic adjacency matrix
  #  n_sp = number os species
  #  antprob = probability of a link become antagonist
  #  interaction = interaction to be ended
  #
  #  Return:
  #    In case of lines full of zeros, make new adjacency matrices
  #
  #check the lines
  source("Antagonize.R")
  source("EndInteraction.R")
  
  for(z in 1:1000000){
    m = rowSums(M)
    v = rowSums(V)
    
    if(any(m == 0) | any(v == 0)){
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
    
    else{
      break
    }
  }
  
  mats = list(M, V)
  return(mats)
}