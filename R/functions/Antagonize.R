#----------------------------------------------------------------------------------------------#
Antagonize = function(M, antprob){
  # Following a certain probability, transform interactions outcomes from positive to negative. 
  #
  # Args:
  #   M = matrix of positive interaction outcomes
  #   antprob = probability of a certain link became a negative link
  #
  # Return:
  #   mats = list of positive and negative outcomes, respectively
  
  # create V matrix with the same size of mat
  V = M * 0 
  
  for(i in 1:dim(M)[1]){ 
    for(j in 1:dim(M)[2]){
      # for each element of M equal 1
      if(M[i,j] == 1){
        # sample a random value 
        p = runif(1, 0, 1)
        # if this value is lower than antprob
        if(p <= antprob){
        # transform positive effect in negative
        M[j,i] = 0
        V[j,i] = 1
        }
      }
    }
  }

  # no intraespecific interactions
  diag(V) = 0
  
  # create and return a list with the positive and negative matrices
  mats = list(M, V) 
  return(mats)
}
#---------------------------------------------------------------------------------------------------#
