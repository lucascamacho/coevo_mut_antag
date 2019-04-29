#---------------------------------------------------------------------------------------------------#
Antagonize = function(M, antprob){
  # Following a certain probability, transform interactions outcomes from positive to negative. 
  #
  # Args:
  #   M = matrix of positive interaction outcomes
  #   antprob = probability of a certain link became a negative link
  #
  # Return:
  #   mats = list of positive and negative outcomes, respectively
  # 
  V = M * 0 # create V matrix with the same size of mat
  
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      if(M[i,j] == 1){
        p = runif(1, 0, 1)
        
        if(p <= antprob){
          M[j,i] = 0
          V[j,i] = 1
        }  
      }
    }
  }
  #P = matrix(runif(dim(M)[1]*dim(M)[2], min = 0, max = 1),  # create matrix of probabilities
  #            ncol = dim(M)[1], nrow = dim(M)[2])

  #M[antprob >= P] = 0 # if the probability is lower or equal than antprob,the link in M in zero
  #V[antprob >= P] = 1 # and the link in V is 1
  diag(V) = 0 # diagonal's matrices must be zero

  mats = list(M, V) # create and return a list with the positive and negative matrices
  return(mats)
}
#---------------------------------------------------------------------------------------------------#
