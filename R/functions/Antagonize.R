#----------------------------------------------------------------------------------------------#
Antagonize = function(M, antprob){
  # Following a certain probability, transform effects from positive to negative. 
  #
  # Args:
  #   M = matrix of positive effects
  #   antprob = probability of a certain link became negative
  #
  # Return:
  #   mats = list of positive and negative effects, respectively
  
  # create V matrix with the same size of mat
  V = M * 0 
  
  # define the matrix upper triangle
  index = which(upper.tri(M), arr.ind = TRUE)
  
  for(i in 1:nrow(index)){ # for each element in index
    if(M[index[i,][1], index[i,][2]] == 1){ # if this element is equal to 1
      p = runif(1, 0, 1) # sample a number between 0 and 1
      if(p <= antprob){ # if this number is equal or lower than antprob...
        M[index[i,][2], index[i,][1]] = 0 # inverse element is zero (M[j,i])
        V[index[i,][2], index[i,][1]] = 1 # inverse element in V is "on" (V[j,i])
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