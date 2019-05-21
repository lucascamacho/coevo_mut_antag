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
  
  # get the indexes of interaction outcomes
  ints = which(lower.tri(M) == 1, arr.ind = TRUE)
  
  # invert index to get the negative outcome position
  ints_2 = ints[ ,c("col", "row")]
  
  # sample to change interaction outcomes
  P = matrix(matrix(runif(nrow(ints), min = 0, max = 1), 
                    ncol = 1, nrow = nrow(ints)))
  
  # the the position of the outcomes that will change
  position = which(P <= antprob)
  l = c(ints_2[position, ][,1])
  c = c(ints_2[position, ][,2])
  
  # change interaction outcomes
  for(i in 1:length(l)){
    M[l[i],c[i]] = 0
    V[l[i],c[i]] = 1
  }
  
  # diagonal's matrices must be zero
  diag(V) = 0
  
  # create and return a list with the positive and negative matrices
  mats = list(M, V) 
  return(mats)
}
#---------------------------------------------------------------------------------------------------#
