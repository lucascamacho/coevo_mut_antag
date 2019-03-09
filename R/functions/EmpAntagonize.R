#-----------------------------------------------------------------------------------------------------#
EmpAntagonize = function(M, antprob){
  # Following a certain probability, transform interactions of empirical networks 
  # from mutualism to antagonisms. 
  # These networks need to be squared, otherwise, these funciton will return an error
  #
  # Args:
  #   M = empirical and square adjacency matrix of mutualisms
  #   antprob = probability of a link become antagonist
  # Return:
  #   mats = list of mutualistic adjacency matrix and antagonistic adjacency matrix, respectively
  # 
  V = M * 0 # create V matrix
  
  int = M == 1 # where are the interactions?
  
  P = matrix(0, ncol = dim(M)[2], nrow = dim(M)[1]) # create matrix of probabilities
  P[int] = runif(length(M[int]), 0, 1) # take a number from 0 to 1 only for the interactions of M

  M[antprob >= P] = 0 # turn mutualism off 

  V[antprob >= P] = 1 # transform links to antagonism
  V[int == FALSE] = 0 # ignore the zero's
  
  mats = list(M, V) # create list with the antagonistic and mutualistic matrices
  return(mats) # return a list with the antagonistic and mutualistic matrices
  
}
#-----------------------------------------------------------------------------------------------------#