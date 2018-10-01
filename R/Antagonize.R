#---------------------------------------------------------------------------------------------------#
Antagonize <- function(M, antprob){
  # Following a certain probability, transform interactions from mutualism to antagonisms. 
  #
  # Args:
  #   M = square binary matrix of mutualistic interactions
  #   antprob = probability of a certain link became a antagonistic link
  #
  # Return:
  #   intr = list of mutualistic and antagonistic networks, respectively
  # 
  V = M * 0 # create V matrix with the same size of mat
  P = matrix(runif(dim(M)[1]*dim(M)[2], min = 0, max = 1),  # create matrix of probabilities
              ncol = dim(M)[1], nrow = dim(M)[2])

  M[antprob >= P] = 0 # if the probability is lower or equal than antprob,the link in M in zero
  V[antprob >= P] = 1 # and the link in V is 1
  diag(V) = 0 # diagonal's matrices must be zero

  mats = list(M, V) # create list with the antagonistic and mutualistic matrices
  return(mats) # return a list with the antagonistic and mutualistic matrices
}
#---------------------------------------------------------------------------------------------------#
