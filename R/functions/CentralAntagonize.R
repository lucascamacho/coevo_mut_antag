CentralAntagonize = function(M){
  # Identify the central species in the network (highest degree),
  # transform all the central species interaction outcomes in negative
  # and count the number of negative outcomes in the network.
  #
  # This funtion is used in empirical networks of interactions. 
  #
  # Args:
  #   M = matrix of positive interaction outcomes
  #
  # Return:
  #   mats = list of positive and negative outcomes, respectively
  #   p_vl = estimated p (antprob) value
  
  # create V matrix with the same size of mat
  V = M * 0
  
  # identify the central species
  d = colSums(M) / (ncol(M) - 1)
  d_zs = (d - mean(d))/sd(d)
  c_sp = which(d_zs > 1)
  
  # change interaction outcomes of central species only
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      if(M[i,j] == 1 & any(c_sp == i)){
        M[j,i] = 0
        V[j,i] = 1
      }
    }
  }
  
  # no intraespecific interactions
  diag(V) = 0
  
  # estimate p as a frequency
  A = M + V
  p_vl = sum(V) / sum(A)
  
  # create and return a list with the positive and negative matrices\
  mats = list(M, V, p_vl)
  return(mats)
}
#---------------------------------------------------------------------------------------------------#