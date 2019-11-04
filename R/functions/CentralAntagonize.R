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
  
  # identify the central species using degree centrality
  d = colSums(M)
  
  # identify the number of species in both sides of bip. matrix
  nomes = colnames(M)
  animals = nomes[grepl("^A.*", nomes)]
  plants = nomes[grepl("^P.*", nomes)]
  
  # calculate degree centrality in bipartite matrix
  for(l in 1:length(d)){
    if(any(names(d)[l] == plants) == TRUE){
      d[l] = d[l] / length(animals)
    }
    else{
      d[l] = d[l] / length(plants)
    }
  }
  
  #degree centrality z-score
  d_zs = (d - mean(d)) / sd(d)
  # choosing species that has > 1 z-score
  c_sp = which(d_zs > 1)
  
  # read only the matrix upper.tri
  index_pos = which(lower.tri(M))
  index_sp = which(lower.tri(M), arr.ind = TRUE)
  
  # change interaction outcomes of central species only
  for(i in 1:length(index_pos)){ # for each element in index
    if(M[index_pos[i]] == 1 & any(c_sp == index_sp[i])){
      M[index_sp[i,][2], index_sp[i,][1]] = 0 # inverse element is zero (M[j,i])
      V[index_sp[i,][2], index_sp[i,][1]] = 1 # inverse element in V is "on" (V[j,i])
    }
  }
  
  # read only the matrix upper.tri
  index_pos = which(upper.tri(M))
  index_sp = which(upper.tri(M), arr.ind = TRUE)
  
  # change interaction outcomes of central species only
  for(i in 1:length(index_pos)){ # for each element in index
    if(M[index_pos[i]] == 1 & any(c_sp == index_sp[i])){
      M[index_sp[i,][2], index_sp[i,][1]] = 0 # inverse element is zero (M[j,i])
      V[index_sp[i,][2], index_sp[i,][1]] = 1 # inverse element in V is "on" (V[j,i])
    }
  }

  # no intraespecific interactions
  diag(V) = 0
  
  # estimate p as a frequency
  A = M + V
  p_vl = sum(V) / (sum(A) / 2)
  
  # create and return a list with the positive and negative matrices\
  mats = list(M, V, p_vl)
  return(mats)
}
#---------------------------------------------------------------------------------------------------#