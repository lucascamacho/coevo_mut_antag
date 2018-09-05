FindInteractors = function(M, V){
  # Identify which species are Mutualists, Antagonists and Cheaters in two adjacency matrix of interactions
  #
  # Args:
  # M: binary adjacency matrix of mutualistic interactions
  # V: binary adjacency matrix of antagonistic interactions
  #
  # Return:
  # single data.frame with index of mutualists, antagonists and cheaters in the adjancecy matrix M or V
  #
  # comparing the adjacency matrix with his Transpose
  v1 = which(V == 1, arr.ind = TRUE)
  v2 = which(t(V) == 1, arr.ind = TRUE)
  aa_index = c(unique(v1[v1 %in% v2]))

  # comparing the adjacency matrix with his Transpose
  am_posit = which(M != t(M), arr.ind = TRUE) # the difference between M and m transpose are the cheaters (we can use V too)
  am_index =c(unique(am_posit[,1])) # get the position of these cheaters to use as an index
  
  # comparing the adjacency matrix with his Transpose
  m1 = which(M == 1, arr.ind = TRUE)
  m2 = which(t(M) == 1, arr.ind = TRUE)
  mm_index = c(unique(m1[m1 %in% m2]))
  
  index = list(aa_index, am_index, mm_index)
  return(index)
  
}