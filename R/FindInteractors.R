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
  v = (V == 1) == (t(V) == 1)
  v[V == 0] = FALSE
  aa_posit = which(v == TRUE, arr.ind = TRUE)
  aa_index = unique(aa_posit[,1])
  
  # comparing the adjacency matrix with his Transpose # 1, 2, 5, 3, 
  am_posit = which(M != t(M), arr.ind = TRUE) # the difference between M and m transpose are the cheaters (we can use V too)
  am_index =c(unique(am_posit[,1])) # get the position of these cheaters to use as an index
  
  # comparing the adjacency matrix with his Transpose
  m = (M == 1) == (t(M) == 1)
  m[M == 0] = FALSE
  mm_posit = which(m == TRUE, arr.ind = TRUE)
  mm_index = unique(mm_posit[,1])
  
  index = list(aa_index, am_index, mm_index)
  index = lapply(index, sort)
  return(index)
  
}