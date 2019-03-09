#-----------------------------------------------------------------------------------------------------#
FindInteractors = function(M, V){
  # Identify which species are Mutualists MM, Antagonists AA and Cheaters AM in M and V matrices.
  #
  # Args:
  # M: binary adjacency matrix of mutualistic interactions
  # V: binary adjacency matrix of antagonistic interactions
  #
  # Return:
  # single data.frame with index of mutualists, antagonists and cheaters 
  # in the adjancency matrix M or V.
  #
  # antagonism
  # comparing the adjacency matrix V with his Transpose
  v = (V == 1) == (t(V) == 1)
  v[V == 0] = FALSE
  aa_posit = which(v == TRUE, arr.ind = TRUE) # get the position of antagonists
  aa_index = unique(aa_posit[,1])
  
  # cheaters
  # comparing the adjacency matrix V with the difference of his Transpose 
  am_posit = which(M != t(M), arr.ind = TRUE) # get the position of cheaters
  am_index =c(unique(am_posit[,1])) 
  
  # mutualism
  # comparing the adjacency matrix M with his Transpose
  m = (M == 1) == (t(M) == 1)
  m[M == 0] = FALSE
  mm_posit = which(m == TRUE, arr.ind = TRUE) # get the position of mutualists
  mm_index = unique(mm_posit[,1])
  
  # create, organize and return a list of species AA, AM and MM
  index = list(aa_index, am_index, mm_index)
  index = lapply(index, sort)
  return(index)
}
#-----------------------------------------------------------------------------------------------------#