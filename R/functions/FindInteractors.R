#-----------------------------------------------------------------------------------------------------#
FindInteractors = function(M, V){
  # Identify which species have double-positive, positive-negative or double-negative effects.
  #
  # Args:
  #   M: matrix of positive effects
  #   V: matrix of negative effects
  #
  # Return:
  #  index: data frame with species indexation separate by type of effects combinations
  
  # AA outcomes
  v = (V == 1) == (t(V) == 1)
  v[V == 0] = FALSE
  aa_posit = which(v == TRUE, arr.ind = TRUE)
  aa_index = unique(aa_posit[,1])# get the position of AA
  
  # AM outcomes
  am_posit = which(M != t(M), arr.ind = TRUE) 
  am_index =c(unique(am_posit[,1])) # get the position of AM
  
  # MM outcomes
  m = (M == 1) == (t(M) == 1)
  m[M == 0] = FALSE
  mm_posit = which(m == TRUE, arr.ind = TRUE) 
  mm_index = unique(mm_posit[,1])# get the position of MM
  
  # create, organize and return a list of species AA, AM and MM
  index = list(aa_index, am_index, mm_index)
  index = lapply(index, sort)
  return(index)
}
#-----------------------------------------------------------------------------------------------------#