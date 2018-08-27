MeanTraitInteract = function(M, V, z_mat){
  # Calculate the mean trait values for certain types of interactions (antagonism, cheaters and mutualism)
  #
  # Args:
  # M: mutualistic binary matrix of interactions
  # V: antagonistic binary matrix of interactions
  # z_mat: matrix with trait values of species. 
  # 
  # Obs: z_mat must have species in collums and trait values in rows
  #
  # Return:
  # vector containing the mean trait value for MA, AM and MM interactions
  #
  
  # mean trait of antagonists
  v1 = which(V == 1) 
  v2 = which(t(V) == 1)
  aa_index = c(intersect(v1, v2)) # the intersect of V and transpose V are the pure antagonists
  aa_diff = mean(abs(sapply(z_mat[aa_index], "-", z_mat[aa_index]))) # calculate the mean difference between these species
  
  # mean trait of cheaters
  am_posit = which(M != t(M), arr.ind = TRUE) # the difference between M and m transpose are the cheaters (we can use V too)
  am_index = c(unique(am_posit[,1])) # get the position of these cheaters to use as an index
  am_diff = mean(abs(sapply(z_mat[am_index], "-", z_mat[am_index]))) # calculate the mean difference between these species
  
  # mean trait of mutualists
  m1 = which(M == 1)
  m2 = which(t(M) == 1)
  mm_index = c(intersect(m1, m2)) # the intersect of M and transpose M are the pure mutualists
  mm_diff = mean(abs(sapply(z_mat[mm_index], "-", z_mat[mm_index]))) # calculate the mean difference between these species 

  diff_interactions = vector(aa_diff, am_diff, mm_diff) # single vector with mean differences for each type of interaction
  
  return(diff_interactions) # return the final vector

}