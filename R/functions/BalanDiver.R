#-----------------------------------------------------------------------------------------------------#
BalanDiver = function(z_mat){
  # Calculate the directionality of species traits, a measure of how much traits fluctuate in time.
  # Also, calculates the variance of species traits for each timestep of simulation.
  #
  # Args:
  #   z_mat: matrix with species traits in time
  # Return:
  #   baladiver: list with 2 vectors. A fluctuations vector and a vars vector.
  #
  # calculate the species fluctuations in time
  balan = abs(z_mat[nrow(z_mat), ] - z_mat[1, ]) / apply(abs(diff(z_mat)), 2, sum)
  
  # traits variance for each timestep
  diver = apply(z_mat, 1, var)
  
  balandiver = list(balan, diver) # define a list of vectors and return the list
  return(balandiver)
}
#-----------------------------------------------------------------------------------------------------#