#-----------------------------------------------------------------------------------------------------#
BalanDiver = function(z_mat){
  # Calculate the metrics of how much the species traits fluctuates in time
  # and calculate the variance between traits for each timestep in simulation.
  #
  # Args:
  #  z_mat: matrix with all species traits in time
  # Return:
  #  baladiver: list with 2 vectors. A fluctuations calc vector and a vars vector.
  #
  # calculate the species fluctuations in time
  balan = abs(z_mat[nrow(z_mat), ] - z_mat[1, ]) / sum(abs(apply(z_mat, 2, 
                                                                 function(x) x[-1] - x[-length(x)])))
  # traits variance for each timestep
  diver = apply(z_mat, 1, var)
  
  balandiver = list(balan, diver) # define a list of vectors and return the list
  return(balandiver)
}
#-----------------------------------------------------------------------------------------------------#