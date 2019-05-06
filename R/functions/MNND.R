#------------------------------------------------------------------------------------------#
MNND = function(z_mat){
  # Function to calculate the Mean Nearest Neighbor Distance metric.
  # This metric is calculated to quantify the trait diversification of species
  # after the coevoluionary process. 
  # Args:
  #   z_matrix: matriz with species in collums and traits varying in time (rows)
  # Return:
  #   mnnd = mean nearest neighbor distance value
  
  # get all the distances between species
  distances = sapply(z_mat, "-" , z_mat)
  # sum the lowest distancies
  min_dist = apply(abs(distances), 1, FUN = function(x) {min(x[x > 0])})
  # normalize by species numbers
  d = sum(min_dist) / ncol(z_mat)
  return(d)
}
#------------------------------------------------------------------------------------------#