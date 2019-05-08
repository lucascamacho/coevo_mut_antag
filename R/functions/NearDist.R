#------------------------------------------------------------------------------------------#
NearDist = function(z_mat){
  # Function to calculate the Mean Nearest and Distant Neighbor Distance metric.
  # This metric is calculated to quantify the trait diversification of species
  # after the coevoluionary process. For that we sum the nearest and distant
  # pair of species traits controlling by species numbers
  #
  # Args:
  #   z_matrix: matriz with species in collums and traits varying in time (rows)
  #
  # Return:
  #   d = list with the nearest and distant neighbor distance metric values.
  
  # get all the distances between species
  distances = sapply(z_mat, "-" , z_mat)
  
  # sum the lowest and highest distancies
  min_dist = apply(abs(distances), 1, FUN = function(x) {min(x[x > 0])})
  max_dist = apply(abs(distances), 1, FUN = function(x) {max(x[x > 0])})
  
  # sum of the distances and control by species numbers
  min_n = sum(min_dist) / ncol(z_mat)
  max_n = sum(max_dist) / ncol(z_mat)
  
  # return results
  d = list(min_n, max_n)
  return(d)
}
#------------------------------------------------------------------------------------------#