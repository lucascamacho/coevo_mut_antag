#---------------------------------------------------------------------------------------#
MeanPairDist = function(z_mat){
  # Calculate the absolute difference of traits between all species.
  # You can also find this metric called Mean Pairwise Distance (Ciamplaglio, 2001)
  # 
  #Args:
  #  z_mat: matrix of species traits. Rows are timesteps and columns are species.
  # 
  # Return:
  #  difs: vector containing the absolute difference of traits.
  
  # difference between all species
  difs = sapply(z_mat,"-", z_mat)
  
  # sum of the absolute values of trait differences
  dif_total = sqrt(sum(difs ** 2))
  
  # calculate average and return result
  mpd = dif_total / (length(z_mat) * (length(z_mat) - 1))
  return(mpd)
}
#---------------------------------------------------------------------------------------#