#---------------------------------------------------------------------------------------#
MeanPairDist = function(z_mat){
  # Calculate the absolute difference of traits between all species.
  # You can also find this metric called Mean Pairwise Distance (Ciamplaglio, 2001)
  #
  # Ciampaglio, Charles N., Matthieu Kemp, e Daniel W. McShea. 
  # “Detecting changes in morphospace occupation patterns in the fossil record: 
  # characterization and analysis of measures of disparity”. 
  # Paleobiology 27, nº 4 (dezembro de 2001): 695–715. 
  # https://doi.org/10.1666/0094-8373(2001)027<0695:DCIMOP>2.0.CO;2.
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