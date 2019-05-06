#------------------------------------------------------------------------------------------#
PartRatio = function(z_mat){
  # Function to calculate the Participation Ratio metric of species traits
  # using the min and max range of traits. This function calculates this metric
  # by quantifying how much of the morphospace are used. The participation ratio range
  # from 1 (species are all simillar) to n_sp, which is the number of species (species 
  # are equally different).
  #  Args:
  #   z_matrix: matriz with species in collums and traits varying in time (rows)
  #  Return:
  #   part_ratio: participation ratio metric (ranging from 1 to n_sp) 
  
  # get the morphospace range, and create the bins for each specie
  max_range = max(z_mat)
  min_range = min(z_mat)
  create_bins = cut(z_mat, ncol(z_mat))
  
  # get the frequency of occupied bins
  freq_bins = table(create_bins)
  freq_bins = freq_bins / ncol(z_mat)
  
  # sum of the squares of all occupied bins values
  bins_sums = sum((freq_bins ** 2))
  
  # calculate participation ratio and return this values
  part_ratio = 1 / bins_sums
  return(part_ratio)
}
#------------------------------------------------------------------------------------------#
