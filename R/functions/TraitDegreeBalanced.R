#-----------------------------------------------------------------------------------------------------#
TraitDegreeBalanced = function(z_mat){
  # Calculate the z-score of mean trait value for cheaters and mutualists balanced by degree
  # 
  #Args:  
  #  z_mat: matriz of species traits after simulations, rows are time steps and columns are species
  #
  #Return:
  #  vector with z-score of mean trait value for cheaters and mutualists species
  #
  # create empty vectors
  to_sum_c = c()
  to_sum_m = c()

  # loop for each specie
  for(g in 1:length(z_mat)){
    
    #cheaters
    trait_k_c = degree[2, g] * z_mat[g] # trait * degree of cheaters
    to_sum_c = append(to_sum_c, trait_k_c) # insert these values in vector
    
    # mutualism
    trait_k_m = degree[3, g] * z_mat[g] # trait * degree of mutualisms
    to_sum_m = append(to_sum_m, trait_k_m) # insert hese values in vector
  }
  
  # check if there are no cheaters in the network
  # insert NA in cheaters z-score and calculate only for mutualists
  if(sum(degree[2,]) == 0){
    t_am = NA
    s_m = sum(to_sum_m) / sum(degree[3,])
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat) # z-score calculation
    
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  # check if there are no mutualists in the network
  # insert NA in mutualists z-score and calculate only for cheaters
  if(sum(degree[3,]) == 0){
    t_mm = NA
    s_c = sum(to_sum_c) / sum(degree[2,])
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat) # z-score calculation
    
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  # in case are mutualists and cheaters in the network
  # calculate z-score
  else{
    s_c = sum(to_sum_c) / sum(degree[2,]) # cheaters z-score
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat)
  
    s_m = sum(to_sum_m) / sum(degree[3,]) # mutualists z-score
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat)
  
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
}
#-----------------------------------------------------------------------------------------------------#