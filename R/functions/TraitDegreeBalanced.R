#-------------------------------------------------------------------------------------------------------#
TraitDegreeBalanced = function(z_mat){
  # Calculate the z-score of mean trait value for species involved in positive or negative outcomes
  # of interactions.
  # 
  #Args:  
  #  z_mat: matriz of species traits after simulations, rows are timesteps and columms are species
  #
  #Return:
  #  t_final: vector with z-score of mean trait value for species involved in MM and AM interactions.
  #
  to_sum_c = vector()
  to_sum_m = vector()

  for(g in 1:length(z_mat)){ # loop to multiply the species trait by their degree (Kmm or Kam)
    # antagonisms
    trait_k_c = degree[g, 2] * z_mat[g] # trait * degree of antagonists
    to_sum_c = append(to_sum_c, trait_k_c) # insert these values in vector
    
    # mutualisms
    trait_k_m = degree[g, 3] * z_mat[g] # trait * degree of mutualisms
    to_sum_m = append(to_sum_m, trait_k_m) # insert hese values in vector
  }
  
  # check if there are no antagonists in the network
  # insert NA in antagonists z-score and calculate only for mutualists
  if(sum(degree[,2]) == 0){
    t_am = NA
    s_m = sum(to_sum_m) / sum(degree[,3])
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat) # z-score calculation
    
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  # check if there are no mutualists in the network
  # insert NA in mutualists z-score and calculate only for antagonists
  if(sum(degree[,3]) == 0){
    t_mm = NA
    s_c = sum(to_sum_c) / sum(degree[,2])
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat) # z-score calculation
    
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  # in case are mutualists and antagonists in the network
  # calculate z-score
  else{
    s_c = sum(to_sum_c) / sum(degree[,2]) # antagonists z-score
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat)
  
    s_m = sum(to_sum_m) / sum(degree[,3]) # mutualists z-score
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat)
  
    # return z-scores
    t_final = c(t_am, t_mm)
    return(t_final)
  }
}
#-----------------------------------------------------------------------------------------------------#