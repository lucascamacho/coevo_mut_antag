TraitDegreeBalanced = function(z_mat){
  #
  #
  
  # create empty vectors
  to_sum_c = c()
  to_sum_m = c()

  # loop for eanch specie
  for(g in 1:length(z_mat)){
    
    #cheaters
    trait_k_c = degree[2, g] * z_mat[g]
    to_sum_c = append(to_sum_c, trait_k_c)
    
    # mutualism
    trait_k_m = degree[3, 2] * z_mat[2]
    to_sum_m = append(to_sum_m, trait_k_m)
  }
  
  if(sum(degree[2,]) == 0){
    t_am = NA
    s_m = sum(to_sum_m) / sum(degree[3,])
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat)
    
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  if(sum(degree[3,]) == 0){
    t_mm = NA
    s_c = sum(to_sum_c) / sum(degree[2,])
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat)
    
    t_final = c(t_am, t_mm)
    return(t_final)
  }
  
  else{
    s_c = sum(to_sum_c) / sum(degree[2,])
    t_am = abs((s_c -  mean(z_mat))) / sd(z_mat)
  
    s_m = sum(to_sum_m) / sum(degree[3,])
    t_mm = abs((s_m - mean(z_mat))) / sd(z_mat)
  
    t_final = c(t_am, t_mm)
    return(t_final)
  }
}