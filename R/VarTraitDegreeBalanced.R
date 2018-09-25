VarTraitDegreeBalanced = function(z_mat){
  # Calculate the absolutte difference of traits between species by interaction type 
  # Args:
  #  z_mat: matrix of species traits. Rows are time steps and columns are species
  # 
  # Return:
  # vector containing the absolutte difference of traits by interaction type (AM and MM)

  # check if there are cheaters in the network
  # if there are no cheaters, insert NA and calculate only for mutualists
  if(sum(index[[2]]) == 0){
    dif_c = NA
    
    # mutualisms calculation
    dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
    dif_mut = abs(dif_mut) # module of differences
    dif_m = sum(dif_mut[lower.tri(dif_mut)]) # sum of absolute difference of mutualists traits
    
    # return the differences
    difs = c(dif_c, dif_m)
    return(difs)
  }
  
  # check if there are mutualists in the network
  # if there are no mutualists, insert NA and calculate only for cheaters
  if(sum(index[[3]]) == 0){
    
    dif_m = NA
    
    # cheaters calculation
    dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
    dif_cheat = abs(dif_cheat) # module of differences
    dif_c = sum(dif_cheat[lower.tri(dif_cheat)]) # sum of absolute difference of cheaters traits
    
    # return the differences
    difs = c(dif_c, dif_m)
    return(difs)
  }
  
  # in case there are mutualists and cheaters
  # calculate the difference for then
  else{
    # cheaters calculation
    dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
    dif_cheat = abs(dif_cheat)
    dif_c = sum(dif_cheat[lower.tri(dif_cheat)])
    
    # mutualisms calculation
    dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
    dif_mut = abs(dif_mut)
    dif_m = sum(dif_mut[lower.tri(dif_mut)])   
 
    #return the differences
    difs = c(dif_c, dif_m)
    return(difs)
  }
}