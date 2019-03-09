VarTraitDegreeBalanced = function(z_mat){
  # Calculate the absolutte difference of traits between species by interaction type 
  # Args:
  #  z_mat: matrix of species traits. Rows are time steps and columns are species
  # 
  # Return:
  # vector containing the absolutte difference of traits by interaction type (AM and MM)

  # check if there are cheaters in the network
  # if there are no cheaters, insert NA and calculate only for mutualists
  if(sum(index[[2]]) == 0 & sum(index[[3]]) == 0){
  # no interactions in cheaters and mutualisms
    dif_c = NA
    dif_m = NA
    
    # return the differences
    difs = c(dif_c, dif_m)
    return(difs)
    }
  else{
    # there are both interactions of cheaters and mutulisms
    if(sum(index[[2]]) != 0 & sum(index[[3]]) != 0){
      # cheaters calculation
      dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
      dif_cheat = abs(dif_cheat)
      dif_c = sum(dif_cheat[lower.tri(dif_cheat)])
      dif_c = dif_c / c[[2]]
      
      # mutualisms calculation
      dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
      dif_mut = abs(dif_mut)
      dif_m = sum(dif_mut[lower.tri(dif_mut)])   
      dif_m = dif_m / c[[3]]
      
      #return the differences
      difs = c(dif_c, dif_m)
      return(difs)
    }
    else{
      # there are only cheaters interactions
      if(sum(index[[3]]) == 0){
        
        dif_m = NA
        
        # cheaters calculation
        dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
        dif_cheat = abs(dif_cheat) # module of differences
        dif_c = sum(dif_cheat[lower.tri(dif_cheat)]) # sum of absolute difference of cheaters traits
        dif_c = dif_c / c[[2]]
        
        # return the differences
        difs = c(dif_c, dif_m)
        return(difs)
      }
      else{
      # there are only mutualisms interactions
        dif_c = NA
        
        # mutualisms calculation
        dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
        dif_mut = abs(dif_mut) # module of differences
        dif_m = sum(dif_mut[lower.tri(dif_mut)]) # sum of absolute difference of mutualists traits
        dif_m = dif_m / c[[3]]
        
        # return the differences
        difs = c(dif_c, dif_m)
        return(difs)
      }
    }
  }
}  