#-----------------------------------------------------------------------------------------------------#
VarTraitDegreeBalanced = function(z_mat){
  # Calculate the absolute difference of traits between species by interaction outcome type 
  # Args:
  #  z_mat: matrix of species traits. Rows are timesteps and columns are species.
  # 
  # Return:
  #  difs: vector containing the absolute difference of traits by interaction type (AM and MM).
  #
  #  check if there are antagonisms in the network
  if(sum(index[[2]]) == 0 & sum(index[[3]]) == 0){
  # no interactions in positive and negative outcomes
    dif_c = NA
    dif_m = NA
    
    # return the differences
    difs = c(dif_c, dif_m)
    return(difs)
    }
  else{
    # there are both positive and negative outcomes
    if(sum(index[[2]]) != 0 & sum(index[[3]]) != 0){
      # antagonism calculation
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
      # there are only negative outcomes
      if(sum(index[[3]]) == 0){
        
        dif_m = NA
        
        # negative outcomes calculation
        dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
        dif_cheat = abs(dif_cheat)
        dif_c = sum(dif_cheat[lower.tri(dif_cheat)])
        dif_c = dif_c / c[[2]]
        
        # return the differences
        difs = c(dif_c, dif_m)
        return(difs)
      }
      else{
      # there are only mutualisms interactions
        dif_c = NA
        
        # positive outcomes calculation
        dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
        dif_mut = abs(dif_mut)
        dif_m = sum(dif_mut[lower.tri(dif_mut)])
        dif_m = dif_m / c[[3]]
        
        # return the differences
        difs = c(dif_c, dif_m)
        return(difs)
      }
    }
  }
}
#-----------------------------------------------------------------------------------------------------#