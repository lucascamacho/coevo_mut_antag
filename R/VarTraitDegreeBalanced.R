VarTraitDegreeBalanced = function(z_mat){
  if(sum(index[[2]]) == 0){
    dif_c = NA
    
    #mutualisms
    dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
    dif_mut = abs(dif_mut)
    dif_m = sum(dif_mut[lower.tri(dif_mut)])   
    
    difs = c(dif_c, dif_m)
    return(difs)
  }
  
  if(sum(index[[3]]) == 0){
    
    dif_m = NA
    
    #cheaters
    dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
    dif_cheat = abs(dif_cheat)
    dif_c = sum(dif_cheat[lower.tri(dif_cheat)])
    
    difs = c(dif_c, dif_m)
    return(difs)
  }
  
  else{
    #cheaters
    dif_cheat = sapply(z_mat[index[[2]]],"-", z_mat[index[[2]]])
    dif_cheat = abs(dif_cheat)
    dif_c = sum(dif_cheat[lower.tri(dif_cheat)])
    
    #mutualisms
    dif_mut = sapply(z_mat[index[[3]]],"-", z_mat[index[[3]]])
    dif_mut = abs(dif_mut)
    dif_m = sum(dif_mut[lower.tri(dif_mut)])   
 
    difs = c(dif_c, dif_m)
    return(difs)
  }
}