#-----------------------------------------------------------------------------------------------------#
EndInteraction = function(M, V, interaction){
  # Choose a type of outcome combination to ignore in your matrices. 
  #  The combination of outcomes choosed will be turned to 0.
  #
  #Args:
  #   M: matrix of positive interaction outcomes
  #   V: matrix of negative interaction outcomes
  #   interaction: type of combination of outcomes you want to ignore
  #
  #Return:
  #  list of matrices M and V with 0 in interactions MM, AA or AM
  
  # check the choosed parameter
  if(interaction != "mutualism" & interaction != "antagonism" & interaction != "interference"){
    stop("Choose between interference, antagonism and mutualism in interaction parameter")
  }
  
  # end outcomes combinations based on "interaction" parameter
  if(interaction == "interference"){ # end AA
    index = V == t(V)
    V[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  if(interaction == "antagonism"){ # end AM
    index = V != t(V)
    V[index] = 0
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  if(interaction == "mutualism"){ # end MM
    index = M == t(M)
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
}
#-----------------------------------------------------------------------------------------------------#