#-----------------------------------------------------------------------------------------------------#
EmpEndInteraction = function(M, V, interaction){
  # Choose a type of outcome to ignore in your empirical adjacency matrices
  #
  #Args:
  # M: empirical matrix of mutualism
  # V: matrix of negative effects
  # interaction: type on interaction outcome you want to ignore
  #  You can choose between Interferences, Antagonisms and Mutualisms
  #
  #Return:
  # list of matrices M and V with 0 in effects that you choose
  #
  if(interaction != "mutualism" & interaction != "interference" & interaction != "antagonism"){
    stop("Choose a valid interaction type: mutualism, antagonism or interference")
  }
  
  if(interaction == "interference"){ # end interferences
    index = V == t(V)
    index[V == 0] = FALSE
    V[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  if(interaction == "antagonism"){ # end antagonisms
    index = V != t(V)
    V[index] = 0
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  if(interaction == "mutualism"){ # end mutualism
    index = M == t(M)
    index[M == 0] = FALSE
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
}
#-----------------------------------------------------------------------------------------------------#