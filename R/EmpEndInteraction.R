EmpEndInteraction = function(M, V, interaction){
  # Choose a type of interaction to ignore in your empirical adjacency matrices
  # The interactions are identified based on the index of species which has a 1
  # in V and the Transpose of V (in case of antagonism). 
  # These index are used to insert 0 on the matrices.
  #
  #Args:
  # M: empirical mutualistic adjacency matrix
  # V: empirical antagonistic adjacency matrix
  # interaction: type on interaction you want to ignore (mutualism MM, antagonism AA or cheaters AM)
  #
  #Return:
  # list of matrices M and V with 0 in interactions that you choose
  #
  if(interaction != "mutualism" & interaction != "antagonism" & interaction != "cheaters"){
    stop("Choose a valid interaction type")
  }
  
  if(interaction == "antagonism"){ # end antagonism based on index of V and V transpose
    index = V == t(V)
    index[V == 0] = FALSE
    V[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  if(interaction == "cheaters"){ # end cheaters based on index of V and diff V transpose
    index = V != t(V)
    V[index] = 0
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
  
  else{
    index = M == t(M) # end mutualism based on index of M and M transpose
    index[M == 0] = FALSE
    M[index] = 0
    mats = list(M, V)
    return(mats)
  }
}