#----------------------------------------------------------------------------------------------#
Counting = function(M, V){
  # Count the number of double-negative (AA), negative-positive (AM) 
  # and double-positive (MM) effects between species.
  #
  # Args:
  #   M: matrix of positive effects
  #   V: matrix of negative effects
  #
  # Return:
  #   A list with the number of AA, AM and MM interactions
  n_cheat = sum(V) #counting AM

  m1 = which(M == 1)
  m2 = which(t(M) == 1)
  n_mut = sum(table(m1[m1 %in% m2])) / 2 #counting MM

  v1 = which(V == 1)
  v2 = which(t(V) == 1)
  n_inter = sum(table(v1[v1 %in% v2])) / 2 # counting AA

  numbers = list(n_inter, n_cheat, n_mut) #create and return a list with numbers of AA, AM and MM
  return(numbers)

}
#----------------------------------------------------------------------------------------------#