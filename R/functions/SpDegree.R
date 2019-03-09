#-----------------------------------------------------------------------------------------------------#
SpDegree = function(M, V){
  # Count the number of AA, AM and MM interactions for all the species in the network
  # Args:
  #  M: Mutualistic network of interactions
  #  V: Antagonistic network of interactions
  # Return:
  #  degree: matrix with species in columns and the degree fot each interaction type in rows
  #
  # create a matrix with species in columns and 3 rows for degree of each interaction type (AA, AM and MM)
  #
  degree = matrix(NA, ncol = ncol(M), nrow = 3) # create matrix to insertthe results
  rownames(degree) = c("AA", "AM", "MM")
  colnames(degree) = c(seq(1, ncol(M), 1))
  
  # using M and V to find the interaction classes
  v = (V == 1) == (t(V) == 1) # antagonism
  v[V == 0] = FALSE # ignore the zero's
  
  c = (V == 1) != (t(V) == 1) # cheaters
  
  m = (M == 1) == (t(M) == 1) # mutualism
  m[M == 0] = FALSE # ignore the zero's

  #fill the degree matrix with the sum of interactions of certain class
  degree[1,] = apply(v, 2, sum)
  degree[2,] = apply(c, 2, sum)
  degree[3,] = apply(m, 2, sum)
  
  # return the matrix with degree of species separating by interaction class
  return(degree)
}
#-----------------------------------------------------------------------------------------------------#