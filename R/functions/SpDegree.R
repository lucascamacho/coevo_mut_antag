#-----------------------------------------------------------------------------------------------------#
SpDegree = function(M, V){
  # Count the number of AA, AM and MM effects for all the species in the matrix
  #
  # Args:
  #   M: matrix of positive effects
  #   V: matrix of negative effects
  #
  # Return:
  #  degree: matrix with species in rows and number of AA, AM and MM in collums.
  degree = matrix(NA, ncol = ncol(M), nrow = 3) # create matrix to insert results
  rownames(degree) = c("AA", "AM", "MM")
  colnames(degree) = c(seq(1, ncol(M), 1))
  
  # using M and V to find the effects classes
  v = (V == 1) == (t(V) == 1) # AA
  v[V == 0] = FALSE # ignore the zero's
  
  c = (V == 1) != (t(V) == 1) # AM
  
  m = (M == 1) == (t(M) == 1) # MM
  m[M == 0] = FALSE # ignore the zero's

  #fill the degree matrix with the sum of signs for classes
  degree[1,] = apply(v, 2, sum)
  degree[2,] = apply(c, 2, sum)
  degree[3,] = apply(m, 2, sum)
  
  # return the matrix with degree of species separating by effects
  return(t(degree))
}
#-----------------------------------------------------------------------------------------------------#