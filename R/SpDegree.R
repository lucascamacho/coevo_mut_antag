SpDegree = function(M, V){
  # Count the number of AA, AM and MM interactions for all the species in the network
  # Args:
  #    M: Mutualistic network of interactions
  #    V: Antagonistic network of interactions
  # Return:
  #    matrix with species in collums and the degree fot each interaction type in rows
  #
  
  # create a matrix with species in columms and 3 rows for degree of each interaction type
  degree = matrix(NA, ncol = ncol(M), nrow = 3)
  rownames(degree) = c("AA", "AM", "MM")
  
  # using M and V to find the interaction classes

which(V == 1)
which(t(V) == 1)
  
  c = (V == 1) != (t(V) == 1) # cheaters
  m = (M == 1) == (t(M) == 1) # mutualism

  #fill the degree matrix with the sum of interactions of certain class
  degree[1,] = apply(v, 2, sum)
  degree[2,] = apply(c, 2, sum)
  degree[3,] = apply(m, 2, sum)
  
  # transform the matrix in data frame and return
  degree = data.frame(degree)
  return(degree)
  
}
which(V == t(V))
V
