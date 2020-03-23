#------------------------------------------------------------------------------------------#
Disconnected = function(M, V){
  # Identify and count which species are disconnected of the network based on the M and V.
  # If the line and collum are zero, the species are exclude of the network.
  # Args:
  #   M: matrix of positive effects
  #   V: matrix of negative effects
  # Return:
  #   res: list with which species are disconnected and how many are disconnected
  
  # get which species has rows and columms with 0's
  A = M + V
  r = which(rowSums(A) == 0)
  c = which(colSums(A) == 0)
  
  # get the overlap between r and c
  index = r[r %in% c]
  
  # the number of overlaps are the number of species disconnected in the matrix
  num = length(index)
  
  # index are which species are disconnnected
  # num is the number of species disconnected
  res = list(index, num)
  return(res)
}
#------------------------------------------------------------------------------------------#