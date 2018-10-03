EmpZeroLines = function(M, V){
  
  source("SquareMatrix.R")
  
  for(z in 1:1000000){
    A = M + V
    r = apply(A, 1, sum)
    c = apply(A, 1, sum)
    
    if(any(r == 0) & any(c == 0)){
     row_zeros = which(r == 0)
     col_zeros = which(c == 0)
     M = M[-row_zeros, -col_zeros]
     V = V[-row_zeros, -col_zeros]
     warning(paste("Species", row_zeros, "have been removed from the matrix. "))

    }
    else{
      break
    }
  }
  
  M = SquareMatrix(M)
  V = SquareMatrix(V)
  mats = list(M, V)
  return(mats)

}