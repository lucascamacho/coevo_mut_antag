SquareMatrix <- function(rec.mat){
  # Mirrors the rectangular matrix to transform it into a square matrix
  #
  #Args:
  #      rec.mat = Rectangular adjacency matrix
  #Return:
  #      cor.sq = Square adjacency matrix
  #
  rownames(rec.mat) <- paste("A", 1:nrow(rec.mat), sep="")  
  sq.1 <- matrix(0, nrow(rec.mat), nrow(rec.mat))
  sq.1 <- rbind(sq.1, t(rec.mat))
  dimnames(sq.1)[[1]][1:nrow(rec.mat)] <- colnames(sq.1)
  sq.2 <- matrix(0, ncol(rec.mat), ncol(rec.mat))
  dimnames(sq.2) <- list(colnames(rec.mat), colnames(rec.mat))
  sq.2 <- rbind(rec.mat, sq.2)
  
  cor.sq <- cbind(sq.1, sq.2)
  
  return(cor.sq)
  
}
