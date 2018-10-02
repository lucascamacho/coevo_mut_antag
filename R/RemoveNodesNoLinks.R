#-----------------------------------------------------------------------------------------------------#

RemoveNodesNoLinks = function(mat) {
  # Removes nodes with no links (degree = 0) from a network.
  #
  # Args:
  #   mat: bipartite adjacency matrix
  #           
  # Returns:
  #   A new matrix without the row(s) and/or column(s) corresponding to nodes with no links.
  #   A warning is printed indicating which nodes have been removed.

  k_row = apply(mat, 1, sum) # degree of row nodes
  k_col = apply(mat, 2, sum) # degree of columns nodes
  if (any(k_row == 0)) { 
    k_row_zeros = which(k_row == 0) # node positions
    mat = mat[-k_row_zeros, ] # remove nodes
    warning(paste("Row(s)", k_row_zeros, "have been removed from the matrix. "))
  }
  if (any(k_col == 0)) { 
    k_col_zeros = which(k_col == 0) # node positions
    mat = mat[ ,-k_col_zeros] # remove nodes
    warning(paste("Column(s)", k_col_zeros, "have been removed from the matrix. "))
  }
  return(mat)
}

#-----------------------------------------------------------------------------------------------------#