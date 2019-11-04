GiantComponent = function(A){
  # This function provides a measure of giant component in ecological networks
  # apllying the condition fow undirected graphs:
  # E[kˆ2] - 2E[k] > 0 the giant component exists
  # E[kˆ2] - 2E[k] < 0 the giant component don't exist
  #
  # Args:
  #  A = adjancency matrix of interactions
  # Return:
  #  the value of the condition applied 
  
  # get the full adjancency matrix
  V = as.matrix(V)
  d_1 = mean(rowSums(A) ** 2) # E[kˆ2]
  d_2 = 2 * mean(rowSums(A)) # 2E[k]
  
  res_degree = d_1 - d_2 # E[kˆ2] - 2E[k]
  return(res_degree) # return the condition value
}