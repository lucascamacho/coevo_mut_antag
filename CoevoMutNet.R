#-----------------------------------------------------------------------------------------------------#

CoevoMutNet = function(n_sp, f, h, alpha, theta, init, m, epsilon, t_max) {
  # Simulates the coevolutionary dynamics on a mutualistic network of species interactions
  #
  # Args:
  #   n_sp: total number of species in the mutualistic network
  #   f: square adjacency matrix representing the mutualistic network
  #   h: vector of heritability values
  #   alpha: parameter alpha, sensitivity of selection to trait matching
  #   theta: vector of environmental optimum values
  #   init: vector of initial trait values
  #   m: vector of proportion of selection due to mutualism
  #   epsilon: value to determine when equilibrium is reached
  #   t_max: maximum number of timesteps allowed
  #   
  # Obs: 
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   A matrix containing, in each row t, the trait values (z) of all species at time t.
  
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat[1, ] = init # initial trait values   
  for (t in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    z = z_mat[t, ] # current z values
    z_dif = t(f*z) - f*z # matrix with all trait differences
    q = f*(exp(-alpha * (z_dif^2))) # calculating matrix q
    q_n = q / apply(q, 1, sum) # normalizing the matrix
    q_m = q_n * m # multiplying each row i of matrix q by m[i]
    sel_dif = q_m * z_dif # calculating selection differentials
    r_mut = h * apply(sel_dif, 1, sum) # response to selection related to mutualism
    r_env = h * (1 - m) * (theta - z) # response to selection related to the environment
    z_mat[t+1, ] = z + r_mut + r_env # updating z values
    dif = mean(abs(z - z_mat[t+1, ])) # computing the mean difference between old and new z values
    if (dif < epsilon)
      break
  }
  return(z_mat[1:(t+1), ])
}

#-----------------------------------------------------------------------------------------------------#
