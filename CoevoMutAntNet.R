#-----------------------------------------------------------------------------------------------------#

CoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
  # Simulates the coevolutionary dynamics of mutualists and antagonists in a network
  #
  # Args:
  #   n_sp: total number of species in the network
  #   M: square adjacency matrix representing the mutualistic interactions
  #   V: square adjacency matrix representing the antagonistic interactions
  #   phi: vector of phi values (additive genetic variance*heritability)
  #   alpha: alpha value (sensitivity of selection to trait matching)
  #   theta: vector of environmental optimum values
  #   init: vector of initial trait values
  #   p: vector of strength of selection due to the environment
  #   epsilon: barrier value for antagonistic interactions
  #   eq_dif: value to determine when equilibrium is reached
  #   t_max: maximum number of timesteps allowed
  #   
  # Obs: 
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   A matrix containing, in each row t, the trait values (z) of all species at time t.
  
  # matrix to store z values
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp)
  # initial trait values   
  z_mat[1, ] = init
  # simulation runs for a maximum of t_max timesteps
  for (x in 1:(100 - 1)) { 
    # current z values
    z = z_mat[x, ]
    # matrix with all interactions (mutualistic and antagonistic)
    A = M + V
    # matrix with all trait differences 
    z_dif = t(A*z) - A*z
    # matrix Q
    Q = A*(exp(-alpha * (z_dif^2)))
    # normalizing the matrix
    Q_n = Q / apply(Q, 1, sum) 
    # multiplying each row i of matrix Q by (1 - p[i])
    Q_m = Q_n * (1 - p) 
    
    ### response to selection - environment ###
    # response to selection related to the environment
    r_env = phi * p * (theta - z)
    
    ### response to selection - mutualism ###
    # calculating selection differentials 
    sel_dif = Q_m * z_dif 
    # response to selection related to mutualism
    r_mut = phi * apply(sel_dif, 1, sum)
 
    ### response to selection - antagonism ###
    # excluding values that are larger than the barrier
    z_dif[abs(z_dif) > epsilon] = 0
    # matrix with barrier (epsilon) values 
    epsilon_plus = (z_dif < 0)*matrix(epsilon, n_sp, n_sp)
    epsilon_minus = (z_dif > 0)*matrix(-epsilon, n_sp, n_sp)
    # adding barrier values to trait differences
    z_dif = z_dif + epsilon_plus + epsilon_minus
    # calculating selection differentials
    sel_dif = Q_m * z_dif
    # response to selection related to antagonisms 
    r_ant = phi * apply(sel_dif, 1, sum)
    # updating z values
    z_mat[x+1, ] = z + r_env + r_mut + r_ant
    # computing the mean difference between old and new z values
    dif = mean(abs(z - z_mat[x+1, ])) 
    if (dif < eq_dif)
      break
  }
  return(z_mat[1:(x+1), ])
}

#-----------------------------------------------------------------------------------------------------#