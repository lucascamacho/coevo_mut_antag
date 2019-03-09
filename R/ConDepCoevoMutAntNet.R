#-----------------------------------------------------------------------------------------------------#
ConDepCoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
  # Simulates the coevolutionary dynamics of mutualists and antagonists in a network
  #
  # Args:
  #   n_sp: total number of species in the network
  #   M: square adjacency matrix representing the mutualistic interactions
  #   V: square adjacency matrix representing the antagonistic interactions
  #   phi: vector of phi values (additive genetic variance * heritability)
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
  source("MutualizeAntagonize.R")
  
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat[1, ] = init # initial trait values
  w_t_vm = vector() #vectors to track in which timestep occurs the interaction changes
  w_t_mv = vector()
  for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    mutualizeantagonize = MutualizeAntagonize(M, V, r) # changing interactions and track the timestep
    M = mutualizeantagonize[[1]]
    V = mutualizeantagonize[[2]]
    w_t_vm = append(w_t_vm, mutualizeantagonize[[3]])
    w_t_mv = append(w_t_mv, mutualizeantagonize[[4]])
    
    z = z_mat[r, ] # current z values
    A = M + V # matrix with all interactions (mutualistic and antagonistic)
    z_dif = t(A * z) - A * z # matrix with all trait differences
    Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
    diag(Q) = 0 # intraespecific effects are not allowed
    Q_n = Q / apply(Q, 1, sum) # normalizing the matrix
    Q_n[is.nan(Q_n)] = 0 # transform NaN values to 0 when a species don't have interactions
    Q_m = Q_n * (1 - p) # multiplying each row i of matrix Q by (1 - p)
    
    r_env = phi * p * (theta - z) # response to selection related to the environment
    
    sel_dif_mut = M * Q_m * z_dif # calculating selection differentials to mutualism
    r_mut = phi * apply(sel_dif_mut, 1, sum) # response to selection related to mutualism
    
    V[abs(z_dif) > epsilon] = 0 # excluding interactions of traits that are larger than the barrier
    epsilon_plus = (z_dif < 0) * matrix(epsilon, n_sp, n_sp) # matrix with barrier (epsilon) values
    epsilon_minus = (z_dif > 0) * matrix(-epsilon, n_sp, n_sp) # matrix wih -epsilon values
    z_dif = z_dif + epsilon_plus + epsilon_minus # adding barrier values to trait differences
    sel_dif_ant = V * Q_m * z_dif # calculating selection differentials
    r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to antagonisms
    
    z_mat[r+1, ] = z + r_env + r_mut + r_ant # updating z values
    
    dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
    if (dif < eq_dif) # if the difference is lower than eq_dif...
      break # stop simulation
    
  }
  
  return(list(z_mat[1:(r+1), ], w_t_vm, w_t_mv)) # return final matrix with species traits
  
}

#-----------------------------------------------------------------------------------------------------#