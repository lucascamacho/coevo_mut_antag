#-------------------------------------------------------------------------------------------------------#
DiscoCoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max, bar) {
  # Simulates the coevolutionary dynamics of mutualists and exploitative interactions in a network taking
  # into account the trait difference between species. If they are too much different, the interactions
  # is a Forbidden Link.
  #
  # Args:
  #   n_sp: total number of species in the network
  #   M: square matrix representing the mutualistic interactions
  #   V: square matrix representing the exploitative interactions
  #   phi: vector of phi values (additive genetic variance * heritability)
  #   alpha: alpha value (sensitivity of selection to trait matching)
  #   theta: vector of environmental optimum values
  #   init: vector of initial trait values
  #   p: vector of strength of selection due to the environment
  #   epsilon: barrier value for exploitative interactions happend
  #   eq_dif: value to determine when equilibrium is reached
  #   t_max: maximum number of timesteps allowed
  #   bar: trait matching barrier (for both interaction outcomes)
  #   
  # Obs: 
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   z_mat: a matrix containing, in each row t, the trait values (z) of all species at time t.
  #   m_final: list with all the adjacency matrix of interactions in cronological time
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat[1, ] = init # initial trait values
  m_final = list()
  m_final[[1]] = M + V
  
  for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    z = z_mat[r, ] # current z values
    A = M + V # matrix with all interactions (mutualistic and exploitative)
    z_dif = t(A * z) - A * z # matrix with all trait differences
    d = which(abs(z_dif) > bar) # which species trait diff is higher than bar?
    A[d] = 0 # transform this interactions in Forbidden Links
    #d_inverse = which(abs(z_dif) < bar) # which species trait diff is lower than bar?
    #ind = sample(d_inverse, ceiling(((sum(A) * 5) / 100))) # sample 5% of this interactions
    #A[ind] = 1 # reconnect species with trait similarity
    Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
    diag(Q) = 0 # intraespecific effects are not allowed
    degree = rowSums(A)
    abnd = rlnorm(n_sp, meanlog = 1, sdlog = 1) #ligar/desligar quando quiser
    abnd = abnd[order(match(rank(abnd), rank(degree)))] # cor abnd and degree
    Q = t(abnd * Q)
    Q_n = Q / apply(Q, 1, sum) # normalizing the matrix
    Q_n[is.nan(Q_n)] = 0 # transform NaN values to 0 when a species don't have interactions
    Q_m = Q_n * (1 - p) # multiplying each row i of matrix Q by (1 - p)
    
    r_env = phi * p * (theta - z) # response to selection related to the environment
    
    sel_dif_mut = M * Q_m * z_dif # calculating selection differentials to mutualism
    r_mut = phi * apply(sel_dif_mut, 1, sum) # response to selection related to mutualism
    
    V_m = V # create V_m to use in exploitative selection differential 
    V_m[abs(z_dif) > epsilon] = 0 # excluding interactions of traits that are larger than the barrier
    epsilon_plus = (z_dif < 0) * matrix(epsilon, n_sp, n_sp) # matrix with barrier (epsilon) values
    epsilon_minus = (z_dif > 0) * matrix(-epsilon, n_sp, n_sp) # matrix wih -epsilon values
    z_dif_a = z_dif + epsilon_plus + epsilon_minus # adding barrier values to trait differences
    sel_dif_ant = V_m * Q_m * z_dif_a # calculating selection differentials
    r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to exploitative interactions
    
    z_mat[r+1, ] = z + r_mut + r_ant + r_env # updating z values
    m_final[[r+1]] = A # insert the current adjancency matrix in a list
        
    dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
    if (dif < eq_dif) # if the difference is lower than eq_dif...
      break # stop simulation
    
  }
  
  f_coevo = list(z_mat[1:(r+1), ], m_final)
  return(f_coevo) # return final matrix with species traits
  
}
#-------------------------------------------------------------------------------------------------------#