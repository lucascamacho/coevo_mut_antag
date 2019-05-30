#-------------------------------------------------------------------------------------------------------------#
ConDepCoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max, prob_change){
  # Simulates the coevolutionary dynamics of mutualists and antagonists outcomes in a network
  # with context dependent interactions. The interactions change in time following a 
  # certain probability prob_change.
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
  #   prob_change: probability that in each timestep the interaction shift between M and V
  #   
  # Obs: 
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   A matrix containing, in each row t, the trait values (z) of all species at time t.
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MutAntag.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
  
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat[1, ] = init # initial trait values
  w_time = vector() #vectors to track in which timestep occurs the interaction changes
  c = list()
  
  for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    z = z_mat[r, ] # current z values
    
    mutantag = MutAntag(M, V, r, prob_change) # run the outcomes changer
    M = mutantag[[1]] # define new M matrix
    V = mutantag[[2]] # define new V matrix
    w_time = append(w_time, mutantag[[3]]) # vector with timesteps of change
    
    c[[r]] = Counting(M, V)
    
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
    
    V_m = V # create V_m to use in antagonism selection differential 
    V_m[abs(z_dif) > epsilon] = 0 # excluding interactions of traits that are larger than the barrier
    epsilon_plus = (z_dif < 0) * matrix(epsilon, n_sp, n_sp) # matrix with barrier (epsilon) values
    epsilon_minus = (z_dif > 0) * matrix(-epsilon, n_sp, n_sp) # matrix wih -epsilon values
    z_dif_a = z_dif + epsilon_plus + epsilon_minus # adding barrier values to trait differences
    sel_dif_ant = V * Q_m * z_dif_a # calculating selection differentials
    r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to antagonisms
    
    z_mat[r+1, ] = z + r_mut + r_ant + r_env # updating z values
    
    dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
    if (dif < eq_dif) # if the difference is lower than eq_dif...
      break # stop simulation
    
  }
  
  return(list(z_mat[1:(r+1), ], w_time, c)) # return final matrix with species traits
  
}

#-----------------------------------------------------------------------------------------------------#