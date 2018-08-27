setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
library(ggplot2)
source("MeanTraitInteract.R")

#-----------------------------------------------------------------------------------------------------#
CoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
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
  z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat[1, ] = init # initial trait values
  for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    z = z_mat[r, ] # current z values
    A = M + V # matrix with all interactions (mutualistic and antagonistic)
    z_dif = t(A * z) - A * z # matrix with all trait differences
    Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
    diag(Q) = 0 # intraespecific effects are not allowed
    Q_n = Q / apply(Q, 1, sum) # normalizing the matrix
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
    if (dif < eq_dif)
      break
    
    diff_interactions = MeanTraitInteract(M, V, z_mat[r+1])
    data[r+1, ] = diff_interactions
    
  }
  
  return(data)
}

# initial conditions
n_sp = 10 
antprob = 0.5

# creating matrix of interactions
M = matrix(1, ncol = n_sp, nrow = n_sp)
diag(M) = 0
V = M * 0

P = matrix(runif(n_sp*n_sp, min = 0, max = 1),
           nrow = n_sp, ncol = n_sp)
V[antprob >= P] = 1
M[antprob >= P] = 0
diag(V) <- 0

# coevolution model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# building data frame to use in ggplot
data = matrix(NA, nrow = t_max, ncol = 3)
colnames(data) = c("AA", "AM", "MM")

# simulate coevolution
traits = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

#prepare data for plot
