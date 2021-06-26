#### Figura 1 ####
# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(gridExtra)

# initial parameters
n_sp = 2 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive effects
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
V = M * 0
M[2,1] = 0
V[2,1] = 1

data = vector()
q_s = seq(from = 1, to = 0.01, by = -0.01)

for(k in 1:length(q_s)){  
  # coevolutionary model parameters
  phi = 0.2
  alpha = 0.2
  theta = c(2, 2) 
  init = c(2, 6)
  p = q_s[k]
  epsilon = 5
  eq_dif = 0.0001
  t_max = 1000
  
  CoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
    z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
    z_mat[1, ] = init # initial trait values
    
    for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
      to_ep = runif(1, 0, 1)
      z = z_mat[r, ] # current z values
      A = M + V # matrix with all interactions (mutualistic and exploitative)
      z_dif = t(A * z) - A * z # matrix with all trait differences
      Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
      diag(Q) = 0 # intraespecific effects are not allowed
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
      epsilon_equal = (z_dif == 0) * matrix(ifelse(to_ep > 0.5, epsilon, -epsilon), n_sp, n_sp)
      z_dif_a = z_dif + epsilon_plus + epsilon_minus + epsilon_equal
      sel_dif_ant = V_m * Q_m * z_dif_a # calculating selection differentials
      r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to exploitative
      
      z_mat[r+1, ] = z + r_mut + r_ant + r_env # updating z values
      
      dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
      if (dif < eq_dif) # if the difference is lower than eq_dif...
        break # stop simulation
      
    }
    
    return(z_mat[1:(r+1), ]) # return final matrix with species traits
    
  }
  
  # running coevolution simulation
  z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  
  disp = z_mat[nrow(z_mat),][1] - z_mat[nrow(z_mat),][2]
  data[k] = abs(disp)
} 

#### Forca ambiental é igual
x = seq(0.01, 1, 0.01)
e = 5
y = vector()

for(i in 1:length(x)){
  y[i] = (x[i]*e)  / (1 + x[i])
}

plot_data = data.frame(x,y,data)

plot_1 = ggplot(data = plot_data) +
  geom_line(aes(x = x, y = y)) +
  geom_point(aes(x = x, y = data)) +
  theme_classic() +
  ylab("Disparity") +
  xlab("Evolutive effects (qij)")
  #ggtitle("Considerando thetas iguais para i e j")

plot_1

#### Figura 2 ####
# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(gridExtra)
### Forca ambiental diferente
x = seq(0.01, 1, 0.01)
e = 5
y = vector()

theta_i = 4
theta_j = 2
d_theta = abs(theta_j - theta_i)

for(i in 1:length(x)){
  y[i] = ((x[i]*e) + ((1-x[i])*d_theta)) / (1+x[i])
}

# initial parameters
n_sp = 2 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive effects
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
V = M * 0
M[2,1] = 0
V[2,1] = 1

data = vector()
q_s = seq(from = 1, to = 0.01, by = -0.01)

for(k in 1:length(q_s)){  
  # coevolutionary model parameters
  phi = 0.2
  alpha = 0.2
  theta = c(2, 4) # AQUI
  init = c(2, 6)
  p = q_s[k]
  epsilon = 5
  eq_dif = 0.0001
  t_max = 1000
  
  CoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
    z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
    z_mat[1, ] = init # initial trait values
    
    for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
      to_ep = runif(1, 0, 1)
      z = z_mat[r, ] # current z values
      A = M + V # matrix with all interactions (mutualistic and exploitative)
      z_dif = t(A * z) - A * z # matrix with all trait differences
      Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
      diag(Q) = 0 # intraespecific effects are not allowed
      Q_n = Q / apply(Q, 1, sum) # normalizing the matrix
      Q_n[is.nan(Q_n)] = 0 # transform NaN values to 0 when a species don't have interactions
      Q_m = Q_n * (1 - p) # multiplying each row i of matrix Q by (1 - p)
      
      r_env = phi * p * (theta - z) # response to selection related to the environment
      
      sel_dif_mut = M * Q_m * z_dif # calculating selection differentials to mutualism
      r_mut = phi * apply(sel_dif_mut, 1, sum) # response to selection related to mutualism
      
      V_m = V # create V_m to use in exploitative selection differential 
      V_m[abs(z_dif) > epsilon] = 0
      epsilon_plus = (z_dif < 0) * matrix(epsilon, n_sp, n_sp) # matrix with barrier (epsilon) values
      epsilon_minus = (z_dif > 0) * matrix(-epsilon, n_sp, n_sp)
      epsilon_equal = (z_dif == 0) * matrix(ifelse(to_ep <= 0.5, epsilon, -epsilon), n_sp, n_sp)
      z_dif_a = z_dif + epsilon_plus + epsilon_minus + epsilon_equal 
      
      sel_dif_ant = V_m * Q_m * z_dif_a # calculating selection differentials
      r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to exploitative
      
      z_mat[r+1, ] = z + r_mut + r_ant + r_env # updating z values
      
      dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
      if (dif < eq_dif) # if the difference is lower than eq_dif...
        break # stop simulation
      
    }
    
    return(z_mat[1:(r+1), ]) # return final matrix with species traits
    
  }
  
  # running coevolution simulation
  z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  
  disp = z_mat[nrow(z_mat),][1] - z_mat[nrow(z_mat),][2]
  data[k] = abs(disp)
} 

plot_data = data.frame(x, y, data)

plot_2 = ggplot(data = plot_data) +
  geom_line(aes(x = x, y = y)) +
  geom_point(aes(x = x, y = data)) +
  theme_classic() +
  ylab("Disparity") +
  xlab("Evolutive effects (qij)")
  #ggtitle("Considerando thetas diferentes para i e j")

plot_2

#### Figura 3 ####
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(gridExtra)

# initial parameters
n_sp = 2 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive effects
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
V = M * 0
M[2,1] = 0
V[2,1] = 1

data = vector()
m = rev(seq(0.1, 0.8, 0.01))
v = seq(0.1, 0.8, 0.01)

for(k in 1:length(m)){  
  # coevolutionary model parameters
  phi = 0.2
  alpha = 0.2
  theta = c(2, 2)
  init = c(2, 4)
  p = 0.1
  epsilon = 5
  eq_dif = 0.0001
  t_max = 1000
  to_ep = runif(1, 0, 1)
  
  CoevoMutAntNet = function(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max) {
    z_mat = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
    z_mat[1, ] = init # initial trait values
    
    for (r in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
      z = z_mat[r, ] # current z values
      A = M + V # matrix with all interactions (mutualistic and exploitative)
      z_dif = t(A * z) - A * z # matrix with all trait differences
      Q = A * (exp(-alpha * (z_dif ^ 2))) # matrix Q
      diag(Q) = 0 # intraespecific effects are not allowed
      Q_n = Q / apply(Q, 1, sum) # normalizing the matrix
      Q_n[is.nan(Q_n)] = 0 # transform NaN values to 0 when a species don't have interactions
      #Q_m = Q_n * (1 - p) # multiplying each row i of matrix Q by (1 - p)
      
      r_env = phi * p * (theta - z) # response to selection related to the environment
      
      Q_m = Q_n * m[k] # define o mij
      sel_dif_mut = M * Q_m * z_dif # calculating selection differentials to mutualism
      r_mut = phi * apply(sel_dif_mut, 1, sum) # response to selection related to mutualism
      
      Q_m = Q_n * v[k] # define o vij
      V_m = V # create V_m to use in exploitative selection differential 
      V_m[abs(z_dif) > epsilon] = 0 # excluding interactions of traits that are larger than the barrier
      epsilon_plus = (z_dif < 0) * matrix(epsilon, n_sp, n_sp) # matrix with barrier (epsilon) values
      epsilon_minus = (z_dif > 0) * matrix(-epsilon, n_sp, n_sp) # matrix wih -epsilon values
      epsilon_equal = (z_dif == 0) * matrix(ifelse(to_ep > 0.5, epsilon, -epsilon), n_sp, n_sp) # matrix with epsilons from z_dif = 0  
      z_dif_a = z_dif + epsilon_plus + epsilon_minus + epsilon_equal # adding barrier values to trait differences
      sel_dif_ant = V_m * Q_m * z_dif_a # calculating selection differentials
      r_ant = phi * apply(sel_dif_ant, 1, sum) # response to selection related to exploitative
      
      z_mat[r+1, ] = z + r_mut + r_ant + r_env # updating z values
      
      dif = mean(abs(z - z_mat[r+1, ])) # computing the mean difference between old and new z values
      if (dif < eq_dif) # if the difference is lower than eq_dif...
        break # stop simulation
      
    }
    
    return(z_mat[1:(r+1), ]) # return final matrix with species traits
    
  }
  
  # running coevolution simulation
  z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  
  disp = z_mat[nrow(z_mat),][1] - z_mat[nrow(z_mat),][2]
  data[k] = abs(disp)
}

plot_data = data.frame(m, v, data)

plot_3 = ggplot(data = plot_data) +
  geom_point(aes(x = m, y = data)) +
  theme_classic() +
  ylab("Disparity") +
  xlab("Importance of mutualism (mij)")
  #ggtitle("Considerando thetas iguais e mij e vij diferentes")

plot_3



plot_final = grid.arrange(plot_1, plot_2, plot_3, nrow = 3)

ggsave(plot_final, filename = "analitical_tests.png", dpi = 300,
       width = 18, height = 21, units = "cm",  bg = "transparent")