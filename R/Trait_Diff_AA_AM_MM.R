# This code create matrix of interactions of mutualism and antagonism,
# simulate the coevolutionary process of the species and
# calculate the mean trait value for cheaters, pure mutualistis and antagonists
#
# Load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
source("CoevoMutAntNet.R")
source("FindInteractors.R")
library(ggplot2)
library(reshape2)

# initial conditions
n_sp = 10
#antprob = 0.2

# creating matrix of interactions
M = matrix(1, ncol = n_sp, nrow = n_sp)
diag(M) = 0
V = M * 0

# changing interactions from mutualism to antagonism
P = matrix(runif(n_sp*n_sp, min = 0, max = 1),
           nrow = n_sp, ncol = n_sp)
V[antprob >= P] = 1
M[antprob >= P] = 0
diag(V) = 0

# coevolution model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# simulate coevolution and calculate mean trait values for each interaction type
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
index = FindInteractors(M, V)

data = matrix(NA, nrow = nrow(z_mat), ncol = 4)
data[,4] = seq(1,nrow(z_mat), 1)
colnames(data) = c("AA", "AM", "MM", "time")

for(i in 1:nrow(z_mat)){
  data[i,1] = mean(z_mat[i,index[[1]]])
  data[i,2] = mean(z_mat[i,index[[2]]])
  data[i,3] = mean(z_mat[i,index[[3]]])
  
}

# prepare to plot and plotting the results
#data = data.frame(data)
#test_data_long = melt(data, id = "time")  # convert to long format
#
#plotar = ggplot(data = test_data_long,
#                aes(x = time, y = value, colour = variable)) +
#         geom_point(alpha = 0.4) +
#         theme_bw()
