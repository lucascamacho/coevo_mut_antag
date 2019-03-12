# Calculate the mean trait value for cheaters, pure mutualists and pure antagonists.
# We use the FindInteractors function to separate species in 3 groups (AA, AM and MM).
# Then, we calculate the mean trait for each group and plot these means

# Load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/FindInteractors.R")

library(ggplot2)
library(reshape2)

# initial conditions
n_sp = 10
antprob = 0.25

# creating matrix of interactions
M = matrix(1, ncol = n_sp, nrow = n_sp)
diag(M) = 0

# Antagonize M (transform links in antagonisms)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# coevolution model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# simulate coevolution and find the cheaters, mutualists and antagonists
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
index = FindInteractors(M, V)

# create a data matrix
data = matrix(NA, nrow = nrow(z_mat), ncol = 4)
data[,4] = seq(1,nrow(z_mat), 1)
colnames(data) = c("AA", "AM", "MM", "time")

# loop to calculate the mean trait value of each group of species
for(i in 1:nrow(z_mat)){
  data[i,1] = mean(z_mat[i,index[[1]]])
  data[i,2] = mean(z_mat[i,index[[2]]])
  data[i,3] = mean(z_mat[i,index[[3]]])
  
}

# prepare to plot and plotting the results
data = data.frame(data)
test_data_long = melt(data, id = "time")  # convert to long format

plotar = ggplot(data = test_data_long,
                aes(x = time, y = value, colour = variable)) +
         geom_point(alpha = 0.4) +
         theme_bw()
plotar