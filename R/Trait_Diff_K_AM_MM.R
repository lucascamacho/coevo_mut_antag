# Simulate a coevolution process without the AA interactions
# then, compute the mean trait value for groups of interaction types
# balancing this mean by the quantity of interactions

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
source("CoevoMutAntNet.R")
source("Antagonize.R")
source("EndInteraction.R")
source("Counting.R")
source("FindInteractors.R")
source("SpDegree.R")

library(ggplot2)
library(cowplot)

# initial conditions
antprob = 0.25  # current probability value
n_sp = 5   # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End antagonism AA
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

# measure the degree of AM and MM for each specie
degree = SpDegree(M, V)

# separate the species in groups by type of interaction
index = FindInteractors(M, V)

# coevolutionary model parameters
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

#create matriz for data
data = matrix(NA, nrow = nrow(z_mat), ncol = 5)
data[,5] = seq(1,nrow(z_mat), 1)
colnames(data) = c("AV_AM", "VAR_AM", "AV_MM", "VAR_MM", "time")

# for each line of z_mat
for(i in 1:nrow(z_mat)){
  # average and variance of all traits in the network
  av = mean(z_mat[i,])
  var = var(z_mat[i,])
  
  #calculating the average balanced by the degree and filling the data matriz
  data[i,1] = 2 # average of cheaters
  data[i,2] = 1 # variance of cheaters
  data[i,3] = 2 # average of mutualisms
  data[i,4] = 2 # variance of mutualisms

}

# building data frame to use in ggplot
data = data.frame(data)
test_data_long = melt(data, id = "time")  # convert to long format

# plot data
plotar = ggplot(data = test_data_long,
                aes(x = time, y = value, color = variable)) +
  geom_point(alpha = 0.7) +
  theme_bw()

plotar