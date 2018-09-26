# Simulate a coevolution process without the AA interactions
# then, compute the mean trait value for groups of interaction types
# balancing this mean and variance of traits by the quantity of interactions.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
source("Antagonize.R")
source("EndInteraction.R")
source("ZeroLines.R")
source("SpDegree.R")
source("FindInteractors.R")
source("CoevoMutAntNet.R")
source("TraitDegreeBalanced.R")
source("VarTraitDegreeBalanced.R")

library(ggplot2)
library(reshape)
library(cowplot)

# initial conditions
#antprob = 0.25  # current probability value
n_sp = 50   # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform links in antagonisms)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End the AA interactions
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

# check for zero lines (otherwise the simulation returns an error)
zero = ZeroLines(M, V, n_sp, antprob)
M = zero[[1]]
V = zero[[2]]

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

# simulate coevolution and calculate mean and variance trait values for each interaction type
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# create matriz for data
data = matrix(NA, nrow = nrow(z_mat), ncol = 5)
data[,5] = seq(1,nrow(z_mat), 1)
colnames(data) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM", "time")

# apply function to get the mean and variance of traits balanced by degree
traits = apply(z_mat, 1, TraitDegreeBalanced)
var_traits = apply(z_mat, 1, VarTraitDegreeBalanced)

# allocate the results in data matrix
data[,1] = traits[1,]
data[,2] = var_traits[1,]
data[,3] = traits[2,]
data[,4] = var_traits[2,]

# plot the results
data = data.frame(data)
par(mfrow=c(2,2))
plot(data$time, data$MEAN_AM, col="red", pch = 19, xlab = "time", ylab = "Mean Trait for Cheaters")
plot(data$time, data$VAR_AM, col="red", pch = 19, xlab = "time", ylab = "Mean Trait for Mutualism")
plot(data$time, data$MEAN_MM, col="red", pch = 19, xlab = "time", ylab = "Delta Trait for Cheaters")
plot(data$time, data$VAR_MM, col="red", pch = 19, xlab = "time", ylab = "Delta Trait for Cheaters")
title("Traits Dynamics of Cheaters and Mutualism (Balanced by degree Kmm and Kam)", line = -2, outer = TRUE)