# For empirical networks:
# Script to simulate the coevolution process without the AA interactions outcomes,
# and compute the mean trait value for groups of interaction outcomes types.
#
# These means are balanced by their frequency in the network Kam and Kmm.
#
# This script returns a plot with the mean trait value and variance of traits
# for reach interaction outcome type AM and MM.
#
# The script probably will have some # symbols in parameters like antprob.
# If you will use only this code without a loop, be sure to get the #'s off.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EmpAntagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SpDegree.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/FindInteractors.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/TraitDegreeBalanced.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/VarTraitDegreeBalanced.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial conditions
#antprob = 0.4 # current probability value
net = "~/Dropbox/Master/Code/coevo_mut_antag/data/B_NS-PS-Galetti&Pizo-1996-SG2Wgt.txt"
M = as.matrix(read.table(net)) # read empirical matrix
M = SquareMatrix(M) # Square the M matrix
n_sp = dim(M)[1] # defining number of species

# Antagonize M (transform positive links in negative)
empantagonize = EmpAntagonize(M, antprob)
M = empantagonize[[1]]
V = empantagonize[[2]]

# measure the degree of positive and negative outcomes for each specie
degree = SpDegree(M, V)

# count the number of double-interactions
c = Counting(M, V)

# separate the species in groups by type of interaction outcome
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

# simulate coevolution
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# create matriz for data
data = matrix(NA, nrow = nrow(z_mat), ncol = 5)
data[,5] = seq(1,nrow(z_mat), 1)
colnames(data) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM", "time")

# apply function to get the mean and variance of traits balanced by degree Kam and Kmm
traits = apply(z_mat, 1, TraitDegreeBalanced)
var_traits = apply(z_mat, 1, VarTraitDegreeBalanced)

# allocate the results in data matrix
data[,1] = traits[1,]
data[,2] = var_traits[1,]
data[,3] = traits[2,]
data[,4] = var_traits[2,]

# plot the results
#data = data.frame(data)
#par(mfrow=c(2,2))
#plot(data$time, data$MEAN_AM, col="blue", pch = 19, xlab = "time", ylab = "Mean Trait for Antagonism")
#plot(data$time, data$MEAN_MM, col="blue", pch = 19, xlab = "time", ylab = "Mean Trait for Mutualism")
#plot(data$time, data$VAR_AM, col="red", pch = 19, xlab = "time", ylab = "Delta Trait for Antagonism")
#plot(data$time, data$VAR_MM, col="red", pch = 19, xlab = "time", ylab = "Delta Trait for Mutualisms")
#title("Traits Dynamics of Antagonism and Mutualism (Balanced by degree Kam and Kmm)", line = -2, outer = TRUE)