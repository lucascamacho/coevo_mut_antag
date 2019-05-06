# Simulate a coevolution process without the AA interactions and
# calculate the mean trait value for groups of interaction outcomes.
#
# The mean trait value are balanced by the number of interaction outcomes
# (AM and MM) that a certain species have.
#
# This scripts returns 2 graphs. One with the average species traits balanced by degree
# and another of species traits variance. Both varying in time. The blue dots are the MM
# and the red dots are the AM species.
#
# The script probably will have some # symbols in parameters like antprob.
# If you will use only this code without a loop, be sure to get the #'s off.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SpDegree.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/FindInteractors.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)

# initial conditions
#antprob = 0.9 # current probability value
n_sp = 50 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# measure the degree of AM and MM for each specie
degree = SpDegree(M, V)

# count the number of double-outcome interactions
c = Counting(M, V)

# separate the species in groups by type of interaction outcomes
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

# apply the calculations of our metrics in z_mat
dif_traits = apply(z_mat, 1, MeanPairDist)
var_traits = apply(z_mat, 1, var)
time_traits = seq(1, nrow(z_mat), 1)

# create data frames to each variable
var_time = data.frame(time_traits, var_traits)
dif_time_AM = data.frame(time_traits, dif_traits[1,])
dif_time_MM = data.frame(time_traits, dif_traits[2,])

# plot the results
#var_plot = ggplot(var_time) + 
#           geom_point(aes(x = time_traits, y = var_traits)) + 
#           ggtitle(paste("proportion antagonists = ", antprob)) +
#           xlab("Time") + 
#           ylab("Variance of species traits (z)") +
#           theme_bw()
#
#dif_plot = ggplot() +
#           geom_point(aes(x = dif_time_AM[,1], y = dif_time_AM[,2]), color = "red") +
#           geom_point(aes(x = dif_time_MM[,1], y = dif_time_MM[,2]), color = "blue") +
#           ggtitle(paste("proportion antagonists = ", antprob)) +
#           xlab("Time") + 
#           ylab("Mean Pairwise Distance separated by signal effect of interaction") +
#          theme_bw()
#
#var_plot
#dif_plot