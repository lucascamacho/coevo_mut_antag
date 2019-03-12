# Simulate a coevolution process without the AA interactions
# then, compute the mean trait value for groups of interaction types
# balancing this mean and variance of traits by the quantity of interactions.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EndInteraction.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SpDegree.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/FindInteractors.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/VarTraitDegreeBalanced.R")
library(ggplot2)

# initial conditions
#antprob = 0.9 # current probability value
n_sp = 30 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform links in antagonisms)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End the AA interactions
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

# measure the degree of AM and MM for each specie
degree = SpDegree(M, V)

# count the number of interactions
c = Counting(M, V)

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

# apply the calculations of our metrics in z_mat
dif_traits = apply(z_mat, 1, VarTraitDegreeBalanced)
var_traits = apply(z_mat, 1, var)
time_traits = seq(1, nrow(z_mat), 1)

# create data frames to each variable
var_time = data.frame(time_traits, var_traits)
dif_time_AM = data.frame(time_traits, dif_traits[1,])
dif_time_MM = data.frame(time_traits, dif_traits[2,])

# ggplots to variables (off the #'s for tests)
#var_plot = ggplot(data = var_time) + 
#           geom_point(aes(x = time_traits, y = var_traits), alpha = 0.7) + 
#           theme_bw()

#dif_plot = ggplot() +
#           geom_point(aes(x = time_traits, y = dif_time_AM[,2]), color = "red") +
#           geom_point(aes(x = time_traits, y = dif_time_MM[,2]), color = "blue") +
#           theme_bw()

#var_plot
#dif_plot