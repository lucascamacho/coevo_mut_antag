# Script to run the coevolution in a theoretical network of interactions
# with context-dependent interactions, plotting the species traits in time
# Also, this script calculates our 4 metrics of trait disparity for each timestep:
#  1. Variance
#  2. Mean Pairwise Distance
#  3. Participation Ratio
#  4. Mean Nearest and Distant Neighbor Distance
#
# This script returns one graph and the four measures of trait disparity
# in the same order as described above.
#
# The script probably will have some # symbols in parameters like antprob.
# If you will use only this code without a loop, be sure to get the #'s off.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/PartRatio.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/NearDist.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
antprob = 0.2 # current probability value
prob_change = 0.2 # current probability of interaction outcome shift
n_sp = 10 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 5)
init = runif(n_sp, 0, 5)
p = 0.1
epsilon = 3
eq_dif = 0.0001
t_max = 1000

# running coevolution simulation
traits = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, 
                              eq_dif, t_max, prob_change)

# calculate our 4 trait disparity metrics #APPLY DOESN'T WORK
variance = apply(traits, 1, var)
meanpairdist = apply(traits, 1, MeanPairDist)
partratio = apply(as.matrix(t(traits)), 1, PartRatio)
neardist = apply(traits, 1, NearDist)

#set the times where the interaction oucomes shift occurs
colnames(w_time) = "xplace"
w_time$yplace = 1

# plot the Divergency in time
time = seq(1, nrow(traits), 1)
diver = data.frame(diver, time)

diver_plot = ggplot(data = diver) +
  geom_path(aes(x = time, y = diver)) +
  ggtitle(paste("Q =", prob_change, ", initial proportion of antagonists = ", antprob)) +
  geom_text(data = w_time, aes(x = xplace, y = yplace),label = "*", size = 7) +
  theme_bw()

#ggsave(diver_plot, filename = "Diver_Plot.png", width = 19, height = 11, units = "cm")
diver_plot
#dev.off()