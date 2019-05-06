# This script will run the coevolutionary process with context-dependency of interactions.
# Then, will calculate the Divergency and Directionality of species traits.
# Finally, the script will plot the Divergency by Time.

# The script probably will have some # symbols in parameters like antprob or prob_change.
# If you will use only this code without a loop, be sure to get the #'s off.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/BalanDiver.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
#antprob = 0.1 # current probability value
#prob_change = 0.2 # current probability of interaction outcome shift
n_sp = 50 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# running coevolution simulation
simulation = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, 
                                  theta, init, p, epsilon, eq_dif, t_max, prob_change)
traits = simulation[[1]]
w_time = as.data.frame(simulation[[2]])

# apply the BalanDiver function
balandiver = BalanDiver(traits)
balan = balandiver[[1]]
diver = balandiver[[2]]

#set the times where the interaction oucomes shift occurs
#colnames(w_time) = "xplace"
#w_time$yplace = 1

# plot the Divergency in time
#time = seq(1, nrow(traits), 1)
#diver = data.frame(diver, time)

#diver_plot = ggplot(data = diver) +
#  geom_path(aes(x = time, y = diver)) +
#  ggtitle(paste("Q =", prob_change, ", initial proportion of antagonists = ", antprob)) +
#  geom_text(data = w_time, aes(x = xplace, y = yplace),label = "*", size = 7) +
#  theme_bw()

#ggsave(diver_plot, filename = "Diver_Plot.png", width = 19, height = 11, units = "cm")
#diver_plot
#dev.off()