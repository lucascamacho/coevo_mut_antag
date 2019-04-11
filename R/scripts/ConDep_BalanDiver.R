# This script will run the coevolutionary process with context-dependency
# Then, will calculate the Balançancia and Diverjância metrics
# Balançancia is a value for each specie and Diverjâncis for each timestep

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EndInteraction.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/BalanDiver.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
#antprob = 0.3 # current probability value
n_sp = 10 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform links in antagonisms)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End pure antagonism AA
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 7)
init = runif(n_sp, 0, 7)
p = 0.1
epsilon = 4
eq_dif = 0.0001
t_max = 1000
#prob_change = 0.2

# running coevolution simulation
simulation = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, 
                                  theta, init, p, epsilon, eq_dif, t_max, prob_change)
traits = simulation[[1]]
w_time = as.data.frame(simulation[[2]])

# apply the BalanDiver function
balandiver = BalanDiver(traits)
balan = balandiver[[1]]
diver = balandiver[[2]]

#set the times where the interaction shift occurs
#colnames(w_time) = "xplace"
#w_time$yplace = 1

# plot and save the diverjância by time
#time = seq(1, nrow(traits), 1)
#diver = data.frame(diver, time)

#diver_plot = ggplot(data = diver) +
#  geom_path(aes(x = time, y = diver)) +
#  geom_smooth(aes(x = time, y = ))
#  ggtitle(paste("Q =", prob_change, ", initial proportion of antagonists = ", antprob)) +
#  geom_text(data = w_time, aes(x=xplace, y=yplace),label = "*", size = 7) +
#  theme_bw()

#ggsave(diver_plot, filename = "Diver_Plot.png", width = 19, height = 11, units = "cm")
#diver_plot
#dev.off()