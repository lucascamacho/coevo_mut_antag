### For empirical networks of interactions only
# Script to run the coevolution in a empirical network of interactions
# with interactions shifts in time, plotting the species traits in time
# Also, this script calculates our 2 metrics of trait disparity for each timestep:
#  1. Variance
#  2. Mean Pairwise Distance
#
# This script returns one graph and the four measures of trait disparity
# in the same order as described above.
#
# The script probably will have some # symbols in parameters like antprob.
# If you will use only this code without a loop, be sure to get the #'s off.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CentralAntagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")


library(ggplot2)
library(cowplot)

# initial parameters
antprob = 0.5 # current probability value
prob_change = 0.01 # current probability of interaction outcome shift

# read and square the empirical network
net = as.matrix(read.table(
  "~/Dropbox/Master/Code/coevo_mut_antag/data/B_SY-AP-IzzomtcCOLpen.txt"))
M = SquareMatrix(net)
n_sp = ncol(M)

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
simulation = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, 
                                  eq_dif, t_max, prob_change)
traits = as.matrix(simulation[[1]])
w_time = as.matrix(simulation[[2]])

# an apply for each line of z_mat
variance = apply(traits, 1, var)
meanpairdist = apply(traits, 1, MeanPairDist)

# set the times where the interaction oucomes shift occurs
colnames(w_time) = "xplace"
yplace = rep(1, nrow(w_time))
w_time = cbind(w_time, yplace)
w_time = as.data.frame(w_time)

# plot the Divergency in time
time = seq(1, nrow(traits), 1)
diver = data.frame(variance, time)

diver_plot = ggplot(data = diver) +
  geom_path(aes(x = time, y = variance)) +
  theme_bw()


# prepare data frame to plot
traits = as.data.frame(traits)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                      time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)

# plotting traits through time
plotar = ggplot() +
  geom_path(data=traits_df, aes(x = time, y = trait, group=species, 
                                color = species),size = 1.8, alpha = 0.7) +
  geom_text(data = w_time, aes(x=xplace, y=yplace),label = "*", size = 7) +
  ggtitle(paste("G =", prob_change, ", P = ", antprob)) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))

plotar

#save plot
#ggsave(plotar, filename = "ConDep_Basic_Traits.pdf", width = 19, height = 11, units = "cm")