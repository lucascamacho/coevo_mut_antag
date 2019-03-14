# Basic script to run the coevolutionary model with context dependendy of interactions.
# this script returns a simple graph with species traits changing in time
# due to coevolution. The stars in the graph shows the timesteps in which the interactions
# shift occurs.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EndInteraction.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
antprob = 0.1 # current probability value
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
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# running coevolution simulation
simulation = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, 
                                  theta, init, p, epsilon, eq_dif, t_max)
traits = simulation[[1]]
w_time = as.data.frame(simulation[[2]])

#prepare data frame with tracked timesteps
colnames(w_time) = "xplace"
w_time$yplace = 1

# building data frame to use in ggplot2
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
#  geom_text(data = w_time, aes(x=xplace, y=yplace),label = "*", size = 7) +
  ggtitle(paste("Q = 1%, initial proportion of antagonists = ", antprob)) +
  geom_text(data = w_time, aes(x=xplace, y=yplace),label = "*", size = 7) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))

#ggsave(plotar, filename = "10_Basic_Traits.png", width = 19, height = 11, units = "cm")
plotar
#dev.off()