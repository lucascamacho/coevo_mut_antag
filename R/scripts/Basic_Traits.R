# Basic script to test the R base-functions to simulate the coevolutionary process
# of networks with positive and negative effects in a single network.
# This script returns a simple graph with species traits changing in time.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")

library(ggplot2)
library(cowplot)
library(reshape2)

# initial parameters
antprob = 0.8 # current probability value
n_sp = 10 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive effects
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
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# building data frame to plot the results
traits = as.data.frame(z_mat)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                       time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)

# plotting traits through time
plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
  geom_path(size = 1.8, alpha = 0.7, show.legend = FALSE) + 
  ggtitle(paste("proportion of exploitation = ", antprob)) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))
plotar

save(z_mat, file = "dados_plot.RData")

# save plot
ggsave(plotar, filename = "Basic_Traits.pdf", dpi = 600, width = 16, height = 12, 
       units = "cm", bg = "transparent")

#####

# Basic script to test the R base-functions to simulate the coevolutionary process
# of networks with positive and negative effects in a single network.
# This script returns a simple graph with species traits changing in time.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")

library(ggplot2)
library(cowplot)
library(reshape2)

# initial parameters
n_sp = 2 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive effects
diag(M) = 0 # no intraespecific interactions
M[1,2] = 0

V = M * 0
V[1,2] = 1

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.5
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# running coevolution simulation
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# building data frame to plot the results
traits = as.data.frame(z_mat)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                       time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)

# plotting traits through time
plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
  geom_path(size = 1.8, alpha = 0.7, show.legend = FALSE) + 
  #  ggtitle(paste("proportion of exploitation = ", antprob)) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))
plotar


