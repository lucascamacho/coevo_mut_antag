#-----------------------------------------------------------------------------------------------------#
rm(list=ls())
### testing coevolutionary model ##

# loading packages and functions
library(ggplot2)
library(cowplot)
source("CoevoMutAntNet.R")

# defining probability of link becoming antagonist
antprob_vec = seq(0.1, 0.9, by = 0.1)

for (i in 1:length(antprob_vec)) {
  # current probability value
  antprob = antprob_vec[i]
  # defining number of species
  n_sp = 5
  # building matrix M (mutualisms)
  M = matrix(1, ncol = n_sp, nrow = n_sp)
  diag(M) = 0
  # building matrix V (antagonisms)
  V = M * 0
  P = matrix(runif(n_sp*n_sp, min = 0, max = 1),
             nrow = n_sp, ncol = n_sp)
  V[antprob >= P] = 1
  M[antprob >= P] = 0
  # coevolutionary model parameters
  phi = 0.2
  alpha = 0.2
  theta = c(3, 7, 4, 4, 8)
  init = c(10, 2, 8, 9, 5)
  p = 0.8
  epsilon = 4
  eq_dif = 0.00001
  t_max = 100
  # running coevolution simulation
  traits = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  print(tail(traits))
  # building data frame to use in ggplot
  traits = as.data.frame(traits)
  n_sp = ncol(traits)
  traits_vec = c(as.matrix(traits))
  traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                         time = rep(1:nrow(traits), times = n_sp),
                         trait = traits_vec)
  # plotting traits through time
  p = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
    geom_point(size = 1.8, shape = 19, alpha = 0.6) + 
    ggtitle(paste("proportion antagonists = ", antprob)) +
    xlab("Time") + 
    ylab("Mean species trait (z)") +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 14), 
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size = 12))

  print(p)
}

#-----------------------------------------------------------------------------------------------------#