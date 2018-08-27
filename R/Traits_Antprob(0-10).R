# testing coevolutionary model
# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

library(ggplot2)
library(cowplot)
source("CoevoMutAntNet.R")

# defining probability of link becoming antagonist
antprob_vec = seq(0.1, 1, 0.1)

for (i in 1:length(antprob_vec)) {
  # creating matrices
  antprob = antprob_vec[i]  # current probability value
  n_sp = 5   # defining number of species
  M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
  diag(M) = 0 # no intraespecific interactions
  V = M * 0   # building matrix V (antagonisms)
  
  P = matrix(runif(n_sp*n_sp, min = 0, max = 1),
             nrow = n_sp, ncol = n_sp) # P is a matrix of sort values to compare with antprob
  V[antprob >= P] = 1 # if value of P is higher, interaction is antagonistic
  M[antprob >= P] = 0 # and no mutualistic
  diag(V) <- 0 # no intraespecific interactions
  
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
  traits = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

  # building data frame to use in ggplot
  traits = as.data.frame(traits)
  n_sp = ncol(traits)
  traits_vec = c(as.matrix(traits))
  traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                         time = rep(1:nrow(traits), times = n_sp),
                         trait = traits_vec)
  # plotting traits through time
  plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
    geom_line(size = 1.8, alpha = 0.6) + 
    ggtitle(paste("proportion antagonists = ", antprob)) +
    xlab("Time") + 
    ylab("Mean species trait (z)") +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 14), 
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size = 12))

  print(plotar)
}