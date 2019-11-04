# Run the coevolutionary model with interaction outcomes that are context-dependent.
# Based on a probability prob_change (Q), in each timestep of the simulation, an interaction
# outcome shifts (AA -> AM and AM -> MM).
#
# This script returns a simple graph with species traits changing in time due to coevolution.
# The asteriscs in the graph shows the timesteps in which the interactions shift occurs.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")

library(ggplot2)
library(reshape2)
library(cowplot)

antprob_vec = seq(0.001, 0.1, 0.001)
dados = data.frame()

for(i in 1:length(antprob_vec)){
  n_sp = 100
  M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  antprob = antprob_vec[i]
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob)
  M = antagonize[[1]]
  V = antagonize[[2]]
  
  s_init = rowSums(V)
  s = length(which(s_init > 1))
  c = mean(rowSums(V))

  results = data.frame(c, s)
  dados = rbind(dados, results)
}

ggplot(dados) +
  geom_point(aes(x = c, y = s))