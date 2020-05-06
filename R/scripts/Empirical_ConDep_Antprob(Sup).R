# Code to make one of the figures of Sup. Materials
# In this figure, I show the influence of interaction shifts
# during the coevolutionary process. The interaction shift happend
# in time don't change our main result of exploitation
# causing higher trait disparity in mutualistic networks.
#
# For each of our 24 empirical networks we gonna run 
# 100 simulations and calculate the MPD and the number of
# species trait clusters in the last timestep of simulations.
#
# In each simulation we have a new probability G that an 
# positive effect pass to negative. Also, one negative effect 
# will pass to positive.
#
# Here we are using the biggest of our networks
# rede[6] = "B_NS-PP-vazquez_CerroLopez_2002"

# set WD
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

# load packages and functions
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  return(n_cl)
}

library(ggplot2)
library(cowplot)
library(NbClust)
library(gridExtra)

# read all mutualism networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# define the vector of p and g and combine 2-by-2 these values
antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0.01, 0.1, 0.5)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 100), ]

# create sheet to alocate the simulations results
condep_data = data.frame()

for(a in 1:nrow(combs)){
  print(a)
  
  M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
  M[which(M > 1)] = 1 # if there are any error, correct that
  M = SquareMatrix(M) # square the adjancency matrix
  n_sp = ncol(M) # define the species number
  
  # define p  and g values
  antprob = combs[a,1]
  prob_change = combs[a,2]
   
  # load functions
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")
  
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
  
  #calculate number of species trait clusters
  df = scale(t(traits))
  opt_clust = camacho(df)
  
  # calculate MPD
  mpd = MeanPairDist(traits[nrow(traits), ])
  
  # allocate results in data.frame
  results = data.frame(antprob, prob_change, mpd, opt_clust)
  condep_data = rbind(condep_data, results) 
}

#save or load ou data
save(condep_data, file = "sup_condep_data.RData")
#load("sup_condep_data.RData")

# plot the results
plot_mpd = ggplot(data = condep_data) +
  geom_violin(aes(x = as.character(prob_change), y = mpd), fill = "grey80") +
  geom_point(aes(x = as.character(prob_change), y = mpd), size = 0.4, 
             shape = 21, alpha = 0.4) +
  facet_grid(prob_change~antprob) +
  theme_bw(base_size = 17) +
  labs(y = "Mean Pairwise Distance between species traits", 
       x = "")

plot_clust = ggplot(data = condep_data) +
  geom_violin(aes(x = as.character(prob_change), y = opt_clust), fill = "grey80") +
  geom_point(aes(x = as.character(prob_change), y = opt_clust), size = 0.4, 
             shape = 21, alpha = 0.4) +
  facet_grid(prob_change~antprob) +
  theme_bw(base_size = 17) +
  labs(y = "Optimized number of species traits cluster", 
       x = "")

# arrange the plots and save the final plot
plot_total = grid.arrange(plot_mpd, plot_clust, nrow = 2)

#ggsave(plot_total, filename = "Sensibility_ConDep.png", dpi = 600,
#       width = 21, height = 29, units = "cm")

#ggsave(plot_total, filename = "Sensibility_ConDep.pdf", dpi = 600,
#       width = 21, height = 29, units = "cm")