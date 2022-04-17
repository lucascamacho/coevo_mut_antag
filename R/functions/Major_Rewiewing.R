# Major Review
# Permitir rewiewring nas interacoes cortadas com b_ij
#
# Code of figure 4, but inside a for to test different barrier values
# Figure 4 of the paper.
#
# Code to explore the influence of exploitative interactions in the structure
# of the adjancency matrix of interactions. We are considering a second
# trait barrier (b) that will cut interactions off if this barrier is tranpassed
# This will allow us to see what's the effect of exploitation in coevolution and
# how this change's the structure of ecological interacions in the community.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the nestedness and modularity
# of out adjancency matrix of interactions.

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(gridExtra)
library(bipartite)
library(igraph)
library(dplyr)
library(schoolmath)

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

# create data.frame to store all my results
antprob_vec = seq(0.1, 1, 0.1)

final_fl = data.frame()

for(k in 1:length(redes)){ # loop to each matrix of interactions
  print(k)
  
  for(a in 1:length(antprob_vec)){ # 100 loops to each matrix
    
    for(q in 1:10){ 
      M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
      M[which(M > 1)] = 1 # if there are any error, correct that
      M = SquareMatrix(M) # square the adjancency matrix
      n_sp = ncol(M) # define the species number
      
      # sample an antprob value
      antprob = antprob_vec[a]
      
      # load functions
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/DiscoCoevoMutAntNet.R")
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
      
      # antagonize matrices (transform positive in negative effects)
      antagonize = Antagonize(M, antprob)
      M = antagonize[[1]]
      V = antagonize[[2]]
      
      # DesQuad matrices
      li = 1:dim(redes[[k]])[1]
      co_1 = dim(redes[[k]])[1] + 1
      co_2 = dim(redes[[k]])[1] + dim(redes[[k]])[2]
      
      # init_m and W are non squared adjacency matrices
      init_m = as.matrix(M + V)[li, co_1:co_2]
      W = as.matrix(M + V)[li, co_1:co_2]
      
      # coevolutionary model parameters
      phi = 0.2
      alpha = 0.2
      theta = runif(n_sp, 0, 10)
      init = runif(n_sp, 0, 10)
      p = 0.1
      epsilon = 5
      eq_dif = 0.000001
      t_max = 1000
      bar = 7
      
      # simulate coevolution
      simulation = DiscoCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max, bar)
      z_mat = simulation[[1]] # trait value matrix
      matrices = simulation[[2]] # all adjacency matrices from coevo function
      last_m = matrices[[length(matrices)]] # last adjancency matrix of the coevo process 
      last_m = as.matrix(last_m)[li, co_1:co_2] # DeSquad the last adjancency matrix
      
      net = paste(names(redes[k]), "_", a, sep = "") # name of the current empirical network
      rich = as.numeric(ncol(M)) # richness
      
      # run the control process
      dif_int = sum(init_m) - sum(last_m) # how much interactions coevolution breaks? 
      index = which(init_m == 1) # where are the interactions in the adjancency matrix?
      if(is.positive(dif_int) == TRUE){
        ints = sample(index, dif_int) # sample this interactions and break the same number
        W[ints] = 0 # of interactions that coevolution breaks
      }
      
      # get nestedness measures
      dnest_control = nested(W, method = "NODF2") - nested(init_m, method = "NODF2")
      dnest_coevo = nested(last_m, method = "NODF2") - nested(init_m, method = "NODF2")
      
      # get modularity measures
      dmod_control = modularity(cluster_louvain(graph_from_incidence_matrix(W))) - 
        modularity(cluster_louvain(graph_from_incidence_matrix(init_m)))
      dmod_coevo = modularity(cluster_louvain(graph_from_incidence_matrix(last_m))) - 
        modularity(cluster_louvain(graph_from_incidence_matrix(init_m)))
      
      # get all the results
      results = data.frame(net, rich, antprob, bar, dnest_control, dnest_coevo, dmod_control, dmod_coevo)
      final_fl = rbind(final_fl, results) # put results in data.frame
    }
  }
}

# save or load the RData file
save(final_fl, file = "Major_Rewiewing.RData")
#load("Major_Rewiewing.RData")

new_data = final_fl%>%
  group_by(antprob)%>%
  summarise(mean_nest_control = mean(dnest_control),
            mean_nest_coevo = mean(dnest_coevo),
            up_nest_coevo = quantile(dnest_coevo, probs = 0.95),
            down_nest_coevo = quantile(dnest_coevo, probs = 0.05),
            mean_mod_control = mean(dmod_control),
            mean_mod_coevo = mean(dmod_coevo),
            up_mod_coevo = quantile(dmod_coevo, probs = 0.95),
            down_mod_coevo = quantile(dmod_coevo, probs = 0.05))%>%
  as.data.frame()

plot_nest = ggplot(data = new_data, aes(x = antprob, y = mean_nest_coevo)) +
  geom_point(alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Delta nestedness with Rewiewing") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_mod = ggplot(data = new_data, aes(x = antprob, y = mean_mod_coevo)) +
  geom_point(alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Delta modularity with Rewiewing") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_final = grid.arrange(plot_mod, plot_nest, nrow = 2)

ggsave(plot_final, filename = "Major_Rewiewing.png", dpi = 600,
       width = 15, height = 15, units = "cm",  bg = "transparent")