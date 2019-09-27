# For empirical networks
# Code to test the influence of degree centrality of cheaters in
# trait MPD and variance of species in mutualistic networks.
# For that we gonna calculate and normalize the species degree centrality
# use the +1 SD species degree as central species and simulates the
# coevolution process with those species as cheaters.
#
# For each of our 24 empirical networks we gonna run 
# 100 simulations and calculate the SD and MPD of
# species traits in the last timestep of simulations

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)
library(viridis)

# read all mutualism networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# create data.frame for final results
central_results = data.frame()

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:50){ # 50 loops to each matrix
    # Cheater centrality simulation
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CentralAntagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    
    # insert cheaters outcomes based on degree centrality
    centralantagonize = CentralAntagonize(M)
    M = centralantagonize[[1]]
    V = centralantagonize[[2]]
    
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
    
    # counting interactions AA, AM and MM (AA must be zero)
    c = Counting(M, V)
    
    # if are some interference AA, stop the loop
    if(c[[1]] != 0){
      break("AA generated")
    }
    
    # coevolutionary model parameters
    phi = 0.2
    alpha = 0.2
    theta = runif(n_sp, 0, 10)
    init = runif(n_sp, 0, 10)
    p = 0.1
    epsilon = 5
    eq_dif = 0.0001
    t_max = 1000
    
    # simulate coevolution
    z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
    
    # get results
    net = names(redes[k])
    rich = as.numeric(ncol(M))
    antprob = centralantagonize[[3]]
    standev = sd(z_mat[nrow(z_mat), ])
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    c_ch = "Central"
    
    # insert these results in our data.frame
    results = data.frame(net, rich, antprob, standev, mpd, c_ch)
    central_results = rbind(central_results, results)
    
    # Cheater Non-centrality simulation
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    
    # use the same frequency of cheaters calculated above
    antprob = centralantagonize[[3]]
    
    #load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    
    # insert cheaters outcomes in the network
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
    
    # simulate coevolution
    z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
    
    # get our results
    net = names(redes[k])
    rich = as.numeric(ncol(M))
    standev = sd(z_mat[nrow(z_mat), ])
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    c_ch = "Distributed"
    
    # insert these results in our data.frame
    results = data.frame(net, rich, antprob, standev, mpd, c_ch)
    central_results = rbind(central_results, results)
  
  }
}

# save or load the results
#save(central_results, file = "central_results.RData")
load("central_results.RData")

# plot and save our results
plot_standev = ggplot(data = central_results) +
  geom_boxplot(aes(x = as.factor(c_ch), y = standev), fill = "grey90") +
  ylab("Standard deviation of species traits (σ)") +
  xlab("") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

plot_mpd = ggplot(data = central_results) +
  geom_boxplot(aes(x = as.factor(c_ch), y = mpd), fill = "grey90") +
  ylab("MPD - Mean Pairwise Distance") +
  xlab("") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

ggsave(plot_standev, filename = "central_cheater_standev.png", dpi = 600,
       width = 10, height = 14, units = "cm")

ggsave(plot_mpd, filename = "central_cheater_mpd.png", dpi = 600,
       width = 10, height = 14, units = "cm")