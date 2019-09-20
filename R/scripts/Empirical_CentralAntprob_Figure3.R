# For empirical networks
# Code to test the influence of degree centrality of cheaters in
# trait MPD and variance of species in mutualistic networks.
# For that we gonna calculate and normalize the species degree centrality
# use the +1 SD species degree as central species and simulates the
# coevolution process with those species as cheaters

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)

# read all mutualism networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# create data.frame for final results
central_results = data.frame()

for(k in 1:length(redes)){
  print(k)
  
  for(a in 1:50){
    # Cheater centrality simulation
    M = as.matrix(redes[[k]])
    M[which(M > 1)] = 1
    M = SquareMatrix(M)
    n_sp = ncol(M)
    
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CentralAntagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    
    centralantagonize = CentralAntagonize(M)
    M = centralantagonize[[1]]
    V = centralantagonize[[2]]
    
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
    
    results = data.frame(net, rich, antprob, standev, mpd, c_ch)
    central_results = rbind(central_results, results)
    
    # Cheater Non-centrality simulation
    M = as.matrix(redes[[k]])
    M[which(M > 1)] = 1
    M = SquareMatrix(M)
    n_sp = ncol(M)
    
    antprob = centralantagonize[[3]]
    
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    
    empantagonize = Antagonize(M, antprob)
    M = empantagonize[[1]]
    V = empantagonize[[2]]
    
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
    
    net = names(redes[k])
    rich = as.numeric(ncol(M))
    standev = sd(z_mat[nrow(z_mat), ])
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    c_ch = "Distributed"
    
    results = data.frame(net, rich, antprob, standev, mpd, c_ch)
    central_results = rbind(central_results, results)
  
  }
}

#save(central_results, file = "central_results.RData")
load("central_results.RData")

plot_standev = ggplot(data = central_results) +
  geom_jitter(aes(x = as.factor(c_ch), y = standev, colour = rich), 
              position=position_jitter(0.2), alpha = 0.8) +
  ylab("Standart deviation of species traits") +
  xlab("") +
  labs(fill = "Richness") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

plot_mpd = ggplot(data = central_results) +
  geom_jitter(aes(x = as.factor(c_ch), y = standev, colour = rich), 
              position=position_jitter(0.2), alpha = 0.8) +
  ylab("MPD - Mean Pairwise Distance") +
  xlab("") +
  labs(fill = "Richness") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

ggsave(plot_standev, filename = "cheater_standev.png", dpi = 600,
       width = 10, height = 14, units = "cm")

ggsave(plot_mpd, filename = "cheater_mpd.png", dpi = 600,
       width = 10, height = 14, units = "cm")

