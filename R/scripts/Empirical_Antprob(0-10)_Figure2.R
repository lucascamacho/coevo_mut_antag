# Code to generate the Figure 2 of the paper where I
# show how the species trait standard deviation and MPD - 
# Mean Pairwise Distance changes with the antprob (p)
# in different empirical networks.
#
# For each of our 24 empirical networks we gonna run 
# 100 simulations and calculate the SD and MPD of
# species traits in the last timestep of simulations

# change WD and load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

p_data = data.frame() # create data.frame to allocate results

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:100){ # 100 simulations per empirical matrix
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number

    antprob = runif(1, 0, 1) # draw a antprob value from 0 to 1
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
    
    # insert cheaters exploitation outcomes
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
    
    # get my results
    net = names(redes[k])
    rich = as.numeric(ncol(M))
    standev = sd(z_mat[nrow(z_mat), ])
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    
    # insert these results in the data.frame
    results = data.frame(net, rich, antprob, standev, mpd)
    p_data = rbind(p_data, results)
  }
}

# save or load the data created
save(p_data, file = "antprob_var.RData")
load(file = "antprob_var.RData")

# plot and save the results using ggplot2
plot_standev = ggplot(data = p_data) +
  geom_point(aes(x = antprob, y = standev, colour = rich), alpha = 0.8) +
  geom_smooth(aes(x = antprob, y = standev), colour = "red") +
  #  geom_hline(yintercept = 0) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Standart deviation of species trait") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

plot_mpd = ggplot(data = p_data) +
  geom_point(aes(x = antprob, y = mpd, colour = rich), alpha = 0.8) +
  geom_smooth(aes(x = antprob, y = mpd), colour = "red") +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("MPD - Mean Pairwise Distance") +
  labs(fill = "Richness") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

#ggsave(plot_standev, filename = "antprob_cheater_sd.png", dpi = 600,
#       width = 20, height = 14, units = "cm")

#ggsave(plot_mpd, filename = "antprob_cheater_mpd.png", dpi = 600,
#       width = 20, height = 14, units = "cm")