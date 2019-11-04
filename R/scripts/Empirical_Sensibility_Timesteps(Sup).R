# This scripts is a sensibility analisys showing
# the maximum simulation timesteps from all the networks
# We gonna run 50 simulations for each network, count the
# max timestep and plot an histogram for each network

# change WD and load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

time_data = data.frame() # create data.frame to allocate results

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:100){ 
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    
    antprob = runif(1, 0, 1) # draw a antprob value from 0 to 1
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
    
    # insert cheaters exploitation outcomes
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
    
    net = names(redes[k])
    time_simu = nrow(z_mat)
    
    # insert these results in the data.frame
    results = data.frame(net, time_simu)
    time_data = rbind(time_data, results)
  }
}

plot = ggplot(time_data, aes(time_simu, group = net)) + 
  geom_density(show.legend = FALSE) +
  xlab("Length of simulations (in timesteps)") +
  ylab("Frequency") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 18), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))


ggsave(plot, filename = "Sensibility_Timesteps.png", dpi = 600,
                         width = 20, height = 14, units = "cm")