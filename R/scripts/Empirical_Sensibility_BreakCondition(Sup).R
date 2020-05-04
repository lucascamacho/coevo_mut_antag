# Sensiblity analisys of the biggest network that we are using
# to test differences in the break condition from simulations.
# The break condition says that if the difference between traits is
# smaller than 10^4, the simulations stops. We gonna test if the
# trait matching metrics that we use changes with vlues of 10^6 and 10^8.

# change WD and load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

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

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

p_data = data.frame() # create data.frame to allocate results

for(a in 1:50){ 
  # 10^4
  print(a)
  M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
  M[which(M > 1)] = 1 # if there are any error, correct that
  M = SquareMatrix(M) # square the adjancency matrix
  n_sp = ncol(M) # define the species number
  
  antprob = 0.8 # draw a antprob value from 0 to 1
  
  # load functions
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")

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
  
  df = scale(t(z_mat))
  
  opt_clust = camacho(df)
  
  # get my results
  mpd = MeanPairDist(z_mat[nrow(z_mat), ])
  
  # insert these results in the data.frame
  results = data.frame(eq_dif, mpd, opt_clust)
  p_data = rbind(p_data, results)
  
  # 10^6
  M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
  M[which(M > 1)] = 1 # if there are any error, correct that
  M = SquareMatrix(M) # square the adjancency matrix
  n_sp = ncol(M) # define the species number
  
  antprob = 0.8 # draw a antprob value from 0 to 1
  
  # load functions
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
  
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
  eq_dif = 0.000001
  t_max = 1000
  
  # simulate coevolution
  z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  
  df = scale(t(z_mat))
  
  opt_clust = camacho(df)
  
  # get my results
  mpd = MeanPairDist(z_mat[nrow(z_mat), ])
  
  # insert these results in the data.frame
  results = data.frame(eq_dif, mpd, opt_clust)
  p_data = rbind(p_data, results)
  
  # 10^8
  M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
  M[which(M > 1)] = 1 # if there are any error, correct that
  M = SquareMatrix(M) # square the adjancency matrix
  n_sp = ncol(M) # define the species number
  
  antprob = 0.8 # draw a antprob value from 0 to 1
  
  # load functions
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
  
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
  eq_dif = 0.00000001
  t_max = 1000
  
  # simulate coevolution
  z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
  
  df = scale(t(z_mat))
  
  opt_clust = camacho(df)
  
  # get my results
  mpd = MeanPairDist(z_mat[nrow(z_mat), ])
  
  # insert these results in the data.frame
  results = data.frame(eq_dif, mpd, opt_clust)
  p_data = rbind(p_data, results)
}

#save(p_data, file = "sup_break_data.RData")
load("sup_break_data.RData")

plot_mpd = ggplot(data = p_data) +
  geom_point(aes(x = as.factor(eq_dif), y = mpd), 
             position=position_jitter(h=0, w=0.2), fill = "grey90", alpha = 0.4) +
  xlab("") +
  ylab("Mean Pairwise Distance between traits") +
  scale_y_continuous(limits = c(0.5, 0.65)) +
  scale_x_discrete(limits = rev(levels(as.factor(p_data$eq_dif))), labels = c("1e-4", "1e-6", "1e-8")) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 13))

plot_clust = ggplot(data = p_data) +
  geom_point(aes(x = as.factor(eq_dif), y = opt_clust), 
             position=position_jitter(h=0, w=0.2), fill = "grey90", alpha = 0.4) +
  xlab("Minimum trait difference between generations to stop simulation") +
  ylab("Number of species trait clusters") +
  scale_y_continuous(limits = c(2, 3)) +
  scale_x_discrete(limits = rev(levels(as.factor(p_data$eq_dif))), labels = c("1e-4", "1e-6", "1e-8")) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 13))

plot_total = grid.arrange(plot_mpd, plot_clust, nrow = 2)

#ggsave(plot_total, filename = "Sup_Figura_5.pdf", dpi = 600,
#       width = 21, height = 29, units = "cm")

ggsave(plot_total, filename = "Sup_Figura_5.png", dpi = 600,
       width = 21, height = 29, units = "cm")