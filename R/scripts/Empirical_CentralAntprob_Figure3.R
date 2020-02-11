# For empirical networks
# Code to test the influence of degree centrality of cheaters in
# trait MPD and standard deviation of species in mutualistic networks.
# For that we gonna calculate and normalize the species degree centrality
# use the +1 SD species degree as central species and simulates the
# coevolution process with those species as cheaters.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the SD, MPD and cluster
# number of species traits in the last timestep of simulations

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)

# read all mutualism networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# create data.frame for final results
central_results = data.frame()
list_mats = list()

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:1500){ #  loops to each matrix, 2 simulations per loop
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
    
    df = scale(t(z_mat))
    list_mats = list.append(list_mats, df)
    
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
    
    df = scale(t(z_mat))
    list_mats = list.append(list_mats, df)
  }
}

#save(central_results, file = "central_results.RData")
#save(list_mats, file = "~/Google Drive File Stream/Meu Drive/Trabalho/central_list_mats.RData")

# create an empty document to allocate the apply results
write("clustering", "clustering.vec")

# simple function to apply the NbCluster in my results and save in clusterig.vec
camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  write(n_cl, file = "clustering.vec", append = TRUE)
  
  return(n_cl)
}

# NbCluster using 14 computer cores
cl = makeCluster(detectCores())
clusterEvalQ(cl, {
  library(NbClust)
  camacho = function(list_mats){
    maximo = nrow(list_mats) - 1
    
    clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                       min.nc = 2, max.nc = maximo, method = "ward.D2", 
                       index = "gap")
    
    n_cl = clust_an$Best.nc[1]
    
    write(n_cl, file = "clustering.vec", append = TRUE)
    
    return(n_cl)
  }
})
opt_clusters = parallel::parSapply(cl, list_mats, camacho)
stopCluster(cl)

# bind the results with the original data table
central_results = cbind(central_results, opt_clusters)

#save(central_results, file = "central_results.RData")
#load("central_results.RData")

# save or load the results
central_results = aggregate(central_results[ ,4:5], 
                            list(central_results$c_ch, central_results$net), mean)
colnames(central_results)[1] = "c_ch"
colnames(central_results)[2] = "net"

# create and insert and mutualism type sequence
type = c(rep("Pollination", 16), rep("Seed dispersal", 16), rep("Ant-Plant", 16))
central_results = cbind(central_results, type)

# using dplyr to get the averages and prepare the data frame
new_data = central_results%>%
  group_by(type, c_ch)%>%
  summarise(mean_std = mean(standev), mean_mpd = mean(mpd), mean_clust = mean(opt_clusters))%>%
  as.data.frame()

# plot and save our results
plot_standev = ggplot(data = new_data, aes(x = as.factor(c_ch), 
                                           y = mean_std, colour = type, group = type)) +
  geom_point(aes(size = 0.1), show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(palette="Dark2") +
  ylab("Standard deviation of species traits (σ)") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_mpd = ggplot(data = new_data, aes(x = as.factor(c_ch), 
                                          y = mean_mpd, colour = type, group = type)) +
  geom_point(aes(size = 0.1), show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(palette="Dark2") +
  ylab("MPD - Mean Pairwise Distance") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_clust = ggplot(data = new_data, aes(x = as.factor(c_ch), 
                                         y = mean_clust, colour = type, group = type)) +
  geom_point(aes(size = 0.1), show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(palette="Dark2") +
  ylab("Average optimized number of species traits clusters") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))


ggsave(plot_standev, filename = "central_cheater_standev.png", dpi = 600,
       width = 14, height = 16, units = "cm", bg = "transparent")

ggsave(plot_mpd, filename = "central_cheater_mpd.png", dpi = 600,
       width = 14, height = 16, units = "cm", bg = "transparent")

ggsave(plot_clust, filename = "central_cheater_clust.png", dpi = 600,
       width = 14, height = 16, units = "cm", bg = "transparent")