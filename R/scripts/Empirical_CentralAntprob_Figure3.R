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
library(gridExtra)
library(reshape2)
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

load("central_results.RData")
load("~/Google Drive File Stream/Meu Drive/Trabalho/central_list_mats.RData")

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
cl = makeCluster(detectCores() - 4)
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
load("central_results.RData")

# create and insert and mutualism type sequence
type = c(rep("Pollination", 24000), rep("Seed dispersal", 24000), rep("Ant-Plant", 24000))
central_results = cbind(central_results, type)

pol = central_results[which(central_results$type == "Pollination"), ]
seed = central_results[which(central_results$type == "Seed dispersal"), ]
ant = central_results[which(central_results$type == "Ant-Plant"), ]

# using dplyr to get the averages and prepare the data frame
new_data = central_results%>%
  group_by(type, c_ch)%>%
  summarise(mean_std = mean(standev), mean_mpd = mean(mpd), mean_clust = mean(opt_clusters))%>%
  as.data.frame()

new_data = cbind(new_data, min_mpd = NA, max_mpd = NA, min_clust = NA, max_clust = NA)

mut_vec = c("Ant-Plant", "Pollination", "Seed dispersal")
exp_vec = c("Central", "Distributed")

total = expand.grid(exp_vec, mut_vec)

for(i in 1:nrow(total)){
  mut = total[i,2]
  exp = total[i,1]
  
  sub = subset(central_results, central_results$type == mut)
  sub = subset(sub, sub$c_ch == exp)
  
  media_mpd = vector()
  media_clust = vector()
  
  for(j in 1:1000){
    novo = sample(sub$mpd, 12000, replace = TRUE)
    media_mpd[j] = mean(novo)
    
    novo = sample(sub$opt_clusters, 12000, replace = TRUE)
    media_clust[j] = mean(novo)
    }
  
  new_data[i,6] = quantile(media_mpd, probs = c(0.05, 0.95))[1]
  new_data[i,7] = quantile(media_mpd, probs = c(0.05, 0.95))[2]
  
  new_data[i,8] = quantile(media_clust, probs = c(0.05, 0.95))[1]
  new_data[i,9] = quantile(media_clust, probs = c(0.05, 0.95))[2]
  
}

plot_mpd = ggplot(data = new_data, aes(x = as.factor(c_ch), 
                                       y = mean_mpd, colour = type, group = type)) +
  geom_pointrange(aes(x = as.factor(c_ch), 
                      y = mean_mpd, colour = type, group = type, ymin = min_mpd, ymax = max_mpd), 
                      show.legend = FALSE) +
  geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, color = "#1B9E77") +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(palette="Dark2") +
  scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
  ylab("Mean Pairwise Distance between species traits") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_clust = ggplot(data = new_data, aes(x = as.factor(c_ch), 
                                         y = mean_clust, colour = type, group = type)) +
  geom_pointrange(aes(x = as.factor(c_ch), 
                      y = mean_clust, colour = type, group = type, ymin = min_clust, ymax = max_clust), 
                  show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(palette="Dark2") +
  scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
  ylab("Average optimized number of species traits clusters") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

ggsave(plot_mpd, filename = "central_cheater_mpd.png", dpi = 600,
       width = 14, height = 16, units = "cm", bg = "transparent")

ggsave(plot_clust, filename = "central_cheater_clust.png", dpi = 600,
       width = 14, height = 16, units = "cm", bg = "transparent")