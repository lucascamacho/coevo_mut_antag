# For empirical networks
# Code to test the influence of degree centrality of exploiters in
# trait MPD and trait clusters of species in mutualistic networks.
# For that we gonna calculate and normalize the species degree centrality
# use the +1 SD species degree as central species and simulates the
# coevolution process with those species as exploiters.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the MPDA and the
# number of species traits in the last timestep of simulations

# load packages and functions
setwd("E:/Lucas")

source("E:/Lucas/SquareMatrix.R")
source("E:/Lucas/MeanPairDist.R")

# simple function to apply the NbCluster in my results and save in clusterig.vec
# create an empty document to allocate the apply results
write("clustering", "clustering.vec")
camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  write(n_cl, file = "clustering.vec", append = TRUE)
  
  return(n_cl)
}

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
    source("E:/Lucas/CentralAntagonize.R")
    source("E:/Lucas/CoevoMutAntNet.R")
    
    # insert cheaters outcomes based on degree centrality
    centralantagonize = CentralAntagonize(M)
    M = centralantagonize[[1]]
    V = centralantagonize[[2]]
    
    source("E:/Lucas/Counting.R")
    
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
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    c_ch = "Central"
    
    # insert these results in our data.frame
    results = data.frame(net, rich, antprob, c_ch, mpd)
    central_results = rbind(central_results, results)
    
    # get traits to camacho function
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
    source("E:/Lucas/Antagonize.R")
    source("E:/Lucas/CoevoMutAntNet.R")
    
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
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    c_ch = "Distributed"
    
    # insert these results in our data.frame
    results = data.frame(net, rich, antprob, c_ch, mpd)
    central_results = rbind(central_results, results)
    
    #get traits to camacho function
    df = scale(t(z_mat))
    list_mats = list.append(list_mats, df)
  }
}

# save of load our results
save(central_results, file = "central_results.RData")
save(list_mats, file = "central_list_mats.RData")

load("central_results.RData")
load("central_list_mats.RData")

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

# create and insert and mutualism type sequence
type = c(rep("Pollination", 24000), rep("Seed dispersal", 24000), rep("Ant-Plant", 24000))
central_results = cbind(central_results, type)

# save or load our results
save(central_results, file = "central_results.RData")
#load("central_results.RData")



# read all mutualism networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# subset results for each mutualism type
pol = central_results[which(central_results$type == "Pollination"), ]
seed = central_results[which(central_results$type == "Seed dispersal"), ]
ant = central_results[which(central_results$type == "Ant-Plant"), ]

# using dplyr to get the averages and prepare the data frame
new_data = central_results%>%
  group_by(net, c_ch)%>%
  summarise(mean_clust = mean(opt_clusters), mean_mpd = mean(mpd))%>%
  as.data.frame()

# change data.frame to allocate the quantiles
new_data = cbind(new_data, min_mpd = NA, max_mpd = NA, min_clust = NA, max_clust = NA)

# prepare names to guide the quantile function
redes_vec = names(redes)
exp_vec = c("Central", "Distributed")
total = expand.grid(exp_vec, redes_vec)

# getting the quantie for each combination of interaction and scenario
for(i in 1:(length(redes)*2)){
  mut = total[i,2]
  exp = total[i,1]
  
  sub = subset(central_results, central_results$net == mut)
  sub = subset(sub, sub$c_ch == exp)
  
   new_data[i,6] = quantile(sub$mpd, probs = c(0.05, 0.95))[1]
   new_data[i,7] = quantile(sub$mpd, probs = c(0.05, 0.95))[2]
   
   new_data[i,8] = quantile(sub$opt_clusters, probs = c(0.05, 0.95))[1]
   new_data[i,9] = quantile(sub$opt_clusters, probs = c(0.05, 0.95))[2]
}

# MPD plots
# Pollination plot
 plot_mpd_pol = ggplot(data = new_data[c(1:16),], aes(x = as.factor(c_ch), 
                                        y = mean_mpd, group = net)) +
   geom_pointrange(aes(x = as.factor(c_ch), 
                       y = mean_mpd, group = net, ymin = min_mpd, ymax = max_mpd), 
                       show.legend = FALSE, alpha = 0.5, colour = "#D95F02") +
   geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#D95F02") +
   geom_line(show.legend = FALSE, colour = "#D95F02") +
   scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
   ylab("") +
   xlab("") +
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title = element_text(size = 16), 
         legend.key.size = unit(0.9, "cm"),
         legend.text = element_text(size = 13))

# Seed dispersal plot 
 plot_mpd_seed = ggplot(data = new_data[c(17:32),], aes(x = as.factor(c_ch), 
                                                  y = mean_mpd, group = net)) +
   geom_pointrange(aes(x = as.factor(c_ch), 
                       y = mean_mpd, group = net, ymin = min_mpd, ymax = max_mpd), 
                   show.legend = FALSE, alpha = 0.5, colour = "#7570B3") +
   geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#7570B3") +
   geom_line(show.legend = FALSE, colour = "#7570B3") +
   scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
   ylab("") +
   xlab("") +
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title = element_text(size = 16), 
         legend.key.size = unit(0.9, "cm"),
         legend.text = element_text(size = 13))

# Ant-Plant plot
 plot_mpd_ant = ggplot(data = new_data[c(33:48),], aes(x = as.factor(c_ch), 
                                                  y = mean_mpd, group = net)) +
   geom_pointrange(aes(x = as.factor(c_ch), 
                       y = mean_mpd, group = net, ymin = min_mpd, ymax = max_mpd), 
                   show.legend = FALSE, alpha = 0.5, colour = "#1B9E77") +
   geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#1B9E77") +
   geom_line(show.legend = FALSE, colour = "#1B9E77") +
   scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
   ylab("") +
   xlab("") +
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title = element_text(size = 16), 
         legend.key.size = unit(0.9, "cm"),
         legend.text = element_text(size = 13))
 
# arrange mpd plots in a line
plot_mpd = grid.arrange(plot_mpd_ant, plot_mpd_pol, plot_mpd_seed, nrow = 1)

# Clusters plots
# Pollintation plot
plot_cluster_pol = ggplot(data = new_data[c(1:16),], aes(x = as.factor(c_ch), 
                                                     y = mean_clust, group = net)) +
  geom_pointrange(aes(x = as.factor(c_ch), 
                      y = mean_clust, group = net, ymin = min_clust, ymax = max_clust), 
                  show.legend = FALSE, alpha = 0.5, colour = "#D95F02") +
  geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#D95F02") +
  geom_line(show.legend = FALSE, colour = "#D95F02") +
  scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))


# Seed dispersal plot
plot_cluster_seed = ggplot(data = new_data[c(17:32),], aes(x = as.factor(c_ch), 
                                                       y = mean_clust, group = net)) +
  geom_pointrange(aes(x = as.factor(c_ch), 
                      y = mean_clust, group = net, ymin = min_clust, ymax = max_clust), 
                  show.legend = FALSE, alpha = 0.5, colour = "#7570B3") +
  geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#7570B3") +
  geom_line(show.legend = FALSE, colour = "#7570B3") +
  scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Ant-plant plot
plot_cluster_ant = ggplot(data = new_data[c(33:48),], aes(x = as.factor(c_ch), 
                                                      y = mean_clust, group = net)) +
  geom_pointrange(aes(x = as.factor(c_ch), 
                      y = mean_clust, group = net, ymin = min_clust, ymax = max_clust), 
                  show.legend = FALSE, alpha = 0.5, colour = "#1B9E77") +
  geom_point(show.legend = FALSE, alpha = 0.05, size = 0.5, colour = "#1B9E77") +
  geom_line(show.legend = FALSE, colour = "#1B9E77") +
  scale_x_discrete(limits = rev(levels(new_data$c_ch)), labels = c("Random", "Central")) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# arrange the clusters plot and final plot
plot_total = grid.arrange(plot_cluster_ant, plot_cluster_pol, plot_cluster_seed, nrow = 1)
plot_total = grid.arrange(plot_mpd, plot_clust, nrow = 2)

# save the final plot
ggsave(plot_total, filename = "Sup_central_mpd_clust.png", dpi = 600,
       width = 18, height = 12, units = "cm", bg = "transparent")