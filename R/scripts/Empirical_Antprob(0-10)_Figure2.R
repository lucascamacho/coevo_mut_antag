# Code to generate the Figure 2 of the paper where I
# show how the species trait standard deviation, MPD - 
# Mean Pairwise Distance and clustering number
# changes with the antprob (p)
# in different empirical networks.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the SD, MPD and 
# clustering number of species traits in the last timestep 
# of simulations

# load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

antprob_vec = seq(0.01, 1, 0.01)

p_data = data.frame()
list_mats = list()
list_degree = list()

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:length(antprob_vec)){ # loop to each antprob (p) value
    
    for(q in 1:30){ # 30 simulations per p value
      M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
      M[which(M > 1)] = 1 # if there are any error, correct that
      M = SquareMatrix(M) # square the adjancency matrix
      n_sp = ncol(M) # define the species number
      
      antprob = antprob_vec[a] # get the antprob value
      
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
      
      df = scale(t(z_mat))
      list_mats = list.append(list_mats, df)
      
      degree = SpDegree(M, V)
      list_degree = list.append(list_degree, degree)
      
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
}

# save or load the data created
#save(p_data, file = "p_data.RData")
#save(list_mats, file = "list_mats.RData")
#load("p_data.RData")

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
cl = makeCluster(detectCores() - 2)
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

# bind the clustering results with our data table
p_data = cbind(p_data, opt_clusters)

# use aggregate to group data by net and antprob
p_data = aggregate(p_data[ ,4:5], list(p_data$antprob, p_data$net), mean)
colnames(p_data)[1] = "antprob"
colnames(p_data)[2] = "net"

# create and insert an mutualism type sequence
type = c(rep("Pollination", 800), rep("Seed dispersal", 800), rep("Ant-Plant", 800))
p_data = cbind(p_data, type)

# using dplyr to get the averages and prepare the data frame
new_data = p_data%>%
  group_by(type, antprob)%>%
  summarise(mean_std = mean(standev), mean_mpd = mean(mpd), mean_clust = mean(opt_clusters))%>%
  as.data.frame()

plot_standev = ggplot(data = new_data) +
  geom_point(aes(x = antprob, y = mean_std, colour = type), show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean_std, colour = type), stat = "smooth", method = "lm",
            alpha = 0.8, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.1)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Standart deviation of species traits (σ)") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_mpd = ggplot(data = new_data) +
  geom_point(aes(x = antprob, y = mean_mpd, colour = type), show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean_mpd, colour = type), stat = "smooth", method = "lm",
            alpha = 0.8, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("MPD - Mean Pairwise Distance") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_clusters = ggplot(data = new_data) +
  geom_point(aes(x = antprob, y = mean_clust, colour = type), show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean_clust, colour = type), stat="smooth",method = "auto",
            alpha = 0.8, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0,0.1)) +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Average optimized number of species traits clusters") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(1.1, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  guides(color = guide_legend(title = "Mutualism type"))

ggsave(plot_mpd, filename = "antprob_cheater_mpd.png", dpi = 600,
       width = 14, height = 16, units = "cm",  bg = "transparent")

ggsave(plot_standev, filename = "antprob_cheater_sd.png", dpi = 600,
       width = 14, height = 16, units = "cm",  bg = "transparent")

ggsave(plot_clusters, filename = "antprob_clusters.png", dpi = 600,
       width = 16, height = 13, units = "cm", bg = "transparent")