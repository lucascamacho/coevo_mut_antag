# Sensibility analyses for epsilon
#
# Code to generate the Figure 2 of the paper where I
# show how the MPD - Mean Pairwise Distance and 
# trait clustering number changes with the antprob (p)
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
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

antprob_vec = seq(0.1, 1, 0.1)

p_data = data.frame()
list_mats = list()
list_degree = list()

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:length(antprob_vec)){ # loop to each antprob (p) value
    
    for(q in 1:10){ # 10 simulations per e seq
      M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
      M[which(M > 1)] = 1 # if there are any error, correct that
      M = SquareMatrix(M) # square the adjancency matrix
      n_sp = ncol(M) # define the species number
      
      antprob = antprob_vec[a] # get the antprob value
      
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
      seq_p = seq(1, 10, 1) # sequence of p values
      p = 0.1 #chose the p value
      epsilon_seq = seq(1, 10, 1)
      epsilon = rep(epsilon_seq[q], n_sp)
      eq_dif = 0.0001
      t_max = 1000
      
      # simulate coevolution
      z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
      
      df = scale(t(z_mat))
      list_mats = list.append(list_mats, df)
      
      # get my results
      net = names(redes[k])
      rich = as.numeric(ncol(M))
      mpd = MeanPairDist(z_mat[nrow(z_mat), ])
      
      # insert these results in the data.frame
      epsilon = unique(epsilon)
      results = data.frame(net, rich, antprob, epsilon, mpd)
      p_data = rbind(p_data, results)
    }
  }
}

# save or load the data created
save(p_data, file = "Major_sensibility_e.RData")
#save(list_mats, file = "~/Google Drive File Stream/Meu Drive/Trabalho/list_mats.RData")
#load("antprob_var.RData")

# NbCluster using 14 computer cores
cl = makeCluster(detectCores()-1)
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

# save or load the results
#save(p_data, file = "Major_sensibility_e.RData")
load("Major_sensibility_e.RData")

plot(p_data$antprob, p_data$opt_clusters)

new_data = p_data%>%
  group_by(antprob, epsilon)%>%
  summarise(mean_mpd = mean(mpd),
            mean_clusters = mean(opt_clusters))%>%
  as.data.frame()

plot_mpd = ggplot(data = new_data, aes(x = antprob, y = mean_mpd)) +
  geom_point(aes(color = epsilon, fill = epsilon), alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Average MPD") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_clusters = ggplot(data = new_data, aes(x = antprob, y = mean_clusters)) +
  geom_point(aes(color = epsilon, fill = epsilon), alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Averge Number of Clusters of Species Traits") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_final = grid.arrange(plot_mpd, plot_clusters, nrow = 2)

ggsave(plot_final, filename = "Major_Sensibility_epsilon.png", dpi = 600,
       width = 15, height = 15, units = "cm",  bg = "transparent")
