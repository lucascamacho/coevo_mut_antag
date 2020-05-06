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
 
      # get my results
      net = names(redes[k])
      rich = as.numeric(ncol(M))
      mpd = MeanPairDist(z_mat[nrow(z_mat), ])
      
      # insert these results in the data.frame
      results = data.frame(net, rich, antprob, mpd)
      p_data = rbind(p_data, results)
    }
  }
}

# save or load the data created
save(p_data, file = "antprob_var.RData")
save(list_mats, file = "~/Google Drive File Stream/Meu Drive/Trabalho/list_mats.RData")
load("antprob_var.RData")

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

# save or load the results
save(p_data, file = "antprob_var.RData")
load("antprob_var.RData")

# create and insert an mutualism type sequence
type = c(rep("Pollination", 24000), rep("Seed dispersal", 24000), rep("Ant-Plant", 24000))
p_data = cbind(p_data, type)

# subset the results to plot mutualism types
pol = p_data[which(p_data$type == "Pollination"), ]
seed = p_data[which(p_data$type == "Seed dispersal"), ]
ant = p_data[which(p_data$type == "Ant-Plant"), ]

# MPD Plots
# Pollination plot
pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = mpd), show.legend = FALSE, alpha = 0.05, 
             size = 0.5, color = "#D95F02") +
  geom_line(aes(x = antprob, y = mpd), stat = "smooth", method = "lm", 
            show.legend = FALSE) +
  ylim(0, 2) +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
       legend.text = element_text(size = 13))

# Seed dispersal plot
seed_plot = ggplot(data = seed) +
   geom_point(aes(x = antprob, y = mpd), show.legend = FALSE, alpha = 0.05, 
              size = 0.5, color = "#7570B3", size = 2) +
   geom_line(aes(x = antprob, y = mpd), stat = "smooth", method = "lm", 
             show.legend = FALSE) +
   ylim(0, 2) +
   xlab(" ") + ylab(" ") +
   scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title = element_text(size = 16), 
         legend.key.size = unit(0.9, "cm"),
         legend.text = element_text(size = 13))
 
# Ant-Plant plot
ant_plot = ggplot(data = ant) +
   geom_point(aes(x = antprob, y = mpd, size = 0.5), show.legend = FALSE, alpha = 0.05, 
              size = 0.5, color = "#1B9E77") +
   geom_line(aes(x = antprob, y = mpd), stat = "smooth", method = "lm", 
             show.legend = FALSE) +
   ylim(0, 4) +
   xlab(" ") + ylab(" ") +
   scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
   theme(axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title = element_text(size = 16), 
         legend.key.size = unit(0.9, "cm"),
         legend.text = element_text(size = 13))

# arrange the MPD plots in a line
plot_final = grid.arrange(ant_plot, pol_plot, seed_plot, nrow = 3)

# save the MPD plots 
ggsave(plot_final, filename = "antprob_cheater_mpd.pdf", dpi = 600,
        width = 12, height = 24, units = "cm",  bg = "transparent")

# Clusters plots
#Pollination plot
pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = opt_clusters), show.legend = FALSE, alpha = 0.05, 
             size = 0.5, color = "#D95F02") +
  geom_line(aes(x = antprob, y = opt_clusters), stat = "smooth", method = "loess", 
            show.legend = FALSE) +
  xlab(" ") + ylab(" ") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Seed dispersal plot
seed_plot = ggplot(data = seed) +
  geom_point(aes(x = antprob, y = opt_clusters), show.legend = FALSE, alpha = 0.05, 
             size = 0.5, color = "#7570B3") +
  geom_line(aes(x = antprob, y = opt_clusters), stat = "smooth", method = "loess", 
            show.legend = FALSE) +
  xlab(" ") + ylab(" ") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

#Ant-Plant plot
ant_plot = ggplot(data = ant) +
  geom_point(aes(x = antprob, y = opt_clusters), show.legend = FALSE, alpha = 0.05, 
             size = 0.5, color = "#1B9E77") +
  geom_line(aes(x = antprob, y = opt_clusters), stat = "smooth", method = "loess", 
            show.legend = FALSE) +
  xlab(" ") + ylab(" ") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# arrange Clusters plots in a line
plot_final = grid.arrange(ant_plot, pol_plot, seed_plot, nrow = 3)

#save Clusters plot
ggsave(plot_final, filename = "antprob_cheater_cluster.pdf", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")