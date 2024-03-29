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
    
    for(q in 1:5){ # 30 simulations per p value
      M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
      M[which(M > 1)] = 1 # if there are any error, correct that
      M = SquareMatrix(M) # square the adjancency matrix
      n_sp = ncol(M) # define the species number
      
      antprob = antprob_vec[a] # get the antprob value
      
      # load functions
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EcoEvoCoevoMutAntag.R")
      
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
save(p_data, file = "Abnd_Cor_minor_antprob_var.RData")
#save(list_mats, file = "~/Desktop")
#load("Abnd_minor_antprob_var.RData")

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

# bind the clustering results with our data table
p_data = cbind(p_data, opt_clusters)

# create and insert an mutualism type sequence
type = c(rep("Pollination", 400), rep("Seed dispersal", 400), rep("Ant-Plant", 400))
p_data = cbind(p_data, type)

# save or load the results
save(p_data, file = "Abnd_Cor_cluster_minor_antprob_var.RData")
#load("minor_antprob_var.RData")

# MPD Plots
# subset the results to plot mutualism types
pol = p_data[which(p_data$type == "Pollination"), ]
seed = p_data[which(p_data$type == "Seed dispersal"), ]
ant = p_data[which(p_data$type == "Ant-Plant"), ]

# Ant
new_ant = ant%>%
  group_by(antprob)%>%
  summarise(mean_mpd = mean(mpd), up_mpd = quantile(mpd, probs = 0.95),
            down_mpd = quantile(mpd, probs = 0.05))%>%
  as.data.frame()

ant_plot = ggplot(data = new_ant) +
  geom_point(aes(x = antprob, y = mean_mpd), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#1B9E77") +
  geom_pointrange(aes(x = antprob, y = mean_mpd, ymin = down_mpd, ymax = up_mpd), width = .2, 
                  show.legend = FALSE, alpha = 0.5, color = "#1B9E77") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Pol
new_pol = pol%>%
  group_by(antprob)%>%
  summarise(mean_mpd = mean(mpd), up_mpd = quantile(mpd, probs = 0.95),
            down_mpd = quantile(mpd, probs = 0.05))%>%
  as.data.frame()

pol_plot = ggplot(data = new_pol) +
  geom_point(aes(x = antprob, y = mean_mpd), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#1B9E77") +
  geom_pointrange(aes(x = antprob, y = mean_mpd, ymin = down_mpd, ymax = up_mpd), width = .2, 
                  show.legend = FALSE, alpha = 0.5, color = "#D95F02") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Seed
new_seed = seed%>%
  group_by(antprob)%>%
  summarise(mean_mpd = mean(mpd), up_mpd = quantile(mpd, probs = 0.95),
            down_mpd = quantile(mpd, probs = 0.05))%>%
  as.data.frame()

seed_plot = ggplot(data = new_seed) +
  geom_point(aes(x = antprob, y = mean_mpd), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#7570B3") +
  geom_pointrange(aes(x = antprob, y = mean_mpd, ymin = down_mpd, ymax = up_mpd), width = .2, 
                  show.legend = FALSE, alpha = 0.5, color = "#7570B3") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_final = grid.arrange(ant_plot, pol_plot, seed_plot, nrow = 3)

ggsave(plot_final, filename = "Abnd_Cor_antprob_cheater_mpd_2.pdf", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")

pol = p_data[which(p_data$type == "Pollination"), ]
seed = p_data[which(p_data$type == "Seed dispersal"), ]
ant = p_data[which(p_data$type == "Ant-Plant"), ]

# Ant
new_ant = ant%>%
  group_by(antprob)%>%
  summarise(mean_clusters = mean(opt_clusters), up_clusters = quantile(opt_clusters, probs = 0.95),
            down_clusters = quantile(opt_clusters, probs = 0.05))%>%
  as.data.frame()

ant_plot = ggplot(data = new_ant) +
  geom_point(aes(x = antprob, y = mean_clusters), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#1B9E77") +
  #geom_pointrange(aes(x = antprob, y = mean_clusters, ymin = down_clusters, ymax = up_clusters), width = .2, 
  #                show.legend = FALSE, alpha = 0.5, color = "#1B9E77") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Pol
new_pol = pol%>%
  group_by(antprob)%>%
  summarise(mean_clusters = mean(opt_clusters), up_clusters = quantile(opt_clusters, probs = 0.95),
            down_clusters = quantile(opt_clusters, probs = 0.05))%>%
  as.data.frame()

pol_plot = ggplot(data = new_pol) +
  geom_point(aes(x = antprob, y = mean_clusters), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#D95F02") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Seed
new_seed = seed%>%
  group_by(antprob)%>%
  summarise(mean_clusters = mean(opt_clusters), up_clusters = quantile(opt_clusters, probs = 0.95),
            down_clusters = quantile(opt_clusters, probs = 0.05))%>%
  as.data.frame()

seed_plot = ggplot(data = new_seed) +
  geom_point(aes(x = antprob, y = mean_clusters), show.legend = FALSE, alpha = 0.5,
             size = 2, color = "#7570B3") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_final = grid.arrange(ant_plot, pol_plot, seed_plot, nrow = 3)

ggsave(plot_final, filename = "Abnd_Cor_antprob_cheater_clusters_2.pdf", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")

# Supplementary Figures
summary(lm(mpd ~ antprob, data = ant))

# MPD Plots
# Pollination plot
pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = mpd), show.legend = FALSE, alpha = 0.05,
             size = 0.5, color = "#D95F02") +
  geom_line(aes(x = antprob, y = mpd), stat = "smooth", method = "lm", 
            show.legend = FALSE) +
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

ggsave(plot_final, filename = "antprob_cheater_mpd.png", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")

# Clusters plots
#Pollination plot
pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = opt_clusters), show.legend = FALSE, alpha = 0.05, 
             size = 0.5, color = "#D95F02") +
  geom_line(aes(x = antprob, y = opt_clusters), stat = "smooth", method = "loess", 
            show.legend = FALSE) +
  xlab(" ") + ylab(" ") +
  scale_y_continuous(limits = c(0,6)) +
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
  scale_y_continuous(limits = c(0,6)) +
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
  scale_y_continuous(limits = c(0,7)) +
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

ggsave(plot_final, filename = "antprob_cheater_cluster.png", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")