# Code to test the influence of different epsilon values
#
# epsilon is our trait barrier for the exploitative interactions
# and here, we evaluate the influence of different values of epsilon.
#
# A single big network and a single analisys
#
# Vulnerability of results with different epsilons
# x = epsilon with (1, 2.5, 5, 7.5, 9)
# y = trait disparity and confidence intervals  
# 
# Vulnerability of results with with different epsilons
# with a certain variance of epsilon
# x = episilon com mean = 5, sd = (0.1, 1, 1.5)
# y = disparidade e intervalo de confiança
#
# load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)
library(reshape2)

camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  write(n_cl, file = "clustering.vec", append = TRUE)
  
  return(n_cl)
}

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

epsilon_vec = c(1, 2.5, 5, 7.5, 9)

p_data = data.frame()
list_mats = list()

for(a in 1:length(epsilon_vec)){
  for(q in 1:30){
    M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    
    antprob = 0.3 # get the antprob value
    
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
    epsilon = epsilon_vec[a]
    eq_dif = 0.0001
    t_max = 1000
    
    # simulate coevolution
    z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

    df = scale(t(z_mat))
    list_mats = list.append(list_mats, df)
    
    # get my results
    rich = as.numeric(ncol(M))
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    
    # insert these results in the data.frame
    results = data.frame(rich, antprob, epsilon, mpd)
    p_data = rbind(p_data, results)
  }
}

clusters_results = sapply(list_mats, camacho)
p_data = cbind(p_data, clusters_results)

save(p_data, file = "epsilon_control_sup.RData")
load("epsilon_control_sup.RData")

###

epsilon_sd = c(0.1, 1, 1.5)
p_data_sd = data.frame()
list_mats = list()

for(a in 1:length(epsilon_sd)){
  for(q in 1:30){
    M = as.matrix(redes[[6]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    
    antprob = 0.3 # get the antprob value
    
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
    epsilon = rnorm(1, mean = 5, sd = epsilon_sd[a])
    if(epsilon <= 0){
      epsilon = rnorm(1, mean = 5, sd = epsilon_sd[a])
    }
    eq_dif = 0.0001
    t_max = 1000
    
    # simulate coevolution
    z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
    
    df = scale(t(z_mat))
    list_mats = list.append(list_mats, df)
    
    # get my results
    rich = as.numeric(ncol(M))
    mpd = MeanPairDist(z_mat[nrow(z_mat), ])
    
    # insert these results in the data.frame
    results = data.frame(rich, antprob, epsilon, epsilon_sd[a], mpd)
    p_data_sd = rbind(p_data_sd, results)
  }
}

clusters_results = sapply(list_mats, camacho)
p_data_sd = cbind(p_data_sd, clusters_results)

save(p_data_sd, file = "epsilon_control_sample.RData")
load("epsilon_control_sample.RData")

### first plots from epsilon controlled
# load data
load("epsilon_control_sup.RData")
load("epsilon_control_sample.RData")

#reshape mpd data
new_data = p_data%>%
  group_by(epsilon)%>%
  summarise(mean_mpd = mean(mpd), mean_clust = mean(clusters_results))%>%
  as.data.frame()
modelo1 = lm(mean_mpd ~ epsilon, data = new_data)

plot_1 = ggplot(data = new_data, aes(x = epsilon, y = mean_mpd)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  ggtitle("B_NS-PP-vazquez_CerroLopez_2002, p = 0.3") +
  ylab("Average MPD") +
  xlab("Epsilon") +
  scale_x_continuous("Epsilon", labels = as.factor(new_data$epsilon), breaks = new_data$epsilon) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))

new_data = p_data%>%
  group_by(epsilon)%>%
  summarise(var_mpd = var(mpd), mean_clust = mean(clusters_results))%>%
  as.data.frame()
modelo2 = lm(var_mpd ~ epsilon, data = new_data)

plot_2 = ggplot(data = new_data, aes(x = epsilon, y = log(var_mpd))) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  ggtitle("B_NS-PP-vazquez_CerroLopez_2002, p = 0.3") +
  ylab("Log of Variance of MPD") +
  xlab("Epsilon") +
  scale_x_continuous("Epsilon", labels = as.factor(new_data$epsilon), breaks = new_data$epsilon) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))


plot_final = grid.arrange(plot_1, plot_2, nrow = 2)

ggsave(plot_final, filename = "Sup_epsilon.png", dpi = 600,
       width = 18, height = 15, units = "cm", bg = "transparent")

###
sub = p_data[,3:4]
plot_1 = ggplot(data = sub, aes(x = as.factor(epsilon), y = mean_mpd)) + 
  geom_point(alpha = 0.2) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red"
  ) +
  ylab("Mean Pairwise Distance") +
  xlab("Epsilon") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))

sub = p_data[,c(3,5)]
plot_2 = ggplot(data = sub, aes(x = as.factor(epsilon), y = clusters_results)) + 
  geom_point(alpha = 0.2) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red"
  ) +
  ylab("Number of trait clusters") +
  xlab("Epsilon") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))

plot_final = grid.arrange(plot_1, plot_2, nrow = 2)

ggsave(plot_final, filename = "Sup_epsilon.pdf", dpi = 600,
       width = 18, height = 15, units = "cm", bg = "transparent")

## second plots
sub = p_data_sd[,4:5]
plot_1 = ggplot(data = sub, aes(x = as.factor(epsilon_sd.a.), y = mpd)) + 
  geom_point(alpha = 0.5) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red"
  ) +
  ggtitle("B_NS-PP-vazquez_CerroLopez_2002, p = 0.3") +
  ylab("Mean Pairwise Distance") +
  xlab("Epsilon Standard Deviation") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))

sub =p_data_sd[,c(4,6)]
plot_2 = ggplot(data = sub, aes(x = as.factor(epsilon_sd.a.), y = clusters_results)) + 
  geom_point(alpha = 0.5) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red"
  ) +
  ylab("Number of trait clusters") +
  xlab("Epsilon Standard Deviation") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 10))

plot_final = grid.arrange(plot_1, plot_2, nrow = 2)

ggsave(plot_final, filename = "Sup_epsilon_SD.pdf", dpi = 600,
       width = 18, height = 15, units = "cm", bg = "transparent")