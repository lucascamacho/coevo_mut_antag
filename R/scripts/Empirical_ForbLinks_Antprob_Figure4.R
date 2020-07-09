# Figure 4 of the paper.
#
# Code to explore the influence of exploitative interactions in the structure
# of the adjancency matrix of interactions. We are considering a second
# trait barrier (b) that will cut interactions off if this barrier is tranpassed
# This will allow us to see what's the effect of exploitation in coevolution and
# how this change's the structure of ecological interacions in the community.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the nestedness and modularity
# of out adjancency matrix of interactions.

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(gridExtra)
library(bipartite)
library(igraph)
library(dplyr)

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

# create data.frame to store all my results
antprob_vec = seq(0.01, 1, 0.01)
final_fl = data.frame()

for(k in 1:length(redes)){ # loop to each matrix of interactions
  print(k)
  
  for(a in 1:length(antprob_vec)){ # 100 loops to each matrix
    
    for(q in 1:30){ # 30 simulations per p value
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number

    # sample an antprob value
    antprob = antprob_vec[a]
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/DiscoCoevoMutAntNet.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
    
    # antagonize matrices (transform positive in negative effects)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
    
    # DesQuad matrices
    li = 1:dim(redes[[k]])[1]
    co_1 = dim(redes[[k]])[1] + 1
    co_2 = dim(redes[[k]])[1] + dim(redes[[k]])[2]
    
    # init_m and W are non squared adjacency matrices
    init_m = as.matrix(M + V)[li, co_1:co_2]
    W = as.matrix(M + V)[li, co_1:co_2]
    
    # coevolutionary model parameters
    phi = 0.2
    alpha = 0.2
    theta = runif(n_sp, 0, 10)
    init = runif(n_sp, 0, 10)
    p = 0.1
    epsilon = 5
    eq_dif = 0.000001
    t_max = 1000
    bar = 7
    
    # simulate coevolution
    simulation = DiscoCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max, bar)
    z_mat = simulation[[1]] # trait value matrix
    matrices = simulation[[2]] # all adjacency matrices from coevo function
    last_m = matrices[[length(matrices)]] # last adjancency matrix of the coevo process 
    last_m = as.matrix(last_m)[li, co_1:co_2] # DeSquad the last adjancency matrix
    
    net = paste(names(redes[k]), "_", a, sep = "") # name of the current empirical network
    rich = as.numeric(ncol(M)) # richness
    
    # run the control process
    dif_int = sum(init_m) - sum(last_m) # how much interactions coevolution breaks? 
    index = which(init_m == 1) # where are the interactions in the adjancency matrix?
    ints = sample(index, dif_int) # sample this interactions and break the same number
    W[ints] = 0 # of interactions that coevolution breaks
   
    # get nestedness measures
    dnest_control = nested(W, method = "NODF2") - nested(init_m, method = "NODF2")
    dnest_coevo = nested(last_m, method = "NODF2") - nested(init_m, method = "NODF2")
    
    # get modularity measures
    dmod_control = modularity(cluster_louvain(graph_from_incidence_matrix(W))) - 
      modularity(cluster_louvain(graph_from_incidence_matrix(init_m)))
    dmod_coevo = modularity(cluster_louvain(graph_from_incidence_matrix(last_m))) - 
      modularity(cluster_louvain(graph_from_incidence_matrix(init_m)))
    
    # get all the results
    results = data.frame(net, rich, antprob, dnest_control, dnest_coevo, dmod_control, dmod_coevo)
    final_fl = rbind(final_fl, results) # put results in data.frame
    }
  }
}

# save or load the RData file
save(final_fl, file = "data_structure.RData")
load("data_structure.RData")

# create and insert an mutualism type sequence
type = c(rep("Pollination", 24000), rep("Seed dispersal", 24000), rep("Ant-Plant", 24000))
final_fl = cbind(final_fl, type)

# subset for plotting
pol = final_fl[which(final_fl$type == "Pollination"), ]
seed = final_fl[which(final_fl$type == "Seed dispersal"), ]
ant = final_fl[which(final_fl$type == "Ant-Plant"), ]

new_data = final_fl%>%
  group_by(type)%>%
  summarise(mean_nest_control = mean(dnest_control),
            mean_nest_coevo = mean(dnest_coevo),
            mean_mod_control = mean(dmod_control),
            mean_mod_coevo = mean(dmod_coevo))%>%
  as.data.frame()

sd_new_data = final_fl%>%
  group_by(type)%>%
  summarise(sd_nest_control = sd(dnest_control),
            sd_nest_coevo = sd(dnest_coevo),
            sd_mod_control = sd(dmod_control),
            sd_mod_coevo = sd(dmod_coevo))%>%
  as.data.frame()

# Pollintation plots
mod_pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = dmod_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#D95F02") +
  geom_point(aes(x = antprob, y = dmod_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dmod_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#D95F02") +
#  geom_line(aes(x = antprob, y = dmod_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

nest_pol_plot = ggplot(data = pol) +
  geom_point(aes(x = antprob, y = dnest_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#D95F02") +
  geom_point(aes(x = antprob, y = dnest_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dnest_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#D95F02") +
#  geom_line(aes(x = antprob, y = dnest_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Seed dispersal plots plots
mod_seed_plot = ggplot(data = seed) +
  geom_point(aes(x = antprob, y = dmod_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#7570B3") +
  geom_point(aes(x = antprob, y = dmod_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dmod_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#7570B3") +
#  geom_line(aes(x = antprob, y = dmod_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

nest_seed_plot = ggplot(data = seed) +
  geom_point(aes(x = antprob, y = dnest_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#7570B3") +
  geom_point(aes(x = antprob, y = dnest_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dnest_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#7570B3") +
#  geom_line(aes(x = antprob, y = dnest_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# Ant-plant plots
mod_ant_plot = ggplot(data = ant) +
  geom_point(aes(x = antprob, y = dmod_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#1B9E77") +
  geom_point(aes(x = antprob, y = dmod_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dmod_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#1B9E77") +
#  geom_line(aes(x = antprob, y = dmod_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

nest_ant_plot = ggplot(data = ant) +
  geom_point(aes(x = antprob, y = dnest_coevo), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "#1B9E77") +
  geom_point(aes(x = antprob, y = dnest_control), show.legend = FALSE, alpha = 0.1, 
             size = 0.5, color = "black") +
#  geom_line(aes(x = antprob, y = dnest_coevo), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "#1B9E77") +
#  geom_line(aes(x = antprob, y = dnest_control), stat = "smooth", method = "loess", 
#            show.legend = FALSE, color = "black") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# arrange all the plots
plot_pol = grid.arrange(mod_pol_plot, nest_pol_plot, nrow = 2)
plot_seed = grid.arrange(mod_seed_plot, nest_seed_plot, nrow = 2)
plot_ant = grid.arrange(mod_ant_plot, nest_ant_plot, nrow = 2)

plot_total = grid.arrange(plot_ant, plot_pol, plot_seed, nrow = 1)

# save the final plot
ggsave(plot_total, filename = "antprob_structure.pdf", dpi = 600,
       width = 25, height = 12, units = "cm",  bg = "transparent")
ggsave(plot_total, filename = "antprob_structure.png", dpi = 600,
       width = 25, height = 12, units = "cm",  bg = "transparent")