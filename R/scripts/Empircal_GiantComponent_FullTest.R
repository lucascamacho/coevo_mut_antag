# Run the coevolutionary model with interaction outcomes that are context-dependent.
# Based on a probability prob_change (Q), in each timestep of the simulation, an interaction
# outcome shifts (AA -> AM and AM -> MM).
#
# This script returns a simple graph with species traits changing in time due to coevolution.
# The asteriscs in the graph shows the timesteps in which the interactions shift occurs.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(igraph)
library(dplyr)

# read all the empirical networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

# create data drame to allocate data
dados = data.frame()

# first, showing if the empirical networks has, indeed, a giant component.
# It will have if the average species degree > 1. 
for(a in 1:length(redes)){
  M = as.matrix(redes[[a]]) # M is the adjancency matrix of interactions
  M[which(M > 1)] = 1 # if there are any error, correct that  
  
  deg_av = mean(rowSums(M))
  
  if(deg_av > 1){
    gc = TRUE
  }
  else{
    gc = FALSE
  }
  
  net = names(redes[a])
  
  results = data.frame(net, gc)
  dados = rbind(dados, results)
}

# check if there are any FALSE in dados, if not, you can proceed
any(dados == FALSE)

# next, we gonna relate the antprob value and the average species degree
# in the V matrix of negative effects. With this, we show that increasing
# antprob, we increase average species degree in V matrix.

antprob_vec = seq(0.01, 1, 0.01)
dados = data.frame()

# first graph showing the relation between p and mean(k)
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a] # current probability value
  n_sp = 100 # defining number of species
  M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob)
  M = antagonize[[1]]
  V = antagonize[[2]]
  
  deg_av = mean(rowSums(V))
  
  results = data.frame(antprob, deg_av)
  dados = rbind(dados, results)
}

# plot the graph relating p and average species degree
plot_relate = ggplot(data = dados,
                     aes(x = antprob, y = deg_av)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_line(aes(x = antprob, y = deg_av), stat="smooth", method = "lm",
            alpha = 0.8, col = "red") +
  ylab("Average species degree of negative effects<Kv>") +
  xlab("Value of interaction shift probability (p)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 15), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

plot_relate

#ggsave(plot_relate, filename = "plot_relate_antprob_kv.png", dpi = 600,
#       width = 21, height = 14, units = "cm")

dados = data.frame()
# Finally, we can plot a graph showing the value of p in X axis and the 
# relative size of the biggest component in the Y axis for all the empirical matrices.
for(k in 1:length(redes)){
  print(k)
  for(a in 1:length(antprob_vec)){
    antprob = antprob_vec[a] # define antprob
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number

    # Antagonize M (transform positive links in negative)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
  
    A = graph_from_adjacency_matrix(V)
    comp_emp = max(components(A)$csize) / n_sp
    
    A = graph_from_adjacency_matrix(M)
    max_ng = max(components(A)$csize) / n_sp
    
    M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
    diag(M) = 0 # no intraespecific interactions
  
    # Antagonize M (transform positive links in negative)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
  
    A = graph_from_adjacency_matrix(V)
    comp_teo = max(components(A)$csize) / n_sp
    
    net = names(redes[k])
    rich = n_sp
    
    results = data.frame(net, rich, antprob, comp_emp, comp_teo, max_ng)
    dados = rbind(dados, results)
  }
}

type = c(rep("Pollination", 800), rep("Seed dispersal", 800), rep("Ant-Plant", 800))
dados = cbind(dados, type)

n_dados = dados%>%
  group_by(antprob, type)%>%
  summarise(mean_emp = mean(comp_emp), max_ng = max(max_ng))%>%
  as.data.frame()

plot_component = 
  ggplot(data = n_dados) +
  geom_point(aes(x = antprob, y = mean_emp, colour = type), show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean_emp, colour = type), stat = "smooth", method = "auto",
            alpha = 0.8, show.legend = FALSE) +
   geom_hline(data = n_dados, aes(yintercept = n_dados[1,4])) +
  geom_hline(data = n_dados, aes(yintercept = n_dados[3,4])) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.1)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Relative size of the biggest cheater component component") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

plot_component

ggsave(plot_component, filename = "component_cheaters.png", dpi = 600,
       width = 21, height = 14, units = "cm")