# Basic script to test the R base-functions to simulate the coevolutionary process
# of networks with positive and negative interactions outcomes.
# This script returns a simple graph with species traits changing in time.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
antprob = 0.1 # current probability value
n_sp = 10 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
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

# running coevolution simulation
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
z_mat = t(z_mat[1:10,])

z_mat_stand = scale(z_mat[-1])
k.means.fit = kmeans(z_mat_stand, 3)

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(z_mat_stand, nc=6) 
clusplot(z_mat_stand, k.means.fit$cluster, main='2D representation of the Cluster solution',
         color=TRUE, shade=TRUE,
         labels=2, lines=0)
