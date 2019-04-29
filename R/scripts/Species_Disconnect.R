# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EndInteraction.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Disconnected.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")

library(ggplot2)
library(reshape2)
library(cowplot)

antprob_vec = seq(0.01, 1, 0.01)
results = matrix(NA, nrow = length(antprob_vec), ncol = 2)
results[,1] = antprob_vec

for(a in 1:length(antprob_vec)){
# initial parameters
  antprob = antprob_vec[a]
#  print(antprob)
  n_sp = 10 # defining number of species
  M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions

  V = M * 0
  
#  for(i in 1:dim(M)[1]){
#    for(j in 1:dim(M)[2]){
#      if(M[i,j] == 1){
#        p = runif(1, 0, 1)
#        
#        if(p <= antprob){
#          M[j,i] = 0
#          V[j,i] = 1
#        }  
#      }
#    }
#  }
  
 #Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob)
  M = antagonize[[1]]
  V = antagonize[[2]]

  # End interferences AA
  end = EndInteraction(M, V, "interference")
  M = end[[1]]
  V = end[[2]]

  # parameters to cut interactions and cut
  init = runif(n_sp, 0, 10)
  epsilon = 2
  A = M + V 
  z_dif = t(A * init) - A * init
  d = which(abs(z_dif) > epsilon)
  A[d] = 0

  # Ho and how many species are disconnected
  disco = Disconnected(M, V)
  index = disco[[1]]
  num = disco[[2]]

  print(Counting(M, V)[[1]])
  results[a,2] = num
}

plot(results[,1], results[,2], main = "rede com 10 espécies", xlab = "p", ylab = "numero de espécies desconectadas")