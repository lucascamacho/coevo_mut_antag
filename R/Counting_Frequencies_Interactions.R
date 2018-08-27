# Counting the frequencies of AA, AM and MM in the M and V matrices.
# load functions and packages
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
source("CoevoMutAntNet.R")
source("Counting.R")

library(ggplot2)
library(reshape2)

# Define a antprob sequence and a final results matrix
antprob = seq(0.01, 1, 0.01)
data = matrix(NA, ncol = 7, nrow = length(antprob))

for(i in 1:length(antprob)){ 
# loop to count the frequencies of aa, am and mm

  n_sp = 100 # number of species
  n_int = ((n_sp ** 2) - n_sp) / 2 # number of interactions in the matrix
  M = matrix(1, ncol = n_sp, nrow = n_sp) # Mutualism matrix
  diag(M) = 0 # no intraespecific interactions
  V = M * 0 # create antagonism interactions

  # tranforming links
  P = matrix(runif(n_sp*n_sp, min = 0, max = 1),
             nrow = n_sp, ncol = n_sp) # P is a sort matrix
  V[antprob[i] >= P] = 1 # if P is higher than antprob, transform to antagonism
  M[antprob[i] >= P] = 0 # and turn that mutualism off
  diag(V) <- 0 # no intraespecific interactions

  # counting interactions frequencies
  c = Counting(M, V) #apply Counting function

  # allocate frequencies in data matrix
  data[i,1] = antprob[i] # antprob values
  data[i,2] = antprob[i] * antprob[i] # expected p ^ 2
  data[i,3] = c[[1]] / n_int # observed AA frequencies
  data[i,4] = (2 * antprob[i]) * (1 - antprob[i]) # expected 2p(1 - p)
  data[i,5] = c[[2]] / n_int ## observed AM frequencies
  data[i,6] = (1 - antprob[i]) * (1 - antprob[i]) # expected (1 - p) ^ 2
  data[i,7] = c[[3]] / n_int # observed MM frequencies

}

# prepare to plot and plotting the results
data = data.frame(data)
colnames(data) = c("antprob", "esp_AA", "obs_AA", "esp_AM", "obs_AM", "esp_MM", "obs_MM")
test_data_long <- melt(data, id="antprob")  # convert to long format
plotar = ggplot(data=test_data_long,
         aes(x=antprob, y=value, colour=variable)) +
         geom_point() +
         theme_bw()

plotar