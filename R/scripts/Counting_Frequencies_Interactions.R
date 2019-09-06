# Consistency test for the frequency of antagonism outcomes (AM).
#
# Counting the frequencies of AA, AM and MM in the M and V matrices.
# We define p (or antprob) as the probability of an positive interaction outcome
# become negative.
#
# We expect that the frequencies of AA, AM and MM follow the proportions:
# AA = p ^ 2, AM = 2p(1 - p), MM = (1 - p) ^ 2.
#
# This script returns a plot with the expected and observed frequencies 
# of AA, AM and MM.

# load functions and packages
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")

library(ggplot2)
library(reshape2)

# Define a antprob sequence and a final results matrix
antprob = seq(0.01, 1, 0.01)
data = matrix(NA, ncol = 7, nrow = length(antprob))

for(i in 1:length(antprob)){ 
# loop to count the frequencies of AA, AM and MM

  n_sp = 100 # number of species
  n_int = ((n_sp ** 2) - n_sp) / 2 # number of interactions (double outcomes) in the matrix
  M = matrix(1, ncol = n_sp, nrow = n_sp) # matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob[i])
  M = antagonize[[1]]
  V = antagonize[[2]]

  # counting interactions outcomes frequencies
  c = Counting(M, V)

  # allocate frequencies in data matrix
  data[i,1] = antprob[i] # antprob values
#  data[i,2] = antprob[i] * antprob[i] # expected p ^ 2
#  data[i,3] = c[[1]] / n_int # observed AA frequencies
  data[i,4] = ((2 * antprob[i]) * (1 - antprob[i])) / (((2 * antprob[i]) * (1 - antprob[i])) + ((1 - antprob[i]) ** 2))
#  data[i,4] = (2 * antprob[i]) * (1 - antprob[i]) # expected 2p(1 - p)
  data[i,5] = c[[2]] / n_int ## observed AM frequencies
  data[i,6] = (1 - antprob[i]) * (1 - antprob[i]) # expected (1 - p) ^ 2
  data[i,7] = c[[3]] / n_int # observed MM frequencies

}

# prepare to plot and plotting the results
data = data.frame(data)
colnames(data) = c("antprob", "esp_AA", "obs_AA", "esp_AM", "obs_AM", "esp_MM", "obs_MM")
test_data_long = melt(data, id="antprob")  # convert to long format

plotar = ggplot(data = test_data_long,
         aes(x = antprob, y = value, colour = variable)) +
  geom_point(alpha = 0.7, size = 3) +
  ylab("Frequency of outcomes in adjacency matrix") +
  xlab("Probability of outcomes shift (p)")

# plot and save the plot
plotar
#ggsave(plotar, file = "Counting_Frequencies_Interactions.png")