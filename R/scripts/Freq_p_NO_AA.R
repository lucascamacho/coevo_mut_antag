# Consistency test to check the frequency of antagonism interactions,
# after we exclude the interference interactions.
#
# We gonna count the frequencies of AM and MM in the M and V matrices.
# then check if p = freqAM / freqMM + freqAM as expected.
#
# This script returns a plot showing the xpected and observed values of AM
# for different values of p (called antprob).

# load functions and packages
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")

library(ggplot2)
library(reshape2)

# Define a antprob sequence and a final results matrix
antprob = seq(0.01, 1, 0.01)
data = matrix(NA, ncol = 3, nrow = length(antprob))
data[,3] = antprob

# loop to count the frequencies of AA, AM and MM
for(i in 1:length(antprob)){
  n_sp = 1000 # number of species
  n_int = ((n_sp ** 2) - n_sp) / 2 # number of interactions in the matrix
  M = matrix(1, ncol = n_sp, nrow = n_sp) # matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob[i])
  M = antagonize[[1]]
  V = antagonize[[2]]
  
  # counting interactions AA, AM and MM (AA must be zero)
  c = Counting(M, V)
  
  # if are some interference AA, stop the loop
  if(c[[1]] != 0){
    break
  }
  
  # p is the antprob value in a particular timestep
  # frequency of antagonisms AM without the interference AA outcomes
  
  # expected frequency of antagonisms
  data[i,1] = ((2 * antprob[i]) * (1 - antprob[i])) / 
              (((2 * antprob[i]) * (1 - antprob[i])) + ((1 - antprob[i]) ** 2))
  
  # observed frequency of antagonisms
  data[i,2] = (c[[2]] / n_int) / 
              ((c[[2]] / n_int) + (c[[3]] / n_int))
}

# prepare, plot and save the plot
data = data.frame(data)
colnames(data) = c("esp_p*", "obs_p*", "antprob")
test_data_long = melt(data, id = "antprob")  # convert to long format

plotar = ggplot(data = test_data_long,
                aes(x = antprob, y = value, colour = variable)) +
  geom_point(alpha = 0.7) +
  theme_bw()

plotar
#ggsave(plotar, file = "Freq_p_NO_AA.png")