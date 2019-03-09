# Consistency test
# Counting the frequencies of AM and MM in the M and V matrices.
# then check if p = freqAM / freqMM + freqAM

# load functions and packages
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("Antagonize.R")
source("EndInteraction.R")
source("Counting.R")

library(ggplot2)
library(reshape2)

# Define a antprob sequence and a final results matrix
antprob = seq(0.01, 1, 0.01)
data = matrix(NA, ncol = 3, nrow = length(antprob))
data[,3] = antprob

for(i in 1:length(antprob)){ 
  # loop to count the frequencies of AA, AM and MM
  
  n_sp = 100 # number of species
  n_int = ((n_sp ** 2) - n_sp) / 2 # number of interactions in the matrix
  M = matrix(1, ncol = n_sp, nrow = n_sp) # Mutualism matrix
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform links in antagonisms)
  antagonize = Antagonize(M, antprob[i])
  M = antagonize[[1]]
  V = antagonize[[2]]
  
  # End the AA interactions
  end = EndInteraction(M, V, "antagonism")
  M = end[[1]]
  V = end[[2]]

  # counting interactions AA, AM and MM (AA must be zero)
  c = Counting(M, V)
  
  # the p value is the antprob in a particular timestep
  # frequency of cheaters without the AA interactions
  
  #expected frequency of cheaters
  data[i,1] = ((2 * antprob[i]) * (1 - antprob[i])) / 
              (((2 * antprob[i]) * (1 - antprob[i])) + ((1 - antprob[i]) ** 2))
  
  #observed frequency of cheaters
  data[i,2] = (c[[2]] / n_int) / 
              ((c[[2]] / n_int) + (c[[3]] / n_int))
}

# prepare to plot and plotting the results
data = data.frame(data)
colnames(data) = c("esp_p*", "obs_p*", "antprob")
test_data_long = melt(data, id="antprob")  # convert to long format

plotar = ggplot(data = test_data_long,
                aes(x = antprob, y = value, colour = variable)) +
  geom_point(alpha = 0.7) +
  theme_bw()

ggsave(plotar, file = "Freq_p_NO_AA.pdf")