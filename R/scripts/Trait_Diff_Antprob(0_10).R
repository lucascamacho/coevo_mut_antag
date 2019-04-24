# Loop to describe how the mean trait values of species in equilibrium changes
# when we vary the p (antprob).
#
# We are calculating the mean trait value for different outcomes of interactions.
# We run, for each antprob value in a range from 0.01 to 1 the coevolutionary process
# then, we get the average of species traits for each group of interaction:
# AA, AM and MM. For each simulation, we have 3 values of species trait at "equilibrium", 
# each value representing a interaction type.

# This scripts returns a plot which has the mean trait values for 3 outcomes groups by antprob

# Load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/FindInteractors.R")

library(ggplot2)
library(reshape2)

#sequence of antprob
antprob_vec = seq(0.01, 1, 0.01)

# second data frame for last line of z_mat for each simulation
last_traits = matrix(NA, ncol = 4, nrow = length(antprob_vec))
colnames(last_traits) = c("AA", "AM", "MM", "antprob")
last_traits[,4] = antprob_vec

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a]
  print(antprob)
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Trait_Diff_AA_AM_MM.R")

  last_traits[a,1] = tail(data)[6,1]
  last_traits[a,2] = tail(data)[6,2]
  last_traits[a,3] = tail(data)[6,3]
      
}

# prepare last_traits and plot
last_traits = data.frame(last_traits)
test_data_long = melt(last_traits, id = "antprob")  # convert to long format

# plot last_traits
plotar = ggplot(data = test_data_long,
                aes(x = antprob, y = value, color = variable)) +
         geom_point(alpha = 0.7) +
         theme_bw()

plotar