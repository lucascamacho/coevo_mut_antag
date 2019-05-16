# Code to run the script of coevolution with context-dependent interactions and
# calculate 4 metrics of trait disparity for several values of P and Q.
# - antprob (p) = (0.2, 0.5, 0.8)
# - prob_change (q) = (0, 0.1, 0.5).
#
# This code runs the ConDep_Disparity.R 9000 times and for each time
# calculate 4 metrics of disparity.
# - Variance
# - Mean Pairwise Distance
# - Participation Ratio
# - Near Pairwise Distance with min and max
#
# This script returns a data.Rdata file with a data.frame of all the results
# and a plot for each disparity measure and different values of P and Q.

# set work directory
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts")

# load packages
library(ggplot2)
library(reshape2)

# define the vector of p and q and combine 2-by-2 these values
antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0, 0.1, 0.5)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 1000), ]

# create a results matrix
data = matrix(NA, nrow = 9000, ncol = 7)
colnames(data) = c("variance", "meanpairdist", "partratio", "neardist_min", 
                   "neardist_max", "antprob", "prob_change")

# loop to define the combination of P and Q and run the coevolution model
for(i in 1:nrow(data)){
  print(i)
  antprob = combs[i,1]
  prob_change = combs[i,2]
  
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/ConDep_Disparity.R")
  
  # keep the results in the data matrix
  data[i,1] = varian[length(varian)]
  data[i,2] = mediapar[length(mediapar)]
  data[i,3] = partiratio[length(partiratio)]
  data[i,4] = neardist[[nrow(traits)]][[1]]
  data[i,5] = neardist[[nrow(traits)]][[2]]
  data[i,6] = antprob
  data[i,7] = prob_change
}

# save or load the data file
#save(data, file = "~/Dropbox/Master/Code/coevo_mut_antag/data/ConDep_Disparity.RData")
load(file = "~/Dropbox/Master/Code/coevo_mut_antag/data/ConDep_Disparity.RData")
data = as.data.frame(data)

# prepare and plot the results
boxplot = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = neardist_min, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.3, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = neardist_min, 
                  color  = as.factor(prob_change))) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  scale_fill_discrete(name = "Q") +
  labs(x = "Valores de probabilidade de mudança de interação no tempo", 
       y = "Valores de Discrepância (Variância) das espécies")

boxplot