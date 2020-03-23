# Code to run the script of coevolution with interactions shifting in time and
# calculate 4 metrics of trait disparity for several values of p and g.
# - antprob (p) = (0.2, 0.5, 0.8)
# - prob_change (g) = (0, 0.1, 0.5).
#
# This code runs the ConDep_Disparity.R 9000 times and for each time
# calculate 4 metrics of disparity.
# - Variance
# - Mean Pairwise Distance
# - Participation Ratio
# - Near Pairwise Distance with min and max
#
# Make sure that the ConDep_Disparity funciton has #'s on certain lines
#
# This script returns a data.Rdata file with a data.frame of all the results
# and a plot for each disparity measure and different values of p and q.

# set work directory
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/data/")

# load packages
library(ggplot2)

# define the vector of p and q and combine 2-by-2 these values
antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0, 0.1, 0.5)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 1000), ]

# create a results matrix
data = matrix(NA, nrow = 9000, ncol = 7)
colnames(data) = c("variance", "meanpairdist", "partratio", "neardist_min", 
                   "neardist_max", "antprob", "prob_change")

# loop to define the combination of p and q and run the coevolution model
for(i in 1:nrow(data)){
  print(i)
  antprob = combs[i,1]
  prob_change = combs[i,2]
  
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/ConDep_Disparity.R")
  
  # keep the results in the data matrix
  data[i,1] = variance[length(variance)]
  data[i,2] = meanpairdist[length(meanpairdist)]
  data[i,3] = partratio[length(partratio)]
  data[i,4] = neardist[[nrow(traits)]][[1]]
  data[i,5] = neardist[[nrow(traits)]][[2]]
  data[i,6] = antprob
  data[i,7] = prob_change
}

# save or load the data file
save(data, file = "ConDep_Disparity.RData")
load(file = "ConDep_Disparity.RData")

# data in data.frame to plot
data = as.data.frame(data)

# prepare and plot the results
# var plot
var_boxplot = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = variance, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = variance, 
                  fill = as.factor(prob_change)), alpha = 0.2) +
  facet_wrap(antprob) +
  theme_bw(base_size = 16) +
  labs(x = "Probability of interaction shift in time (g)", 
       y = "Variance")
  
# the mean pairwise distance plot
meanpairdist_boxplot = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = meanpairdist, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = meanpairdist, 
                  fill = as.factor(prob_change)), alpha = 0.2) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  labs(x = "Probability of interaction shift in time (g)", 
       y = "Mean Pairwise Distance")

# the participation ratio plot
partratio_boxplot = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = partratio, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = partratio, 
                  fill = as.factor(prob_change)), alpha = 0.2) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  labs(x = "Probability of interaction shift in time (g)", 
       y = "Participation Ratio")

# Near and Long Pairwise Distance plot
nearlong_boxplot = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = neardist_min, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_point(aes(x = as.character(prob_change), y = neardist_max, 
                  fill = as.factor(prob_change)), alpha = 0.2) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  labs(x = "Probability of interaction shift in time (g)", 
       y = "Near and Longest Pairwise Distance")

# save plots
ggsave(var_boxplot, fil = "var_boxplot.png")
ggsave(var_boxplot, fil = "var_boxplot.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")
ggsave(meanpairdist_boxplot, filename = "meanpairdist_boxplot.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")
ggsave(partratio_boxplot, fil = "partratio_boxplot.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")
ggsave(nearlong_boxplot, filename = "nearlong_boxplot.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")