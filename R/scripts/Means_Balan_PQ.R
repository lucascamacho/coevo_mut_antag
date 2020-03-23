# Code to run the coevolutionary process with different values of g and P.
# This script create all possible combinations between values of p and g
# - antprob (p) = (0.2, 0.5, 0.8)
# - prob_change (g) = (0, 0.01, 0.1).
#
# This combinations will be used to run different scenarios of coevolution
# in interactions shifting in time. Then, we will calculate the 
# species trait disparity
#
# This script returns one graph of directionality of species and another
# of divergency of species. Both graphs with 9 different scenarios of p and g.

# set work directory
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/data/")

# load packages
library(ggplot2)
library(reshape2)

# define the vector of p and q and combine 2-by-2 these values
antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0, 0.1, 0.5)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 1000), ]

# create a results matrix
data = matrix(NA, nrow = 9000, ncol = 4)
colnames(data) = c("mean", "var", "antprob", "prob_change")

# loop to define the combination of p and g and run the coevolution model
for(i in 1:nrow(data)){
  print(i)
  antprob = combs[i,1]
  prob_change = combs[i,2]
  
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/ConDep_BalanDiver.R")
  
  # keep the results in the data matrix
  data[i,1] = mean(balan)
  data[i,2] = diver[length(diver)]
  data[i,3] = antprob
  data[i,4] = prob_change
}

# save or load the data file
save(data, file = "Mean_Balan_PQ.RData")
load(file = "Mean_Balan_PQ.RData")

# plot the measures of Divergency and Directionality by different p and g combinations
data = as.data.frame(data)

box_plot_mean = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = mean, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = mean, 
                  fill = as.factor(prob_change)), alpha = 0.2) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  scale_fill_discrete(name = "Q") +
  labs(x = "Probability of interaction shift in time", 
       y = "Species trait directionality")

box_plot_var = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = var, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.3, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = var, 
                  color  = as.factor(prob_change))) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  scale_fill_discrete(name = "Q") +
  labs(x = "Probability of interaction shift in time", 
       y = "Species trait disparity")

# save the plots
ggsave(box_plot_mean, fil = "boxplot_directionality.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")
ggsave(box_plot_var, filename = "boxplot_disparity.pdf", 
       dpi = 600, width = 12, height = 8, units = "in")