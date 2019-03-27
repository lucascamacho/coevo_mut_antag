# Master code to run the ConDep_BalanDiver.R file 1000 times.
# We will get the mean Balan for each simulation and plot
# the results in boxplots of antprob = (0.2, 0.5, 0.8)
# and prob_change = (0, 0.01, 0.1). 

# set work directory and define antprob sequence
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts")
library(ggplot2)
library(reshape2)

antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0, 0.01, 0.1)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 1000), ]

data = matrix(NA, nrow = 9000, ncol = 4)
colnames(data) = c("mean", "var", "antprob", "prob_change")

for(i in 1:nrow(data)){
  print(i)
  antprob = combs[i,1]
  prob_change = combs[i,2]
  
  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/ConDep_BalanDiver.R")
  
  data[i,1] = mean(balan)
  data[i,2] = diver[length(diver)]
  data[i,3] = antprob
  data[i,4] = prob_change
}

load(file = "~/Dropbox/Master/Code/coevo_mut_antag/data/Mean_Balan_PQ.RData")
data = as.data.frame(data)

box_plot = ggplot(data = data) +
  geom_boxplot(aes(x = as.character(prob_change), y = mean)) +
  facet_wrap(~antprob) +
  theme_classic()

box_plot