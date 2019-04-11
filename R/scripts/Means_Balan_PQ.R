# Master code to run the ConDep_BalanDiver.R file 1000 times.
# We will get the mean Balan for each simulation and plot
# the results in boxplots of antprob = (0.2, 0.5, 0.8)
# and prob_change = (0, 0.01, 0.1). 

# set work directory and define antprob sequence
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts")
library(ggplot2)
library(reshape2)

antprob_vec = c(0.2, 0.5, 0.8)
prob_change_vec = c(0, 0.1, 0.5)
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

#save(data, file = "~/Dropbox/Master/Code/coevo_mut_antag/data/RUNFREE_Mean_Balan_PQ.RData")
#load(file = "~/Dropbox/Master/Code/coevo_mut_antag/data/Mean_Balan_PQ.RData")

data = as.data.frame(data)

box_plot_mean = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = mean, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.4, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = mean, 
                  fill = as.factor(prob_change))) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  scale_fill_discrete(name = "Q") +
  labs(x = "Valores de probabilidade de mudança de interação no tempo", 
       y = "Média dos valores de Direcionalidade das espécies")

box_plot_var = ggplot(data = data) +
  geom_point(aes(x = as.character(prob_change), y = var, 
                 color = as.factor(prob_change)), size = 1, 
             shape = 21, alpha = 0.3, position = position_jitterdodge()) +
  geom_violin(aes(x = as.character(prob_change), y = var, 
                  color  = as.factor(prob_change))) +
  facet_wrap(~antprob) +
  theme_bw(base_size = 16) +
  scale_fill_discrete(name = "Q") +
  labs(x = "Valores de probabilidade de mudança de interação no tempo", 
       y = "Valores de Discrepância (Variância) das espécies")

#ggsave(box_plot_mean, fil = "boxplot_balancancia.pdf", 
#       dpi = 600, width = 12, height = 8, units = "in")
#ggsave(box_plot_var, filename = "boxplot_discrepancia.pdf", 
#       dpi = 600, width = 12, height = 8, units = "in")