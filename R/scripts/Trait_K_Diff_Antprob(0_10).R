# Loop to describe how the mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# in each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.
# For each value of antprob, we are doing 10 simulations

# set work directory and define antprob sequence
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
library(ggplot2)
library(ggridges)
library(reshape2)

antprob_vec = seq(0.01, 1, 0.01)

# data frame for var and time
data_time_var = matrix(NA, nrow = 100, ncol = 5)
data_time_var[,1] = antprob_vec

# data frame for am and mm diff traits
data_diffs = matrix(NA, nrow = 100, ncol = 3)
data_diffs[,1] = antprob_vec
colnames(data_diffs) = c("pvalue", "AM", "MM")

for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a]
  print(antprob)

  data_times = vector()  
  data_vars = vector()

  data_am = vector()
  data_mm = vector()
  
  for(y in 1:100){
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Trait_Diff_K_AM_MM.R")
    
    data_times = append(data_times, tail(var_time)[6,1])
    data_vars = append(data_vars, tail(var_time)[6,2])
    
    data_am = append(data_am, tail(dif_time_AM)[6,2])
    data_mm = append(data_mm, tail(dif_time_MM)[6,2])
  }
  
  data_time_var[a,2] = mean(data_times)
  data_time_var[a,3] = mean(data_vars)
  data_time_var[a,4] = tail(sort(data_times))[2]
  data_time_var[a,5] = head(sort(data_times))[5]
  
  data_diffs[a,2] = mean(data_am)
  data_diffs[a,3] = mean(data_mm)
}

# save the data files
save(data_time_var, file = "Data_Time_Var.RData")
save(data_diffs, file = "Data_Diffs.RData")

# Adicionar colunas 4 e 5 para maior e menor valor de tempo em cada 100 valores
time_plot = ggplot(data = as.data.frame(data_time_var)) + 
            geom_point(aes(x = data_time_var[,1], y = data_time_var[,2]), alpha = 0.7, size = 2) +
            geom_errorbar(data = as.data.frame(data_time_var), mapping = aes(x = data_time_var[,1], 
                                                              ymax = data_time_var[,4],
                                                              ymin = data_time_var[,5])) +
            theme_minimal(base_size = 16) +
            ggtitle("Change of time to reach equilibrium by P") +
            xlab("Frequency of antagonisms (P)") + 
            ylab("Average time to reach equilibrium of simulation")
ggsave(time_plot, file = "time.pdf", dpi = 600, width = 12, height = 8, units = "in")

var_plot = ggplot(data = as.data.frame(data_time_var)) + 
           geom_point(aes(x = data_time_var[,1], y = log(data_time_var[,3])), alpha = 0.7, size = 2) + 
           theme_minimal(base_size = 16) +
           ggtitle("Variance of species traits by P") +
           xlab("Frequency of antagonisms (P)") + 
           ylab("Log of Average variance of species traits")
ggsave(var_plot, file = "variance.pdf", dpi = 600, width = 12, height = 8, units = "in")

d <- melt(data_diffs, id.vars="pvalue")

diffs_plot = ggplot(data = as.data.frame(data_diffs)) + 
             geom_point(aes(x = pvalue, y = AM, col="red"), size = 2) +
             geom_point(aes(x = pvalue, y = MM, col="blue"), size = 2) +
             theme_minimal(base_size = 16) +
             guides(color=guide_legend("Interaction")) +
             ggtitle("Average total difference of species separated by type of interactions by P") +
             scale_color_manual(labels = c("MM", "AM"), values = c("blue", "red")) +
             xlab("Frequency of antagonisms (P)") + 
             ylab("Average total difference of species traits balanced by the degrees")

ggsave(diffs_plot, file = "diffs.png", dpi = 600, width = 12, height = 8, units = "in")