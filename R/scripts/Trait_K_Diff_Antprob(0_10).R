# Loop to describe how the mean trait values in equilibrium change with the antprob.
# For that we simulate the coevolutionary process 100 times for each value of antprob (p)
# ranging from 0.01 to 1 by 0.01.
#
# We are calculating the mean trait value for two different interaction outcomes types
# AM and MM. iNn each simulation, we have 3 values of species trait at "equilibrium", 
# each value representing an interaction type.
#
# This script returns 3 graphs of trait measures by antprob values (p):
#  - Number of timesteps to equilibrium
#  - Log of species trait variance in equilibrium
#  - Mean trait difference of species separated by interaction outcomes type
#
# The script probably will have some # symbols in parameters like antprob.
# If you will use only this code without a loop, be sure to get the #'s off.

# set work directory
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

# load packages
library(ggplot2)
library(reshape2)

# define vector of antprob values
antprob_vec = seq(0.01, 1, 0.01)

# data frame for var and time results
data_time_var = matrix(NA, nrow = 100, ncol = 5)
data_time_var[,1] = antprob_vec

# data frame for AM and MM trait difs
data_diffs = matrix(NA, nrow = 100, ncol = 3)
data_diffs[,1] = antprob_vec
colnames(data_diffs) = c("pvalue", "AM", "MM")

for(a in 1:length(antprob_vec)){
  # loop for every antprob value
  antprob = antprob_vec[a]
  print(antprob)

  # vectors to get the results
  data_times = vector()  
  data_vars = vector()

  data_am = vector()
  data_mm = vector()
  
  for(y in 1:100){
    # loop for repeat simulation at each antprob value
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

# save or load the data files
#save(data_time_var, file = "Data_Time_Var.RData")
#save(data_diffs, file = "Data_Diffs.RData")

#load("~/Dropbox/Master/Code/coevo_mut_antag/data/Data_Time_Var.RData")
#load("~/Dropbox/Master/Code/coevo_mut_antag/data/Data_Diffs.RData")

# plot the time, the variance and the difs varying antprob (p)
time_plot = ggplot(data = as.data.frame(data_time_var)) + 
            geom_point(aes(x = data_time_var[,1], y = data_time_var[,2]), 
                       alpha = 0.7, size = 2) +
            geom_errorbar(data = as.data.frame(data_time_var), 
                          mapping = aes(x = data_time_var[,1], ymax = data_time_var[,4], 
                                        ymin = data_time_var[,5])) +
            theme_minimal(base_size = 16) +
            xlab("Frequency of antagonisms (p)") +
            ylab("Number of timesteps to equilibrium")

var_plot = ggplot(data = as.data.frame(data_time_var)) + 
           geom_point(aes(x = data_time_var[,1], y = log(data_time_var[,3])), 
                      alpha = 0.7, size = 2) + 
           theme_minimal(base_size = 16) +
           xlab("Frequency of antagonisms (p)") + 
           ylab("Log of species trait variance in equilibrium")

d <- melt(data_diffs, id.vars = "pvalue")
difs_plot = ggplot(data = as.data.frame(data_diffs)) + 
             geom_point(aes(x = pvalue, y = AM, col="red"), size = 2) +
             geom_point(aes(x = pvalue, y = MM, col="blue"), size = 2) +
             theme_minimal(base_size = 16) +
             guides(color = guide_legend("Interação")) +
             scale_color_manual(labels = c("MM", "AM"), values = c("blue", "red")) +
             xlab("Frequency of antagonisms (p)") + 
             ylab("Mean trait difference of species separated by interaction outcomes type")

# save the plots
#ggsave(time_plot, file = "time.png", dpi = 600, 
#       width = 12, height = 8, units = "in")
#ggsave(var_plot, file = "variance.png", dpi = 600, 
#       width = 12, height = 8, units = "in")
#ggsave(difs_plot, file = "diffs.png", dpi = 600, 
#       width = 12, height = 8, units = "in")