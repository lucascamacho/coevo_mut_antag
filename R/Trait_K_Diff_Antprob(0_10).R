# Great loop to describe how the last mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# for each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.

#sequence of antprob
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
library(ggplot2)
library(cowplot)

antprob_vec = seq(0.01, 1, 0.01)

# second data frame for last line of z_mat for each simulation
data2 = matrix(NA, ncol = 5, nrow = length(antprob_vec))
colnames(data2) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM", "antprob")
data2[,5] = antprob_vec

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a]
  print(antprob)
  
  col_variables = matrix(NA, nrow = 10, ncol = 4)
  colnames(col_variables) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM")
  
  for(y in 1:10){
    source("Trait_Diff_K_AM_MM.R")
    col_variables[y,1] = tail(data[6,1])  
    col_variables[y,2] = tail(data[6,2])
    col_variables[y,3] = tail(data[6,3])
    col_variables[y,4] = tail(data[6,4])

  }
  
  means = apply(col_variables, 2, mean)
  
  data2[a,1] = means[1]
  data2[a,2] = means[2]
  data2[a,3] = means[3]
  data2[a,4] = means[4]
  
  
}

# prepare data2 and plot
data2 = data.frame(data2)
par(mfrow=c(2,2))
plot(data2$antprob, data2$MEAN_AM,col="red")
plot(data2$antprob, data2$VAR_AM,col="red")
plot(data2$antprob, data2$MEAN_MM,col="red")
plot(data2$antprob, data2$VAR_MM,col="red")