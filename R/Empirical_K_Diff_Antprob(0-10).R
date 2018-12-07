# For empirical networks
# Great loop to describe how the mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# in each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.
# For each value of antprob, we are doing 10 simulations

# set work directory and define antprob sequence
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
antprob_vec = seq(0.01, 0.99, 0.01)

# second data frame for last line of z_mat for each simulation
last_traits = matrix(NA, ncol = 5, nrow = length(antprob_vec))
colnames(last_traits) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM", "antprob")
last_traits[,5] = antprob_vec

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a] # define a value fo antprob
  print(antprob) # print this value to follow the simulation process
  
  # create a small data matrix to get the mean of several simulations
  col_variables = matrix(NA, nrow = 10, ncol = 4)
  colnames(col_variables) = c("MEAN_AM", "VAR_AM", "MEAN_MM", "VAR_MM")
  
  # for each value of antprob, simulate 10 times
  for(y in 1:10){
    setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
    source("Empirical_Trait_Diff_K_AM_MM.R")
    col_variables[y,1] = tail(data[6,1])  
    col_variables[y,2] = tail(data[6,2]) / c[[2]]
    col_variables[y,3] = tail(data[6,3])
    col_variables[y,4] = tail(data[6,4]) / c[[3]]
    
  }
  
  means = apply(col_variables, 2, mean) # calculate the mean of these 10 simulations
  
  last_traits[a,1] = means[1] # allocate these means in our final data matrix
  last_traits[a,3] = means[3]
  last_traits[a,2] = means[2]
  last_traits[a,4] = means[4]
  
}

# prepare final data and plot in a single window
data = data.frame(last_traits)
par(mfrow = c(2,2))
plot(data$antprob, data$MEAN_AM, pch = 19, col = "blue", xlab = "antprob (p)", ylab = "Mean Trait for Cheaters")
plot(data$antprob, data$MEAN_MM, pch = 19, col = "blue", xlab = "antprob (p)", ylab = "Mean Trait for Mutualism")
plot(data$antprob, data$VAR_AM, pch = 19, col = "red", xlab = "antprob (p)", ylab = "Delta Trait for Cheaters")
plot(data$antprob, data$VAR_MM, pch = 19, col = "red", xlab = "antprob (p)", ylab = "Delta Trait for Mutualism")
title("Traits of Cheaters and Mutualism (Balanced by degree Kmm and Kam)", line = -2, outer = TRUE)