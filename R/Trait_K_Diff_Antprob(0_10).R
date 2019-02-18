# Great loop to describe how the mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# in each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.
# For each value of antprob, we are doing 10 simulations

# set work directory and define antprob sequence
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
antprob_vec = seq(0.01, 1, 0.01)
library(dplyr)

# second data frame for last line of z_mat for each simulation
last_traits = matrix()
colnames(last_traits) = c("mean", "var", "type")

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){ #
  antprob = antprob_vec[a] # define a value fo antprob
  print(antprob) # print this value to follow the simulation process

  # OTIMIZAR -----------  
  col_variables = matrix()
  colnames(col_variables) = NULL
  
  # for each value of antprob, simulate 10 times
  for(y in 1:100){ 
    source("Trait_Diff_K_AM_MM.R")
    g.data = data[!duplicated(data$X3, fromLast=TRUE),]
    g.data[1,2] = g.data[1,2] / c[[2]]
    g.data[2,2] = g.data[2,2] / c[[3]]
    
    col_variables = bind_rows(as.data.frame(col_variables), g.data)
    
  }
  
  col_variables = col_variables[-1,-1]
  colnames(col_variables) = c("mean", "var", "type")
  
  final_m_AM = mean(col_variables$mean[which(col_variables$type == "AM")])
  final_var_AM = mean(col_variables$var[which(col_variables$type == "AM")])
  final_m_MM = mean(col_variables$mean[which(col_variables$type == "MM")])
  final_var_MM = mean(col_variables$var[which(col_variables$type == "MM")])
  
  quase = matrix(c(final_m_AM, final_m_MM, final_var_AM, final_var_MM), nrow = 2, ncol = 2)
  
  last_traits = bind_rows(as.data.frame(last_traits), as.data.frame(quase))
  
}

last_traits = last_traits[-1,]
last_traits[,3] = rep(c("AM", "MM"), times = 100)
last_traits[,4] = rep(1:100, each=2) / 100


## ggplot
plot = ggplot(data = last_traits) + geom_points(aes(x = last_traits$V4, y = last_traits$V1, 
                                                    color = last_traits$V3))