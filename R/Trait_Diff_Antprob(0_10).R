# Great loop to describe how the last mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# for each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.

#sequence of antprob
antprob_vec = seq(0.01, 1, 0.01)

# second data frame for last line of z_mat for each simulation
data2 = matrix(NA, ncol = 4, nrow = length(antprob_vec))
colnames(data2) = c("AA", "AM", "MM", "antprob")
data2[,4] = antprob_vec

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a]

  source("Trait_Diff_AA_AM_MM.R")

  data2[a,1] = tail(data)[6,1]
  data2[a,2] = tail(data)[6,2]
  data2[a,3] = tail(data)[6,3]
      
}

# prepare data2 and plot
data2 = data.frame(data2)
test_data_long = melt(data2, id = "antprob")  # convert to long format

# plot data2
plotar = ggplot(data = test_data_long,
                aes(x = antprob, y = value, color = variable)) +
         geom_point(alpha = 0.7) +
         theme_bw()

plotar
