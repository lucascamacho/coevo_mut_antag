# Great loop to describe how the last mean trait values in equilibrium change with the antprob
# remembering that we are calculating the mean trait value for different interaction types
# for each simulation, we have 3 values of species trait at "equilibrium", each value
# representing a interaction type.

#sequence of antprob
antprob_vec = seq(0.01, 1, 0.01)

# second data frame for last line of z_mat for each simulation
last_traits = matrix(NA, ncol = 4, nrow = length(antprob_vec))
colnames(last_traits) = c("AA", "AM", "MM", "antprob")
last_traits[,4] = antprob_vec

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){
  antprob = antprob_vec[a]

  source("Trait_Diff_AA_AM_MM.R")

  last_traits[a,1] = tail(data)[6,1]
  last_traits[a,2] = tail(data)[6,2]
  last_traits[a,3] = tail(data)[6,3]
      
}

# prepare last_traits and plot
last_traits = data.frame(last_traits)
test_data_long = melt(last_traits, id = "antprob")  # convert to long format

# plot last_traits
plotar = ggplot(data = test_data_long,
                aes(x = antprob, y = value, color = variable)) +
         geom_point(alpha = 0.7) +
         theme_bw()

plotar