# Simulate a coevolution process without the AA interactions
# then, compute the mean trait value for groups of interaction types
# balancing this mean by the quantity of interactions

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
source("CoevoMutAntNet.R")
source("Antagonize.R")
source("EndInteraction.R")
source("Counting.R")
source("FindInteractors.R")
source("SpDegree.R")

library(ggplot2)
library(cowplot)

antprob = 0.25  # current probability value
n_sp = 10   # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End antagonism AA
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

degree = SpDegree(M, V)

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# simulate coevolution and calculate mean trait values for each interaction type
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
index = FindInteractors(M, V)

data = matrix(NA, nrow = nrow(z_mat), ncol = 5)
data[,5] = seq(1,nrow(z_mat), 1)
colnames(data) = c("AV_AM", "VAR_AM", "AV_MM", "VAR_MM", "time")

for(i in 1:nrow(z_mat)){
  av = mean(z_mat[i,])
  var = var(z_mat[i,])
  
  # nova funcao contando AM e MM de cada espécie
  # usar essa contagem 
  
}
# building data frame to use in ggplot
traits = as.data.frame(traits)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                       time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)
# plotting traits through time
plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
  geom_line(size = 1.8, alpha = 0.6) + 
  ggtitle(paste("proportion antagonists = ", antprob)) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))

print(plotar)