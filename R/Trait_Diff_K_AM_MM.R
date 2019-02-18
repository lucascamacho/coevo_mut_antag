# Simulate a coevolution process without the AA interactions
# then, compute the mean trait value for groups of interaction types
# balancing this mean and variance of traits by the quantity of interactions.

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("Antagonize.R")
source("EndInteraction.R")
source("SpDegree.R")
source("Counting.R")
source("FindInteractors.R")
source("CoevoMutAntNet.R")
source("TraitDegreeBalanced.R")
source("VarTraitDegreeBalanced.R")

library(ggplot2)
library(reshape)
library(cowplot)
library(dplyr)

# initial conditions
#antprob = 0.5  # current probability value
n_sp = 100 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp)   # building matrix M (mutualisms)
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform links in antagonisms)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End the AA interactions
end = EndInteraction(M, V, "antagonism")
M = end[[1]]
V = end[[2]]

# measure the degree of AM and MM for each specie
degree = SpDegree(M, V)

# count the number of interactions
c = Counting(M, V)

# separate the species in groups by type of interaction
index = FindInteractors(M, V)

# coevolutionary model parameters
phi = 0.2
alpha = 0.2
theta = runif(n_sp, 0, 10)
init = runif(n_sp, 0, 10)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# simulate coevolution and calculate mean and variance trait values for each interaction type
z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# APPLY FUNCTIONS NEW
traits = t(apply(z_mat, 1, TraitDegreeBalanced))
colnames(traits) = c("AM", "MM")
var_traits = t(apply(z_mat, 1, VarTraitDegreeBalanced))
colnames(var_traits) = c("AM", "MM")

#OTIMIZAR ------
traitsAM = traits[,1]
type = rep("AM", times = length(traitsAM))
traitsAM = data.frame(traitsAM, type)

traitsMM = traits[,2]
type = rep("MM", times = length(traitsMM))
traitsMM = data.frame(traitsMM, type)

vartraitsAM = var_traits[,1]
type = rep("AM", times = length(vartraitsAM))
vartraitsAM = data.frame(vartraitsAM, type)

vartraitsMM = var_traits[,2]
type = rep("MM", times = length(vartraitsMM))
vartraitsMM = data.frame(vartraitsMM, type)

traitsAM = as.matrix(traitsAM)
traitsMM = as.matrix(traitsMM)

colnames(traitsAM) = NULL
colnames(traitsMM) = NULL

new_traits = rbind(traitsAM, traitsMM)

vartraitsAM = as.matrix(vartraitsAM)
vartraitsMM = as.matrix(vartraitsMM)

colnames(vartraitsAM) = NULL
colnames(vartraitsMM) = NULL

new_vartraits = rbind(vartraitsAM, vartraitsMM)

data = cbind(new_traits, new_vartraits)
data = data[,-2]

data = data.frame(data)

data[,1] = as.numeric(as.character(data[,1]))
data[,2] = as.numeric(as.character(data[,2]))

#