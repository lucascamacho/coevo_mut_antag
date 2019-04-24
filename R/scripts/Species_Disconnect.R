# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/EndInteraction.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Disconnected.R")

library(ggplot2)
library(reshape2)
library(cowplot)

# initial parameters
antprob = 0.9 # current probability value
n_sp = 10 # defining number of species
M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
diag(M) = 0 # no intraespecific interactions

# Antagonize M (transform positive links in negative)
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

# End interferences AA
end = EndInteraction(M, V, "interference")
M = end[[1]]
V = end[[2]]

# Ho and how many species are disconnected
disco = Disconnected(M, V)
index = disco[[1]]
num = disco[[2]]

# 1. cortar especies pelo trait diff
# 2. aplicar disconnected para diferente valores de p
# 3. 
