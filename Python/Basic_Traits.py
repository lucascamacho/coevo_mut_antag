## Lucas A. Camacho
## Basic_Traits.R
## Create matrices of interactions and simulate the coevolution process

# Load packages and define workdirectory
import os
import numpy as np
from Antagonize import Antagonize
from EndInteraction import EndInteraction
from ZeroLines import ZeroLines

os.chdir("/Users/glyptodon/Dropbox/Master/Code/coevo_mut_antag/Python/")

# initial parameters
antprob = 0.7
n_sp = 5

# create M and V matrices
M = np.ones((n_sp, n_sp))
np.fill_diagonal(M, 0)

# antagonize M and generate V matrix
antagonize = Antagonize(M, antprob)
M = antagonize[0]
V = antagonize[1]
np.fill_diagonal(V, 0)

# end the pure antagonism
end = EndInteraction(M, V, "antagonism")
M = end[0]
V = end[1]

# check for zero lines
zero = ZeroLines(M, V, n_sp, antprob)
M = zero[0]
V = zero[1]

# coevolutionary models parameters
phi = 0.2
alpha = 0.2
theta = np.random.uniform(0.0, 10.0, n_sp)
init = np.random.uniform(0.0, 10.0, n_sp)
p = 0.1
epsilon = 5
eq_dif = 0.0001
t_max = 1000

# running the coevolutionary process
traits = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# generate plot

