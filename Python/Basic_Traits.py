## Lucas A. Camacho
## Basic_Traits.R
## Create matrices of interactions and simulate the coevolution process

# Load packages
import os
import numpy as np

os.chdir("/Users/glyptodon/Dropbox/Master/Code/coevo_mut_antag/R/")

# Define functions
def Antagonize(M, antprob):
    V = M * 0
    P = np.array(np.random.rand(n_sp, n_sp))
    
    M[antprob >= P] = 0
    V[antprob >= P] = 1

    return M, V

def End_Interaction(M, V, interaction):
    if(interaction != "antagonism" 
       and interaction != "antagonism" and interaction != "antagonism"):
        print ("Choose a valid interactions type")
    
    if(interaction == "antagonism"):
        index = V == np.transpose(V)
        V[index] = 0
        return M, V

    if(interaction == "cheaters"):
        index = V != np.transpose(V)
        V[index] = 0
        M[index] = 0
        return M, V
    
    else:
        index = M == np.transpose(M) # end mutualism based on index of M and M transpose
        M[index] = 0
        return M, V
        

def Zero_Lines(M, V, n_sp, antprob):
		for z in range(0, 1000000)
				
		
		    
# initial parameters
antprob = 0.1
n_sp = 5

# create M and V matrices
M = np.ones((n_sp, n_sp))

# antagonize M and generate V matrix
antagonize = Antagonize(M, antprob)
M = antagonize[0]
V = antagonize[1]

# end the pure antagonism
end = End_Interaction(M, V, "antagonism")
M = end[0]
V = end[1]

# check for zero lines
zero = Zero_Lines(M, V, n_sp, antprob)
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
coevo = Coevo_Mut_Ant_Net(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

# generate plot



