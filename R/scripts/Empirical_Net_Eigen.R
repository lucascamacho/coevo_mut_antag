# For empirical networks only
# Loop to choose each empirical network in a list of 24.
# 8 of Myrmecophye-Ant; 8 of Frugivory; 8 of Pollination
#
# For each of these networks, run 100 simulations
# and for each simulation, run a value of P and Q
# calculate the variance, richness, connectance, 
# real part of the higher eigenvalue and square root
# of interaction number.

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

# create sheet to alocate the simulations results
data = matrix(NA, nrow = 0, ncol = 8)
colnames(data) = c("net", "p", "q", "var", "rich", "conec", "autoval", "max_autoval")
data = as.data.frame(data)

# loop to simulate 100 times for each empirical network
for(i in 1:length(redes)){
 print(i)
  
  # 100 simulations for each network
  for(j in 1:100){
    M = as.matrix(redes[[i]])
    M[which(M > 1)] = 1
    M = SquareMatrix(M)
    n_sp = ncol(M)
    
    antprob = runif(1, 0, 1)
    prob_change = runif(1, 0, 1)
    
    # Antagonize M (transform positive links in negative)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
    
    # get basic matrix informations
    connect = round((sum(M) / 2) / (ncol(M) ** 2), digits = 3)
    autoval = round(max(Re(eigen(M)$values)), digits = 3)
    max_autoval = round(sqrt(sum(M) / 2), digits = 3)
    
    # coevolutionary model parameters
    phi = 0.2
    alpha = 0.2
    theta = runif(n_sp, 0, 10)
    init = runif(n_sp, 0, 10)
    p = 0.1
    epsilon = 5
    eq_dif = 0.0001
    t_max = 1000
    
    # running coevolution simulation
    simulation = ConDepCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, 
                                      eq_dif, t_max, prob_change)
    traits = as.matrix(simulation[[1]])
    w_time = as.matrix(simulation[[2]])
    
    # an apply for each line of z_mat
    variance = apply(traits, 1, var)
    
    # round variance number
    variance = round(variance[length(variance)], digits = 3)
    rich = as.numeric(ncol(M))
    n = names(redes)[[i]]
    
    # insert number in data frame and bind the results
    results = data.frame(n, antprob, prob_change, variance, rich, connect, 
                         autoval, max_autoval)

    data = rbind(data, results) 
    
  }
}

#save or load out data
save(data, file = "~/Dropbox/Master/Code/coevo_mut_antag/data/dados_planilha.RData")
load("dados_planilha.RData")

# write a table to save the results
write.csv(data, file = "planilha_simulacoes.csv", sep = ",", 
           row.names = FALSE, col.names = TRUE)