# teste de consistencia
# numero de interations shifts muda com a diversidade de outcomes
# espero que redes menos diversas tenham menos trocas
# simular condep 100 redes para cada valor de P
# e guardar o length de w_time
# plotar uma curva de p pela frequencia de trocas

#setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/ConDepCoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MutAntag.R")

antprob_vec = seq(0.01, 1, 0.01)
dados = matrix(NA, nrow = 100, ncol = 2)
colnames(dados) = c("antprob", "m_trocas")

for(i in 1:length(antprob_vec)){
#  print(i)
  antprob = antprob_vec[i]
  prob_change = 0.9
  
  shifts = vector()
  
  for(j in 1:10){
  n_sp = 5 # defining number of species
  M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob)
  M = antagonize[[1]]
  V = antagonize[[2]]
  
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
  w_time = length(as.matrix(simulation[[2]]))
  shifts = append(shifts, w_time)
  }

  new_shifts = mean(shifts)
  dados[i,1] = antprob
  dados[i,2] = new_shifts
    
}

plot(dados[,1], dados[,2], xlab = "antprob (P)", 
     ylab = "Média de shifts de interação nas simulações")