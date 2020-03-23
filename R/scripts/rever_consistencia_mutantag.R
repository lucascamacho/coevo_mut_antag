# teste de consistencia MutAntag
# mantem numero de interacoes? OK
# mantem numero de AM e MM? Quando nao há AA ou AM disponível, so ocorre uma troca
# criou algum AA? OK
# matrizes A sao iguais? OK

setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MutAntag.R")

data = matrix(NA, nrow = 10000, ncol = 4)
colnames(data) = c("n_int", "antag", "mut", "AA")

data_mutantag = matrix(NA, nrow = 10000, ncol = 5)
colnames(data_mutantag) = c("n_int", "antag", "mut", "AA", "igual")

for(i in 1:nrow(data)){
  print(i)
  antprob = 0.5 # current probability value
  prob_change = 0.9 # current probability of interaction outcome shift
  n_sp = 5 # defining number of species
  M = matrix(1, ncol = n_sp, nrow = n_sp) # building matrix M of positive outcomes
  diag(M) = 0 # no intraespecific interactions
  
  # Antagonize M (transform positive links in negative)
  antagonize = Antagonize(M, antprob)
  M = antagonize[[1]]
  V = antagonize[[2]]
  
  A = M + V
  n_int = length(which(A == 1))
  
  n_antag = (length(M[M != t(M)]) / 2)
  
  m1 = which(M == 1)
  m2 = which(t(M) == 1)
  n_mut = sum(table(m1[m1 %in% m2])) / 2
  
  v1 = which(V == 1)
  v2 = which(t(V) == 1)
  n_inter = sum(table(v1[v1 %in% v2])) / 2 
  
  data[i,1] = n_int
  data[i,2] = n_antag
  data[i,3] = n_mut
  data[i,4] = n_inter
  
  
  mutantag = MutAntag(M, V, 1, prob_change)
  M = mutantag[[1]]
  V = mutantag[[2]]
  
  J = M + V
  n_int_mutantag = length(which(J == 1)) # TEM QUE SER TUDO 45
  
  n_antag_mutantag = (length(M[M != t(M)]) / 2) # TEM QUE SER IGUAL ANTES
  
  m1 = which(M == 1)
  m2 = which(t(M) == 1)
  n_mut_mutantag = sum(table(m1[m1 %in% m2])) / 2 # TEM QUE SER IGUAL ANTES
  
  v1 = which(V == 1)
  v2 = which(t(V) == 1)
  n_inter_mutantag = sum(table(v1[v1 %in% v2])) / 2 # TEM QUE SER 0
  
  igual = identical(J, A)
  
  data_mutantag[i,1] = n_int_mutantag
  data_mutantag[i,2] = n_antag_mutantag
  data_mutantag[i,3] = n_mut_mutantag
  data_mutantag[i,4] = n_inter_mutantag
  data_mutantag[i,5] = igual
}

identical(data[,1], data_mutantag[,1]) #NUMERO DE INTERAÇÕES SE MANTÉM?
identical(data[,2], data_mutantag[,2]) #NUMERO DE ANTAGONISMO SE MANTÉM?
identical(data[,3], data_mutantag[,3]) #NUMERO DE MUTUALISMO SE MANTÉM?
any(data_mutantag[,4] > 0) #CRIA INTERFERÊNCIAS?
any(data_mutantag[,5] == 0) #REDES TEM ALGUMA DIFERENÇA?
