setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts")

mlg = matrix(NA, nrow = 10000, ncol = 3)
colnames(mlg) = c("antprob", "prob_change", "variance")

for( i in 1:nrow(mlg)){
  print(i)
  antprob = runif(1, 0, 1)
  prob_change = runif(1, 0, 1)

  source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/ConDep_Disparity.R")

  mlg[i,1] = antprob
  mlg[i,2] = prob_change
  mlg[i,3] = variance[length(variance)]

}


#save(mlg, file = "~/Dropbox/Master/Code/coevo_mut_antag/data/ConDep_MultRegr.RData")
load(file = "~/Dropbox/Master/Code/coevo_mut_antag/data/ConDep_MultRegr.RData")

mlg = as.data.frame(mlg)
modelo = lm(variance ~ antprob + prob_change + antprob:prob_change, data = mlg)
modelo3 = lm(variance ~ antprob, data = mlg)
anova(modelo, modelo2)

summary(modelo)
coefficients(modelo)
vcov(modelo)