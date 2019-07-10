# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  <- gsub(".txt", replacement= "", temp)

# combinations of antprob and prob_change
antprob_vec = seq(0, 1, 0.1)
prob_change_vec = seq(0, 1, 0.1)
combs = expand.grid(antprob_vec, prob_change_vec)
combs = combs[rep(seq_len(nrow(combs)), 20), ]

# create sheet to alocate the simulations results
data = matrix(NA, nrow = 0, ncol = 8)
colnames(data) = c("net", "p", "q", "var", "rich", "conec", "autoval", "max_autoval")
data = as.data.frame(data)

for(i in 1:length(redes)){
  print(i)
  
  # 100 simulations for each network
  for(j in 1:nrow(combs)){
    M = as.matrix(redes[[i]])
    M[which(M > 1)] = 1
    M = SquareMatrix(M)
    n_sp = ncol(M)
    
    antprob = combs[j, 1]
    prob_change = combs[j, 2]
    
    connect = round((sum(M) / 2) / (ncol(M) ** 2), digits = 3)
    autoval = round(max(Re(eigen(M)$values)), digits = 3)
    max_autoval = round(sqrt(sum(M) / 2), digits = 3)
    
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Empirical_ConDep_Disparity.R")
    
    variance = round(variance[length(variance)], digits = 3)
    rich = as.numeric(ncol(M))
    n = names(redes)[[i]]
    
    results = data.frame(n, antprob, prob_change, variance, rich, connect, 
                         autoval, max_autoval)
    
    data = rbind(data, results) 
    
  }
}

save(data, file = "~/Dropbox/Master/Code/coevo_mut_antag/data/planilha_final.RData")

# write a table to save the results
write.table(data, file = "planilha_final.txt", sep = ",", 
            row.names = FALSE, col.names = TRUE)