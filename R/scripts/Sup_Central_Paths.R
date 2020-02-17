# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

library(expm)

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CentralAntagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")


# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

M = as.matrix(redes[[1]]) # M is the adjancency matrix of interactions
M[which(M > 1)] = 1 # if there are any error, correct that
M = SquareMatrix(M) # square the adjancency matrix
n_sp = ncol(M) # define the species number

# insert cheaters outcomes based on degree centrality
centralantagonize = CentralAntagonize(M)
M = centralantagonize[[1]]
V = centralantagonize[[2]]

D = M
P = M

for(k in 2:n_sp){
  P = M %*% P
  for(i in 1:n_sp){
    for(j in 1:n_sp){
      if(P[i,j] > 0){
        if(D[i,j] == 0){
          D[i,j] = k
        }
      }
    }
  }
}

diag(D) = NA

#png("Central_hist_paths.png")
hist(D, main = "Central. Interações negativas estavam focadas em sp centrais")
#dev.off()

###
antprob = centralantagonize[[3]]

M = as.matrix(redes[[1]]) # M is the adjancency matrix of interactions
M[which(M > 1)] = 1 # if there are any error, correct that
M = SquareMatrix(M) # square the adjancency matrix
n_sp = ncol(M) # define the species number

# insert cheaters outcomes in the network
antagonize = Antagonize(M, antprob)
M = antagonize[[1]]
V = antagonize[[2]]

D = M
P = M

for(k in 2:n_sp){
  P = M %*% P
  for(i in 1:n_sp){
    for(j in 1:n_sp){
      if(P[i,j] > 0){
        if(D[i,j] == 0){
          D[i,j] = k
        }
      }
    }
  }
}

diag(D) = NA

#png("Dist_hist_paths.png")
hist(D, main = "Distributed. Interações negativas estavam distribuidas na rede")
#dev.off()