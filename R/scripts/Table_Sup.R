# Code to generate a table with the basic informations
# of all the empirical networks that i used. In this table
# each row is a network nammed and nummered, with informations
# about mutualism type, richness, connectance, where this networks
# was sampled and the paper reference of the matrix

# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

library(bipartite)

# create a single matrix with all networks informations
networks = matrix(NA, nrow = 24, ncol = 4)
colnames(networks) = c("network", "L", "Q", "NODF")

# read all the .txt files
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes) = gsub(".txt", replacement = "", temp)

# loop to get the Name, Richness and Connectance of each network
for(i in 1:length(redes)){
  networks[i,1] = names(redes[i])
  networks[i,2] = (sum(redes[[i]])) / ((dim(redes[[i]])[1] + dim(redes[[i]])[2]) ** 2)
  networks[i,3] = 0
  networks[i,4] = nested(redes[[i]], method = "NODF2")
}

# save the table in a excell file
write.table(networks, file = "infos_2.txt", sep = ",",
            row.names = FALSE, col.names = TRUE)
