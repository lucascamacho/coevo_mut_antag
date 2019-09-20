# Last figure of the paper
# Code to explore the influence of cheaters exploitation in the structure
# of the adjancency matrix of interactions. We are considering a second
# trait barrier that will cut interactions off if this barrier is tranpassed
# This will allow us to see what's the effect of cheaters in coevolution and
# how this change's the strcuture of ecological interacions in the community.

# The nestedness will be calculated using the bipartite package and the modularity
# with the MODULAR program

# Marquitti, F. M. D., P. R. Guimaraes, M. M. Pires, L. F. Bittencourt. 2014. 
# MODULAR: Software for the Autonomous Computation of Modularity in Large Network Sets. 
# Ecography: 37: 221–224.

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(bipartite)

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

# create data.frame to store all my results
final_fl = data.frame()

for(k in 1:length(redes)){
  print(k) # show which network the simulation are using
  
  for(a in 1:30){ # loop to each network
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number
    n_int = sum(M) / 2 # number of interactions in the matrix
    
    antprob = runif(1, 0, 1) # sample an antprob value
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/DiscoCoevoMutAntNet.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
    
    # antagonize matrices (now I have an M and V adjancency matrix)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
    
    c = Counting(M,V)
    antprob = c[[2]] / n_int
    
    # DesQuad matrices
    li = 1:dim(redes[[k]])[1]
    co_1 = dim(redes[[k]])[1]+1
    co_2 = dim(redes[[k]])[1]+dim(redes[[k]])[2]
    
    # init_m and W are non squared adjacency matrices
    init_m = as.matrix(M + V)[li, co_1:co_2]
    W = as.matrix(M + V)[li, co_1:co_2]
    
    # coevolutionary model parameters
    phi = 0.2
    alpha = 0.2
    theta = runif(n_sp, 0, 10)
    init = runif(n_sp, 0, 10)
    p = 0.1
    epsilon = 5
    eq_dif = 0.0001
    t_max = 1000
    bar = 7
    
    # simulate coevolution
    simulation = DiscoCoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max, bar)
    z_mat = simulation[[1]] # trait value matrix
    matrices = simulation[[2]] # all adjacency matrices from coevo function
    last_m = matrices[[length(matrices)]] # last adjancency matrix of the coevo process 
    last_m = as.matrix(last_m)[li, co_1:co_2] # DeSquad the last adjancency matrix
    
    net = paste(names(redes[k]), "_", a, sep = "") # name of the current empirical network
    rich = as.numeric(ncol(M)) # richness
    
    # run the control process
    dif_int = sum(init_m) - sum(last_m) # how much interactions coevolution breaks? 
    index = which(init_m == 1) # where are the interactions in the adjancency matrix?
    ints = sample(index, dif_int) # sample this interactions and break the same number
    W[ints] = 0 # of interactions that coevolution breaks
   
    # get nestedness measures
    dnest_control = nested(W, method = "NODF2") - nested(init_m, method = "NODF2")
    dnest_coevo = nested(last_m, method = "NODF2") - nested(init_m, method = "NODF2")
    
    results = data.frame(net, rich, antprob, dnest_control, dnest_coevo) # get all the results
    final_fl = rbind(final_fl, results) # put results in data.frame
    
    # save the initial adjancency matrix, the control and coevolution adjacency matrix
    write.table(init_m, file = paste(names(redes[k]), "_init_", a, ".txt", sep = ""), 
      row.names = FALSE, col.names = FALSE) # save initial adjacency matrix
    write.table(last_m, file = paste(names(redes[k]), "_final_coevo_", a, ".txt", sep = ""), 
      row.names = FALSE, col.names = FALSE) # save the coevolution adjacency matrix
    write.table(W, file = paste(names(redes[k]), "_final_control_", a, ".txt", sep = ""), 
      row.names = FALSE, col.names = FALSE) # save the control adjacency matrix
  }
}

# save or load the RData file
#save(final_fl, file = "data_nest_mod.RData")
#load("data_nest_mod.RData")

# plot and save the nestedness results graph
plot_nest_coevo = ggplot(data = final_fl) +
  geom_point(aes(x = antprob, y = dnest_coevo, colour = rich), alpha = 0.8) +
  geom_smooth(aes(x = antprob, y = dnest_coevo), colour = "red") +
  geom_smooth(aes(x = antprob, y = dnest_control), colour = "black") +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Delta Nestedness") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

ggsave(plot_nest_coevo, filename = "deltanest_coevo_adj.png", dpi = 600,
       width = 20, height = 14, units = "cm")

# read the MODULAR results
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/matrices/resultsSA/")
mod_results <- read.table("OUT_MOD.txt", header=TRUE)

# separate the initial, final control and final coevo matrices
init = mod_results[grep(mod_results$File, pattern="init"), ]
coevo = mod_results[grep(mod_results$File, pattern="coevo"), ]
control = mod_results[grep(mod_results$File, pattern="control"), ]

# separate the data and prepare the data names
# initia matrices data
init$net <- gsub(init$File, pattern = "_init", replacement = "")
init$net <- gsub(init$net, pattern = ".txt", replacement = "")

final_init <- merge(init, final_fl, by = "net", all.x=TRUE, all.y=TRUE)

# coevo matrices data
coevo$net <- gsub(coevo$File, pattern="_final_coevo", replacement = "")
coevo$net <- gsub(coevo$net, pattern=".txt", replacement = "")

final_coevo <- merge(coevo, final_fl, by="net", all.x=TRUE, all.y=TRUE)

# control matrices data
control$net <- gsub(control$File, pattern="_final_control", replacement = "")
control$net <- gsub(control$net, pattern=".txt", replacement = "")

final_control <- merge(control, final_fl, by="net", all.x=TRUE, all.y=TRUE)

# calculate the delta modularity and get the community richness
dmcontrol = final_control$Modularity - final_init$Modularity
dmcoevo = final_coevo$Modularity - final_init$Modularity
rich = final_init$rich

# final data.frame to plot the results
dados = data.frame(final_init$net, rich, final_init$antprob, dmcontrol, dmcoevo)
final_mod = dados

# save or load the RData file
save(final_mod, file = "data_mod.RData")
#load("data_mod.RData")

# plot and save the delta modularity results
plot_mod_coevo = ggplot(data = dados) +
  geom_point(aes(x = final_init.antprob, y = dmcoevo, colour = rich), alpha = 0.8) +
  geom_smooth(aes(x = final_init.antprob, y = dmcoevo), colour = "red") +
  geom_smooth(aes(x = final_init.antprob, y = dmcontrol), colour = "black") +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Delta Modularity") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

ggsave(plot_mod_coevo, filename = "deltamod_coevo_adj.png", dpi = 600,
       width = 20, height = 14, units = "cm")
