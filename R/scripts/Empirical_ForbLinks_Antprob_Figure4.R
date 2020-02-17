# Last figure of the paper
# Code to explore the influence of cheaters exploitation in the structure
# of the adjancency matrix of interactions. We are considering a second
# trait barrier (b) that will cut interactions off if this barrier is tranpassed
# This will allow us to see what's the effect of cheaters in coevolution and
# how this change's the strcuture of ecological interacions in the community.
#
# The nestedness will be calculated using the bipartite package and the modularity
# with the MODULAR program
#
# Marquitti, F. M. D., P. R. Guimaraes, M. M. Pires, L. F. Bittencourt. 2014. 
# MODULAR: Software for the Autonomous Computation of Modularity in Large Network Sets. 
# Ecography: 37: 221–224.
#
# For each of our 24 empirical networks we gonna run 
# 3.000 simulations and calculate the nestedness and modularity
# of out adjancency matrix of interactions.

# load packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(bipartite)
library(dplyr)

# read all mutualism networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

# create data.frame to store all my results
antprob_vec = seq(0.01, 1, 0.01)
final_fl = data.frame()

for(k in 1:length(redes)){ # loop to each matrix of interactions
  print(k)
  
  for(a in 1:length(antprob_vec)){ # 100 loops to each matrix
    
    for(q in 1:30){ # 30 simulations per p value
    M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
    M[which(M > 1)] = 1 # if there are any error, correct that
    M = SquareMatrix(M) # square the adjancency matrix
    n_sp = ncol(M) # define the species number

    # sample an antprob value
    antprob = antprob_vec[a]
    
    # load functions
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/DiscoCoevoMutAntNet.R")
    source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Counting.R")
    
    # antagonize matrices (now I have an M and V adjancency matrix)
    antagonize = Antagonize(M, antprob)
    M = antagonize[[1]]
    V = antagonize[[2]]
    
    # DesQuad matrices
    li = 1:dim(redes[[k]])[1]
    co_1 = dim(redes[[k]])[1] + 1
    co_2 = dim(redes[[k]])[1] + dim(redes[[k]])[2]
    
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
    #dnest_control = nested(W, method = "NODF2") - nested(init_m, method = "NODF2")
    #dnest_coevo = nested(last_m, method = "NODF2") - nested(init_m, method = "NODF2")
    
    #results = data.frame(net, rich, antprob, dnest_control, dnest_coevo) # get all the results
    #final_fl = rbind(final_fl, results) # put results in data.frame
    
    # save the initial adjancency matrix, the control and coevolution adjacency matrix
    write.table(init_m, file = paste("~/Google Drive File Stream/Meu Drive/Trabalho/matrices/", 
                                     names(redes[k]), "_init_", a,"_", q, ".txt", sep = ""), row.names = FALSE, 
                                     col.names = FALSE) # save initial adjacency matrix
    write.table(last_m, file = paste("~/Google Drive File Stream/Meu Drive/Trabalho/matrices/", 
                                     names(redes[k]), "_final_coevo_", a,"_", q, ".txt", sep = ""), row.names = FALSE, 
                                     col.names = FALSE) # save the coevolution adjacency matrix
    write.table(W, file = paste("~/Google Drive File Stream/Meu Drive/Trabalho/matrices/", 
                                names(redes[k]), "_final_control_", a,"_", q, ".txt", sep = ""), row.names = FALSE, 
                                col.names = FALSE) # save the control adjacency matrix
    }
  }
}

# save or load the RData file
save(final_fl, file = "data_nest.RData")
#load("data_nest.RData")

# aggregate the data using the net and get the measures average
final_fl = aggregate(final_fl[ ,4:5], list(final_fl$net), mean)

# create and insert a antprob sequence
antprob = rep(seq(0.01, 1, 0.01), 24)
final_fl = cbind(final_fl, antprob)

# create and insert an mutualism type sequence
type = c(rep("Pollination", 800), rep("Seed dispersal", 800), rep("Ant-Plant", 800))
final_fl = cbind(final_fl, type)

# forming the data frame to plot, calculating averages
new_data = final_fl%>%
  group_by(antprob, type)%>%
  summarise(mean_nestcontrol = mean(dnest_control), mean_nestcoevo = mean(dnest_coevo))%>%
  as.data.frame()

# plot and save the nestedness results graph
plot_nest_coevo = ggplot(data = new_data) +
  geom_point(aes(x = antprob, y = mean_nestcoevo, colour = type, group = type), show.legend = FALSE, alpha = 0.8) +
  geom_line(aes(x = antprob, y = mean_nestcoevo, colour = type), stat="smooth",method = "lm",
            alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean_nestcontrol, colour = type, group = type), stat= "smooth", method = "lm",
            alpha = 0.8, show.legend = FALSE, linetype = "dashed") +
  scale_color_brewer(palette="Dark2") +
  scale_y_continuous(expand = c(0,1)) +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Delta Nestedness") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

ggsave(plot_nest_coevo, filename = "deltanest_coevo_adj.png", dpi = 600,
       width = 16, height = 12, units = "cm",  bg = "transparent")

# Run MODULAR

# read the MODULAR results
mod_results = read.table("~/Google Drive File Stream/Meu Drive/Trabalho/matrices/resultsSA/OUT_MOD.txt", 
                         header=TRUE)

# separate the initial, final control and final coevo matrices
init = mod_results[grep(mod_results$File, pattern = "init"), ]
coevo = mod_results[grep(mod_results$File, pattern = "coevo"), ]
control = mod_results[grep(mod_results$File, pattern = "control"), ]

# separate the data and prepare the data names
# initia matrices data
init$net = gsub(init$File, pattern = "_init", replacement = "")
init$net = gsub(init$net, pattern = ".txt", replacement = "")

final_init = merge(init, final_fl, by = "net", all.x = TRUE, all.y = TRUE)

# coevo matrices data
coevo$net = gsub(coevo$File, pattern = "_final_coevo", replacement = "")
coevo$net = gsub(coevo$net, pattern = ".txt", replacement = "")

final_coevo = merge(coevo, final_fl, by="net", all.x=TRUE, all.y=TRUE)

# control matrices data
control$net = gsub(control$File, pattern = "_final_control", replacement = "")
control$net = gsub(control$net, pattern = ".txt", replacement = "")

final_control = merge(control, final_fl, by = "net", all.x = TRUE, all.y = TRUE)

# calculate the delta modularity and get the community richness
dmcontrol = final_control$Modularity - final_init$Modularity
dmcoevo = final_coevo$Modularity - final_init$Modularity
rich = final_init$rich

# Rearranje data to plot the results
# final data.frame to plot the results
dados = data.frame(final_init$net, rich, final_init$antprob, dmcontrol, dmcoevo)
final_mod = dados

# save or load the RData file
save(final_mod, file = "data_mod.RData")
#load("data_mod.RData")

# check modularity

novo = final_mod
novo[,1]=gsub('[0-9]+', "", novo[,1])

novo %>% group_by(final_init.net) %>% arrange(final_init.antprob, .by_group = TRUE)


novo = as.data.frame(novo)
final_mod = aggregate(novo[ ,4:5], list(novo$final_init.net), mean)


# create and insert an antprob sequence
antprob = rep(seq(0.01, 1, 0.01), 24)
final_fl = cbind(final_fl, antprob)

# create and insert an mutualism type sequence
type = c(rep("Pollination", 800), rep("Seed dispersal", 800), rep("Ant-Plant", 800))
final_fl = cbind(final_fl, type)

# prepare data frame to plot using the averages
new_data = final_mod%>%
  group_by(antprob, type)%>%
  summarise(mean_modcontrol = mean(dmod_control), mean_modcoevo = mean(dmod_coevo))%>%
  as.data.frame()

#plot and save the modularity data
plot_mod_coevo = ggplot(data = new_data) +
  geom_point(aes(x = final_init.antprob, y = mod_coevo, colour = type, group = type), show.legend = FALSE) +
  geom_line(aes(x = final_init.antprob, y = mod_coevo, colour = type), stat="smooth",method = "lm",
            alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(x = final_init.antprob, y = mod_control, colour = type, group = type), stat= "smooth", method = "lm",
            alpha = 0.8, show.legend = FALSE, linetype = "dashed") +
  scale_color_brewer(palette="Dark2") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  xlab("") +
  ylab("Delta Modularity") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 20), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

ggsave(plot_mod_coevo, filename = "deltamod_coevo_adj.png", dpi = 600,
       width = 16, height = 12, units = "cm",  bg = "transparent")