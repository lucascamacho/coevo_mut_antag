# Sensibility analyses for phi
#

# load packages/functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/MeanPairDist.R")

# simple function to apply the NbCluster in my results and save in clusterig.vec
# create an empty document to allocate the apply results
write("clustering", "clustering.vec")
camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  write(n_cl, file = "clustering.vec", append = TRUE)
  
  return(n_cl)
}

library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)

# read all the empirical networks
temp = list.files(pattern="*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement= "", temp)

antprob_vec = seq(0.1, 1, 0.1)

p_data = data.frame()
list_mats = list()
list_degree = list()


for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)
  
  for(a in 1:length(antprob_vec)){ # loop to each antprob (p) value
    
    for(q in 1:10){ # 10 simulations per e seq
      M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
      M[which(M > 1)] = 1 # if there are any error, correct that
      M = SquareMatrix(M) # square the adjancency matrix
      n_sp = ncol(M) # define the species number
      
      antprob = antprob_vec[a] # get the antprob value
      
      # load functions
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
      source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
      
      # insert cheaters exploitation outcomes
      empantagonize = Antagonize(M, antprob)
      M = empantagonize[[1]]
      V = empantagonize[[2]]
      
      # coevolutionary model parameters
      phi_seq = seq(0.1, 1, 0.1)
      phi = phi_seq[q]
      alpha = 0.2
      theta = runif(n_sp, 0, 10)
      init = runif(n_sp, 0, 10)
      p = 0.1 #chose the p value
      epsilon = 5
      eq_dif = 0.0001
      t_max = 1000
      
      # simulate coevolution
      z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)
      
      #df = scale(t(z_mat))
      #list_mats = list.append(list_mats, df)
      
      # get my results
      net = names(redes[k])
      rich = as.numeric(ncol(M))
      mpd = MeanPairDist(z_mat[nrow(z_mat), ])
      
      # insert these results in the data.frame
      #epsilon = unique(epsilon)
      results = data.frame(net, rich, antprob, phi, mpd)
      p_data = rbind(p_data, results)
    }
  }
}

p_data$mpd = round(p_data$mpd, digits = 2)
#p_data = p_data[!duplicated(p_data), ]

# save or load the data created
save(p_data, file = "Major_sensibility_phi.RData")
#save(list_mats, file = "~/Google Drive File Stream/Meu Drive/Trabalho/list_mats.RData")
#load("antprob_var.RData")

# NbCluster using 14 computer cores
cl = makeCluster(detectCores())
clusterEvalQ(cl, {
  library(NbClust)
  camacho = function(list_mats){
    maximo = nrow(list_mats) - 1
    
    clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                       min.nc = 2, max.nc = maximo, method = "ward.D2", 
                       index = "gap")
    
    n_cl = clust_an$Best.nc[1]
    
    write(n_cl, file = "clustering.vec", append = TRUE)
    
    return(n_cl)
  }
})
opt_clusters = parallel::parSapply(cl, list_mats, camacho)
stopCluster(cl)

# bind the clustering results with our data table
p_data = cbind(p_data, opt_clusters)

# save or load the results
save(p_data, file = "Major_sensibility_phi.RData")
#load("Major_sensibility_e.RData")

new_data = p_data%>%
  group_by(antprob, phi)%>%
  summarise(mean_mpd = mean(mpd),
            mean_clusters = mean(opt_clusters))%>%
  as.data.frame()

plot_mpd = ggplot(data = p_data, aes(x = antprob, y = mpd)) +
  geom_point(aes(color = phi, fill = phi), alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Average MPD") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_clusters = ggplot(data = p_data, aes(x = antprob, y = opt_clusters)) +
  geom_point(aes(color = phi, fill = phi), alpha = 0.5, size = 2,) +
  geom_line(aes(x = antprob, y = opt_clusters), stat = "smooth", method = "loess", 
            show.legend = FALSE) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Averge Number of Clusters of Species Traits") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()


plot_mpd = ggplot(data = new_data, aes(x = antprob, y = mean_mpd)) +
  geom_point(aes(color = phi, fill = phi), alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Average MPD") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_clusters = ggplot(data = new_data, aes(x = antprob, y = mean_clusters)) +
  geom_point(aes(color = phi, fill = phi), alpha = 0.5,size = 2,) +
  xlab("Frequency of cheaters interactions (p)") + ylab("Averge Number of Clusters of Species Traits") +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  theme_classic()

plot_final = grid.arrange(plot_mpd, plot_clusters, nrow = 2)

ggsave(plot_final, filename = "Major_Sensibility_phi.png", dpi = 600,
       width = 15, height = 15, units = "cm",  bg = "transparent")
