# loading packages and functions
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")

library(ggplot2)
library(cowplot)
library(dplyr)
library(NbClust)
library(rlist)
library(parallel)

# read all the empirical networks
temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

antprob_vec = seq(0.01, 1, 0.01)

p_data = data.frame() # create data.frame to allocate results
list_mats = list()
list_degree = list()

for(k in 1:length(redes)){ # loop to each empirical matrix
  print(k)

  for(a in 1:length(antprob_vec)){ # and to each antprob value
    antprob = antprob_vec[a] # define antprob
 
      for(b in 1:30){ # 100 to each antprob value
        set.seed(b)
        M = as.matrix(redes[[k]]) # M is the adjancency matrix of interactions
        M[which(M > 1)] = 1 # if there are any error, correct that
        M = SquareMatrix(M) # square the adjancency matrix
        n_sp = ncol(M) # define the species number
    
        # load functions
        source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/Antagonize.R")
        source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/CoevoMutAntNet.R")
        source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SpDegree.R")

        # insert cheaters exploitation outcomes
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
    
        # simulate coevolution
        z_mat = CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max)

        df = scale(t(z_mat))
        
        list_mats = list.append(list_mats, df)
        
        degree = SpDegree(M, V)
        list_degree = list.append(list_degree, degree)
        
        net = names(redes[k])
        rich = n_sp
        results = data.frame(net, rich, antprob)
        p_data = rbind(p_data, results)
    }
  }
}

# save or load the data created
#save(p_data, file = "p_data.RData")
save(list_mats, file = "list_mats.RData")

# create an empty document to allocate the apply results
write("clustering", "clustering.vec")

# simple function to apply the NbCluster in my results and save in clusterig.vec
camacho = function(list_mats){
  maximo = nrow(list_mats) - 1
  
  clust_an = NbClust(data = list_mats, diss = NULL, distance = "euclidean",
                     min.nc = 2, max.nc = maximo, method = "ward.D2", 
                     index = "gap")
  
  n_cl = clust_an$Best.nc[1]
  
  write(n_cl, file = "clustering.vec", append = TRUE)
  
  return(n_cl)
}

# NbCluster using 14 computer cores
cl = makeCluster(detectCores() - 2)
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

p_data = cbind(p_data, opt_clusters)


# save or load the data created
#save(p_data, file = "clustering_empirical.RData")
load(file = "clustering_empirical.RData")

#type = c(rep("Pollination", 24000), rep("Seed dispersal", 24000), rep("Ant-Plant", 24000))
p_data = cbind(p_data, type)

new_data = p_data%>%
  group_by(type, antprob)%>%
  summarise(mean = mean(opt_clusters))%>%
  as.data.frame()

plot_clusters = ggplot(data = new_data) +
  geom_point(aes(x = antprob, y = mean, colour = type), show.legend = FALSE) +
  geom_line(aes(x = antprob, y = mean, colour = type), stat="smooth",method = "auto",
            alpha = 0.8, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0,0.1)) +
  scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
  xlab("Frequency of cheaters exploitation (p)") +
  ylab("Average optimized number of species traits clusters") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(1.1, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  guides(color = guide_legend(title = "Mutualism type"))

ggsave(plot_clusters, filename = "antprob_clusters.png", dpi = 600,
       width = 16, height = 13, units = "cm", bg = "transparent")