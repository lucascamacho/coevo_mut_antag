  # This is a scratch script to simulate coevolution with exploitation and mutualism interactions.
  #
  # Please, check the original matlab code mutantag.m file in this same directory.
  #
  # This was made in the Workshop of Ecological Networks that happened in April, 2018.
  
  # loading packages
  setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/")
  library(ggplot2)
  library(cowplot)
  
  # initial conditions
  n_sp = 10
  M = matrix(1, ncol = n_sp, nrow = n_sp)
  diag(M) = 0
  t_max = 1000
  z = runif(n_sp, 0, 10)
  theta = runif(n_sp, 0, 10)
  intfor = 0.9
  alpha = 0.2
  resultados = matrix(NA, ncol = n_sp, nrow = 1)
  resultados[1,] = z
  antprob = 0.2
  phi = 0.2
  V = M * 0
  barreira = 5
  
  # Antagonize M (transform positive links in negative)
  for(i in 1:n_sp){
    for(j in 1:n_sp){
      if(M[i,j] == 1){
        l = runif(1, 0, 1)
        if(l <= antprob){
          M[i,j] = 0
          V[i,j] = 1
        }
      }
    }
  }
  
  # 
  Q = M + V
  G = Q
  
  # Calculate Q matrix
  for(r in 1:t_max){
    #
    for(i in 1:n_sp){
      for(j in 1:n_sp){
        if(Q[i,j] > 0){
          Q[i,j] = G[i,j] * exp(-alpha * ((z[j] - z[i]) ** 2))
        }
      }
    }
  
    # Normalize Q matrix
    Q = Q / apply(Q, 1, sum)
    Q = intfor * Q
    S = matrix(0, nrow = 1, ncol = n_sp)
    
    # Calculate selection differencials (SD)
    for(i in 1:n_sp){
      S[i] = S[i] + phi * (1 - intfor) * (theta[i] - z[i]) # environment SD
      for(j in 1:n_sp){
        if(M[i,j] > 0){
          S[i] = S[i] + phi * Q[i,j] * (z[j] - z[i]) # mutualism SD
        }
        else{
          if(V[i,j] > 0){
            if(abs(z[j] - z[i]) < barreira){
                if(z[i] < z[j]){
                  S[i] = S[i] + phi * Q[i,j] * (z[j] + barreira - z[i]) # antagonism SD
                }
                else{
                  S[i] = S[i] + phi * Q[i,j] * (z[j] - barreira - z[i])
              }
            }
          }
        }
      }
    }
    
    #
    z = z + S # update the trait values
    resultados = rbind(resultados, z)
    
  }
  
  # prepare data.frame
  traits = as.data.frame(resultados)
  n_sp = ncol(traits)
  traits_vec = c(as.matrix(traits))
  traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                         time = rep(1:nrow(traits), times = n_sp),
                         trait = traits_vec)
  
  # plot the results
  plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
    geom_path(size = 1.8, alpha = 0.8) + 
    ggtitle(paste("proportion exploitation = ", antprob)) +
    xlab("Time") + 
    ylab("Mean species trait (z)") +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 14), 
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size = 12))
  
  plotar

# save the plot
ggsave(plotar, filename = "Slow_Basic_Traits.pdf", dpi = 600,
       width = 12, height = 24, units = "cm",  bg = "transparent")