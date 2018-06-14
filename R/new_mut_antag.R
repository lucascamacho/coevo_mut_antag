#library(ggplot2)

#
n_sp = 3
M = matrix(1, ncol = n_sp, nrow = n_sp)
diag(M) = 0
t_max = 200
z = runif(n_sp, 0, 10)
theta = runif(n_sp, 0, 10)
intfor = 0.1
alpha = 0.2
resultados = matrix(NA, ncol = n_sp, nrow = 1)
resultados[1,] = z
antprob = 0.2
phi = 0.2
V = M * 0
barreira = 5

#

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

#

for(r in 1:t_max){
  #
  for(i in 1:n_sp){
    for(j in 1:n_sp){
      if(Q[i,j] > 0){
        Q[i,j] = G[i,j] * exp(-alpha * ((z[j] - z[i]) ** 2))
      }
    }
  }

  #
  Q = Q / apply(Q, 1, sum)
  Q = intfor * Q
  S = matrix(0, nrow = 1, ncol = n_sp)
  
  #
  for(i in 1:n_sp){
    S[i] = S[i] + phi * (1 - intfor) * (theta[i] - z[i])
    for(j in 1:n_sp){
      if(M[i,j] > 0){
        S[i] = S[i] + phi * Q[i,j] * (z[j] - z[i])
      }
      else{
        if(V[i,j] > 0){
          if(abs(z[j] - z[i]) < barreira){
            if(z[i] < z[j]){
              S[i] = S[i] + phi * Q[i,j] * (z[j] + barreira - z[i])
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
  z = z + S
  resultados = rbind(resultados, z)
  
}

#matplot(resultados, type = "l", lty = 1, lwd = 2)

# building data frame to use in ggplot
traits = as.data.frame(resultados)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                       time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)
# plotting traits through time
plotar = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
  geom_point(size = 1.8, shape = 19, alpha = 0.6) + 
  ggtitle(paste("proportion antagonists = ", antprob)) +
  xlab("Time") + 
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))

print(plotar)