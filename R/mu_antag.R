# Code to create a matrix, transform a certain proportion of links in antagonists and simulate coevolution
# Lucas Camacho (8604710)

# create matrix
M <- matrix(1, ncol = 3, nrow = 3)
diag(M) <- 0
n_sp <- dim(M)[1]

# define args
max_t <- 1000
init <- runif(n_sp, min = 0, max = 10)
theta <- runif(n_sp, min = 0, max = 10)
intfor <- 0.9
alpha <- 0.2
antprob <- 0.2
A <- M * 0
phi <- 0.2
epsilon <- 2

### change links
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    if(M[i,j] == 1){
      lucky <- runif(1, min = 0, max = 1)
      if(s <= antprob){
        M[i,j] = 0
        A[i,j] = 1
      }
    }  
  }
}

# sum the matrices to get the binary matrix of interactions
Q = M + A

# begin dynamics simulation
#forças de interação
results <- matrix(NA, ncol = n_sp)
results[1,] <- init

for(t in 1:max_t){
  for(i in 1:n_sp){
    for(j in 1:n_sp){
      if(Q[i,j] >= 0){
#        Q[i,j] = exp(-alpha * ((z[j] - z[i])^2))
      }
    }
  }
  
  Q = Q / apply(Q, 1, sum)
  Q = intfor * Q
  S = matrix(0, ncol = n_sp, nrow = 1)
  
  for(i in 1:n_sp){
    S[i] = S[i] + (phi * ((1 - intfor) * (theta[i] - z[i])))
    for(j in 1:n_sp){
      if(M[i,j] > 0){
        S[i] = S[i] + (phi * (Q[i,j] * (z[j] - z[i])))}
      else{
        if(A[i,j] > 0){
          if(abs(z[j] - z[i]) < epsilon){
            if(z[i] > z[j]){
              S[i] = S[i] + (phi * (q[i,j] * (z[j] + barreira - z[i])))}
            else{
              S[i] = S[i] + (phi * (q[i,j] * (z[j] - barreira - z[i])))}
          }
        }
      }
    }
  }
  
  z <- S + z
  resultados <- rbind(resultados, z)
}

traits = as.data.frame(resultados)
n_sp = ncol(traits)
traits_vec = c(as.matrix(traits))
traits_df = data.frame(species = rep(paste("sp", 1:n_sp, sep = ""), each = nrow(traits)),
                       time = rep(1:nrow(traits), times = n_sp),
                       trait = traits_vec)
# plotting traits through time
p = ggplot(traits_df, aes(x = time, y = trait, color = species)) +
  geom_point(size = 1.8, shape = 19, alpha = 0.6) +
  ggtitle(paste("proportion antagonists = ", antprob)) +
  xlab("Time") +
  ylab("Mean species trait (z)") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12))

print(p)