# definindo args
matriz <- matrix(1, ncol = 3, nrow = 3)
diag(matriz) <- 0

tempo <- 300
z <- runif(dim(matriz)[1], min = 0, max = 10)
theta <- runif(dim(matriz)[1], min = 0, max = 10)
intfor <- 0.9
alpha <- 0.2
antprob <- 0.2
vmatriz <- matriz * 0
phi <- 0.2
barreira <- 2

### transformando links
for(i in 1:nrow(matriz)){
  for(j in 1:ncol(matriz)){
    if(matriz[i,j] == 1){
      sorte <- runif(1, min = 0, max = 1)
      if(sorte <= antprob){
        matriz[i,j] = 0
        vmatriz[i,j] = 1
      }
    }  
  }
}

#somando a matriz de mutualismo e antagonismo
q = matriz + vmatriz

#dinamica coevolutiva
#forças de interação
resultados <- matrix(NA, ncol = dim(matriz)[1])
resultados[1,] <- z

for(t in 1:tempo){
  for(i in 1:dim(matriz)[1]){
    for(j in 1:dim(matriz)[1]){
      if(q[i,j] >= 0){
        q[i,j] = exp(-alpha * ((z[j] - z[i])^2))
      }
    }
  }
  
  q = q /apply(q,1,sum)
  q = intfor * q
  s = matrix(0, ncol = dim(matriz)[1], nrow = 1)
  
  for(i in 1:dim(matriz)[1]){
    s[i] = s[i] + (phi * ((1 - intfor) * (theta[i] - z[i])))
    for(j in 1:ncol(matriz)){
      if(matriz[i,j] > 0){
        s[i] = s[i] + (phi * (q[i,j] * (z[j] - z[i])))}
      else{
        if(vmatriz[i,j] > 0){
          if(abs(z[j] - z[i]) < barreira){
            if(z[i] > z[j]){
              s[i] = s[i] + (phi * (q[i,j] * (z[j] + barreira - z[i])))}
            else{
              s[i] = s[i] + (phi * (q[i,j] * (z[j] - barreira - z[i])))}
          }
        }
      }
    }
  }
  
  z <- s + z
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