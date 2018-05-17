# definindo args
matriz <- matrix(1, ncol = 3, nrow = 3)
diag(matriz) <- 0

tempo <- 10
z <- runif(dim(matriz)[1], min = 0, max = 10)
theta <- runif(dim(matriz)[1], min = 0, max = 10)
intfor <- 0.1
alpha <- 0.2
antprob <- 0.1
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
    for(i in 1:dim(matriz)[1]){
      if(q[i,j] > 0){
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
  
  resultados <- rbind(resultados, s)
}

matplot(resultados)
