MutualizeAntagonize = function(M, V, r){
  # probabilidade de ocorrência p e q
#  p = runif(1, 0, 1) # probability of generate mutualisms
  q = runif(1, 0, 1) # probability of generate antagonisms

  if(q <= 0.25 & any(M == TRUE) == TRUE){ # 0, 0.1, 0.10
    w_t_vm = NULL
    w_t_mv = r
    change_m = which(M == 1)
    int_change_m = sample(change_m, 1)
    M[int_change_m] = 0
    V[int_change_m] = 1
    
    change_v = which(V == 1)
    int_change_v = sample(change_v, 1)
    V[int_change_v] = 0
    M[int_change_v] = 1
    
    mats = list(M, V, w_t_vm, w_t_mv) # create list with the antagonistic and mutualistic matrices
    return(mats) # return a list with the antagonistic and mutualistic matrices
  }
  else{
    w_t_vm = NULL
    w_t_mv = NULL
    mats = list(M, V, w_t_vm, w_t_mv) # create list with the antagonistic and mutualistic matrices
    return(mats) # return a list with the antagonistic and mutualistic matrices
  }
}  