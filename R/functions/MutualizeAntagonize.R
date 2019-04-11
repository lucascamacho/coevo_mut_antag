#-----------------------------------------------------------------------------------------------------#
MutualizeAntagonize = function(M, V, r, prob_change){
  # In each timestep of a simulation, transform a random mutualism link in antagonism,
  # also, transform a antagonist link in mutualism. These function tracks in which
  # timestep occurs the interactions shifts.
  #
  # obs: if you use this function in a loop, the w_time object is usefull.
  #
  # Args:
  #  M: adjancency matrix of mutualistic interactions
  #  V: adjacency matrix of antagonistic interactions
  #  r: timestep of simulation
  #  prob_change: Probability that in each timestep the interaction shift between M and V
  #
  # Return:
  #  M: new adjancency matrix of mutualistic interactions
  #  V: new adjancency matrix of antagonistic interactions
  #  w_time: timestep in which the interaction change. 
  #
  p = runif(1, 0, 1) # probability of changing interactions
  
  if(p <= prob_change){ # if the interaction change
     w_time = r # change a mutualism to antagonism
     
     j = runif(1, 0, 1)
     
     if(j <= 0.5){
       dex = sample(1:dim(M)[1], 2)
         if(M[dex[1], dex[2]] == 1){
           M[dex[1], dex[2]] == 0
           diag(M) = 0
           mats = list(M, V, w_time) # create list with results
           return(mats) # return the results
         }
       else{
         M[dex[1], dex[2]] == 1
         diag(M) = 0
         mats = list(M, V, w_time) # create list with results
         return(mats) # return the results
       }
     }
     else{
       dex = sample(1:dim(V)[1], 2)
       if(V[dex[1], dex[2]] == 1){
         V[dex[1], dex[2]] == 0
         diag(V) = 0
         mats = list(M, V, w_time) # create list with results
         return(mats) # return the results
       }
       else{
         V[dex[1], dex[2]] == 1
         diag(V) = 0
         mats = list(M, V, w_time) # create list with results
         return(mats) # return the results
       }
     }
  }   
#    change_m = which(M == 1)
#    int_change_m = sample(change_m, 1)
#    M[int_change_m] = 0
#    V[int_change_m] = 1
#    
#    change_v = which(V == 1) # change a antagonism to a mutualism
#    int_change_v = sample(change_v, 1) # choose an V interaction to change
#    if(int_change_v == int_change_m){ # don`t change the same interaction changed before
#      int_change_v = sample(change_v, 1)
#    }
#    
#    V[int_change_v] = 0
#    M[int_change_v] = 1
#    
#    mats = list(M, V, w_time) # create list with results
#    return(mats) # return the results
#  }
  else{ # if the interaction doesn't change
    w_time = NULL
    mats = list(M, V, w_time) # create list with the antagonistic and mutualistic matrices
    return(mats) # return a list with the antagonistic and mutualistic matrices
  }
     
}
#-----------------------------------------------------------------------------------------------------#