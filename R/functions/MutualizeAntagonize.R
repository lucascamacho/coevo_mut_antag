#-----------------------------------------------------------------------------------------------------#
MutualizeAntagonize = function(M, V, r){
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
  #
  # Return:
  #  M: new adjancency matrix of mutualistic interactions
  #  V: new adjancency matrix of antagonistic interactions
  #  w_time: timestep in which the interaction change. 
  #
  p = runif(1, 0, 1) # probability of changing interactions

  if(p <= 0.01 & any(M == TRUE) == TRUE){ # if the interaction change
    w_time = r # change a mutualism to antagonism
    change_m = which(M == 1)
    int_change_m = sample(change_m, 1)
    M[int_change_m] = 0
    V[int_change_m] = 1
    
    change_v = which(V == 1) # change a antagonism to a mutualism
    if(change_v == change_m){ # don`t change the same interaction changed before`
      change_v = which(V == 1)
    }
    
    int_change_v = sample(change_v, 1)
    V[int_change_v] = 0
    M[int_change_v] = 1
    
    mats = list(M, V, w_time) # create list with results
    return(mats) # return the results
  }
  else{ # if the interaction doesn't change
    w_time = NULL
    mats = list(M, V, w_time) # create list with the antagonistic and mutualistic matrices
    return(mats) # return a list with the antagonistic and mutualistic matrices
  }
}
#-----------------------------------------------------------------------------------------------------#