#-----------------------------------------------------------------------------------------------------#
MutualizeAntagonize = function(M, V, r, prob_change){
  # Follow the probability in prob_change, transform a random mutualism outcome in antagonism,
  # also, transform a antagonist outcome in mutualism. This function tracks in which
  # timestep occurs the interactions shifts.
  #
  # obs: if you use this function in a loop, the w_time object is usefull.
  #
  # Args:
  #  M: matrix of mutualistic outcomes
  #  V: matrix of antagonistic outcomes
  #  r: timestep of simulation
  #  prob_change: probability that in each timestep the interaction shift between M and V
  #
  # Return:
  #  M: new matrix of mutualistic outcomes
  #  V: new matrix of antagonistic outcomes
  #  w_time: timestep in which the outcomes change. 
  #
  p = runif(1, 0, 1) # sample a number between 0 and 1
  
  if(p <= prob_change){ # outcomes shifts will happend
    # identify MM outcome and transform into AM or MA
    inter = (M == 1) == (t(M) == 1)
    inter[M == 0] = FALSE
    mm_posit = which(inter == TRUE)
    ch = sample(mm_posit, 1)
    M[ch] = 0
    V[ch] = 1
    
    # identify AM or MA outcome and transform into MM
    am_posit = which(M != t(M))
    zero = M[am_posit] == 0
    c_trues = table(zero)["TRUE"]
    
    # choose a antagonist outcome if there are more than one...
    if((c_trues > 1) == TRUE){
      mut = sample(am_posit[zero], 1)
      M[mut] = 1
      V[mut] = 0
      
      #return the results
      w_time = r
      mats = list(M, V, w_time)
      return(mats)
    }
    else{ # [...]or use the only avaliable antagonist outcome 
      mut = am_posit[zero]
      M[mut] = 1
      V[mut] = 0
      
      #return the results
      w_time = r
      mats = list(M, V, w_time)
      return(mats)
    }
    
  }
  else{ # outcome shift will not happend
    # return the original matrices
    w_time = NULL
    mats = list(M, V, w_time)
    return(mats)
  }
}  
#-----------------------------------------------------------------------------------------------------#