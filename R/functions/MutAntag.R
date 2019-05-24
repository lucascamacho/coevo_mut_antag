#-----------------------------------------------------------------------------------------------------#
MutAntag = function(M, V, r, prob_change){
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
    ##############
    # MM to AM shift
    # identify MM outcomes
    inter = (M == 1) == (t(M) == 1)
    inter[M == 0] = FALSE
    mm_posit = which(inter == TRUE)
    # there are MM outcomes? If yes...
    if(length(mm_posit) != 0){
      # there's more than 1 MM outcome? If yes...
      ch = sample(mm_posit, 1) # sample 1 outcome and change to AM.
      M[ch] = 0
      V[ch] = 1
    }
    ##############
    # AM to MM shift
    # identify AM outcomes
    am_posit = which(M != t(M))
    zero = M[am_posit] == 0
    c_trues = table(zero)["TRUE"]
    # there are AM outcomes? If yes...
    if(c_trues != 0){
      # there's more than 1 AM outcomes? If yes...
      if(c_trues > 1){
        # sample 1 outcomes and change to MM
        ch = sample(am_posit[zero], 1)
        M[ch] = 1
        V[ch] = 0
      }
      else{ #else, use the available AM outcomes
        ch = am_posit[zero]
        M[ch] = 1
        V[ch] = 0
      }
    }
    #return the results
    w_time = r
    mats = list(M, V, w_time)
    return(mats)
}
  else{
    w_time = NULL
    mats = list(M, V, w_time)
    return(mats)
  }
}  
#-----------------------------------------------------------------------------------------------------#