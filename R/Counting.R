Counting = function(M, V){
  # Count the frequency of Antagonists (AA), Cheaters (AM) and Mutualists (MM) in a network
  #
  # Args:
  # M: Adjacency matrix with mutualistic interactions
  # V: Adjacency matrix with antagonistic interactions
  #
  # Obs: M and V are complementar matrix. The sum of these matrices must give all the interactions in the network
  #
  # Return:
  # A list with the frequency of AA, AM and MM.
  
  n_cheat = (length(M[M != t(M)]) / 2) #counting cheaters

  m1 = which(M == 1)
  m2 = which(t(M) == 1)
  n_mut = sum(table(m1[m1 %in% m2])) / 2 #counting mutualism

  v1 = which(V == 1)
  v2 = which(t(V) == 1)
  n_inter = sum(table(v1[v1 %in% v2])) / 2 # counting antagonism

  freq = list(n_inter, n_cheat, n_mut) #create and return a list with the frequencies
  return(freq)

}