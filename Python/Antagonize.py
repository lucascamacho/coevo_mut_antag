def Antagonize(M, antprob):
		''' Function to antagonize interactions of a adjacency matriz of interactions 
		following a certain probability called antprob.

		Args:
				M: Adjacency matrix of interactions
				antprob: Probabilibty of a link become antagonist

		Return:
				M: Adjacency matrix of interactions with 0's
				V: Adjacency matrix of antagonistic interactions with 1's'''
    
		V = M * 0 # Generate an empty matrix V with equal dimension of M
    P = np.array(np.random.rand(n_sp, n_sp)) # Generate numbers from 0 to 1
    
    M[antprob >= P] = 0 # Off the link if the number generated in P are bigger than antprob
    V[antprob >= P] = 1 # On the link if the number generated in P are bigger than antprob

    return M, V # Return the matrices
