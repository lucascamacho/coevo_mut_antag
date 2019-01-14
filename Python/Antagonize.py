def Antagonize(M, antprob):
    ''' Function to antagonize interactions of a adjacency matriz of 
    interactions following a certain probability called antprob.
    Args:
    M: Adjacency matrix of interactions
    antprob: Probabilibty of a link become antagonist
    Return:
    M: Adjacency matrix of interactions with 0's
    V: Adjacency matrix of antagonistic interactions with 1's'''
    
    # import modules
    import numpy as np
    # Generate an empty matrix V with equal dimension of M
    V = M * 0
    # Generate numbers from 0 to 1
    P = np.array(np.random.rand(len(M), len(M)))
    # Off the link if the number generated in P are bigger than antprob
    M[antprob >= P] = 0
    # On the link if the number generated in P are bigger than antprob
    V[antprob >= P] = 1 
    # matrix diagonal must be zero
    np.fill_diagonal(V, 0)
    # Return the matrices
    return M, V