def ZeroLines(M, V, n_sp, antprob):
    '''Function to check for zero lines in the adjacency matrices.
    If there are some zeroed lines, the Coevolution simulation will stop
    Args:
        M: Adjacency matrix of interactions with 0's
        V: Adjacency matrix of antagonistic interactions with 1's
        n_sp: number of species in the matrix
    antprob: Probabilibty of a link become antagonist
    Return
        new M and V matrices with non-zero lines'''

    import numpy as np
    from Antagonize import Antagonize
    from EndInteraction import EndInteraction

    for z in range(0, 1000000):

        m = M.sum(axis = 1)

        if 0 in m:
            # create M and V matrices
            M = np.ones((n_sp, n_sp))
            np.fill_diagonal(M, 0)
            # antagonize M and generate V matrix
            antagonize = Antagonize(M, antprob)
            M = antagonize[0]
            V = antagonize[1]

            # end the pure antagonism
            end = EndInteraction(M, V, "antagonism")
            M = end[0]
            V = end[1]

        else:
            break

    return M, V