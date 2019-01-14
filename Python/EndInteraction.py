def EndInteraction(M, V, interaction):
    ''' Function to end some type of interaction in adjacency matrix M and V.
    If you choose "antagonism" the elements where V are equal to V transpose
    gonna be changed to zero. The same process goes to "mutualism" and
    cheaters.
    Args:
        M: Adjacency matrix of interactions with 0's
        V: Adjacency matrix of antagonistic interactions with 1's
        interaction: which interaction you wanna end
    Return:
        M and V matrices with a certain interaction ended'''

    # load modules
    import numpy as np

    # message error to a non possible interaction ex: competition
    if(interaction != "antagonism" 
        and interaction != "antagonism" and interaction != "antagonism"):
            print ("Choose a valid interactions type")

    # check and delete the interactions where the element is equal to his
    # transpose. The process is the same for all the interactions choosen
    if(interaction == "antagonism"):
        index = V == np.transpose(V)
        V[index] = 0
        return M, V

    if(interaction == "cheaters"):
        index = V != np.transpose(V)
        V[index] = 0
        M[index] = 0
        return M, V

    else:
        index = M == np.transpose(M) 
        M[index] = 0
        return M, V