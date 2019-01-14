def End_Interaction(M, V, interaction):
    if(interaction != "antagonism" 
       and interaction != "antagonism" and interaction != "antagonism"):
        print ("Choose a valid interactions type")
    
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
        index = M == np.transpose(M) # end mutualism based on index of M and M transpose
        M[index] = 0
        return M, V
