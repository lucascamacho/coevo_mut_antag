def CoevoMutAntNet(n_sp, M, V, phi, alpha, theta, init, p, epsilon, eq_dif, t_max):
    '''Simulates the coevolutionary dynamics of mutualists and antagonists in a network
        Args:
            n_sp: total number of species in the network
            M: square adjacency matrix representing the mutualistic interactions
            V: square adjacency matrix representing the antagonistic interactions
            phi: vector of phi values (additive genetic variance * heritability)
            alpha: alpha value (sensitivity of selection to trait matching)
            theta: vector of environmental optimum values
            init: vector of initial trait values
            p: vector of strength of selection due to the environment
            epsilon: barrier value for antagonistic interactions
            eq_dif: value to determine when equilibrium is reached
            t_max: maximum number of timesteps allowed

       Obs: 
            All vectors need to have first the row species attributes and then 
            the column species attributes (e.g. c(row_sp[1], ..., row_sp[nrow],
                                                  col_sp[1], ..., col_sp[ncol]))

       Returns:
            A matrix containing, in each row t, the trait values (z) of all 
            species at time t'''
    # load modules
    import numpy as np
    import pandas as pd
    # generate an matrix to allocate the species traits
    z_mat = np.zeros((t_max, n_sp))
    z_mat[0,] = init
    #for loop to change the traits in time
    for r in range(t_max):
        z = z_mat[r,]
        A = M + V
        z_dif = np.transpose(A * z) - (A * z)
        Q = A * (np.exp(-alpha * (z_dif ** 2)))
        np.fill_diagonal(Q, 0)
        Q_n = Q / Q.apply(np.sum, axis = 1)
