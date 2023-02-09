"""
Network with two diferent cell types
This script reads the connectivity structure and simulates network activity.
"""

import numpy as np
from network_setup import *
import random


connectivity_path  = './connectivity_structures/'



def act_gen(L,m,pext,ps,conn_type,T):
    """function for simulating network activity.

    Parameters
    -----------
    L : int
        lattice size as L*L.
    m : float
        branching parameter.
    ps : float, array
        self-excitation probability
        for network with two cell types: different probabilities with a mixing coefficient, [ps1, ps2, c1]
        default for c1 = 0.5
    conn_type: string
        connectivity structure ('local', 'random', 'random_sp2', 'random_sp3', 'random_sp5', 'random_sp7')
        'local' is Moore neighborhood, 'random_spR' is 8 randomly selected neighbors within the radius R.
        Connectivity types other than 'local' now only work with L = 100 size.
    T : int
        number of time-steps per each trial
    

    Returns
    -------
    laT : 2d array
        network activity as #units * #time-steps

    """
    
    
    num_neigh = 1 # this is for having 8 neighnors 
    self_exciteP = 1 # cosidering self-exciation        
    
   
    if conn_type == 'random_hetro':
        if len(ps) == 3:
            ps1 = ps[0]
            ps2 = ps[1]
            c1 = ps[2]
        else:
            ps1 = ps[0]
            ps2 = ps[1]
            c1 = 0.5

        pr1 = (m-ps1)/((2*num_neigh+1)**2-1)
        pr2 = (m-ps2)/((2*num_neigh+1)**2-1)
        ps_hetro = np.concatenate((np.ones(int(L*L*c1))*ps1, np.ones(int(L*L*(1-c1)))*ps2))
        pr_hetro = np.concatenate((np.ones(int(L*L*c1))*pr1, np.ones(int(L*L*(1-c1)))*pr2))
        p = [pext, ps_hetro, pr_hetro] 
        neigh_all = np.load(connectivity_path + 'neigh_random_8con_size10000.npy', allow_pickle = True)

    else:
        raise ValueError('The connectivity strcuture is undefined.')


    print('probabilities: ',p)
    do_break = 0
    lattice_activity = []
    save_counter = 0
    do_break = 0

    #initial condition
    average_active = int((pext/(1-m)) * (L*L))
    num_cells = L*L
    print('initial cond (#active units): ', str(average_active))
    id_initial_active = random.sample(range(num_cells), average_active)
    s = np.zeros(L*L)
    s[id_initial_active] = 1


    for i in range(T+1):
        if do_break:
            break
        lattice_activity.append(s)
        save_counter = save_counter+1     	    

        s = np.array([update_network_states_hetrogen(k,L,L,s,p,neigh_all) for k in range(0,len(s))])   

        if save_counter == T:
            do_break = 1
            break
    avg_fr = np.sum(lattice_activity)/(T)
    print('avg_fr: ', avg_fr)
    return np.transpose(lattice_activity)
