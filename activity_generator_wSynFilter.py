"""
This script reads the connectivity structure and simulates network activity.
"""

import numpy as np
from network_setup import *
import random

connectivity_path  = './connectivity_structures/'

def act_gen_synFilt(L,m,pext,ps,conn_type,T, t_s):
    """function for simulating network activity.

    Parameters
    -----------
    L : int
        lattice size as L*L.
    m : float
        branching parameter.
    ps : float
        self-excitation probability
    conn_type: string
        connectivity structure ('local', 'random', 'random_sp2', 'random_sp3', 'random_sp5', 'random_sp7')
        'local' is Moore neighborhood, 'random_spR' is 8 randomly selected neighbors within the radius R.
        Connectivity types other than 'local' now only work with L = 100 size.
    T : int
        number of time-steps per each trial
    t_s: float
        timescale of synaptic filter
    

    Returns
    -------
    laT : 2d array
        network activity as #units * #time-steps

    """
    
    
    num_neigh = 1 # this is for having 8 neighnors 
    self_exciteP = 1 # cosidering self-exciation        
    
   
    if conn_type == 'local':
        neigh_all = find_allneigh_strongSelf(L,L,num_neigh) 
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
    elif conn_type == 'random':
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
        neigh_all = np.load(connectivity_path + 'neigh_random_8con_size10000.npy', allow_pickle=True) 
    elif conn_type == 'random_sp2':
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
        if num_neigh == 1:
                neigh_all = np.load(connectivity_path + 'neigh_randomSpa_8con_size10000_k2.npy', allow_pickle=True) 
    elif conn_type == 'random_sp3':
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
        if num_neigh == 1:
                neigh_all = np.load(connectivity_path + 'neigh_randomSpa_8con_size10000_k3.npy', allow_pickle=True) 
    elif conn_type == 'random_sp5':
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
        if num_neigh == 1:
                neigh_all = np.load(connectivity_path + 'neigh_randomSpa_8con_size10000_k5.npy', allow_pickle=True)                   
    elif conn_type == 'random_sp7':
        p = [pext, ps, (m-ps)/((2*num_neigh+1)**2-1)] 
        if num_neigh == 1:
                neigh_all = np.load(connectivity_path + 'neigh_randomSpa_8con_size10000_k7.npy', allow_pickle=True)
    else:
        raise ValueError('The connectivity strcuture is undefined.')
   

    print('probabilities: ',p)
    do_break = 0
    lattice_activity = []
    filt_activity = []
    save_counter = 0
    do_break = 0

    #initial condition
    average_active = int((pext/(1-m)) * (L*L))
    average_active_f = int((pext/(1-(m-ps))) * (L*L))
    num_cells = L*L
    print('initial cond (#active units): ', str(average_active))
    id_initial_active = random.sample(range(num_cells), average_active)
    id_initial_active_f = random.sample(range(num_cells), average_active_f)
    s = np.zeros(L*L)
    f = np.zeros(L*L)
    s[id_initial_active] = 1
    f[id_initial_active_f] = 1


    for i in range(T+1):
        if do_break:
            break
        lattice_activity.append(s)
        filt_activity.append(f)
        save_counter = save_counter+1         
        # update synaptic filter
        f_new = np.array([update_synFilter(k,L,L,s,p,neigh_all,f,t_s) for k in range(0,len(f))])
        # update unit states
        s = np.array([update_network_states_synFilter(k,L,L,s,p,neigh_all,f) for k in range(0,len(s))])   
        f = f_new

        if save_counter == T:
            do_break = 1
            break
    avg_fr = np.sum(lattice_activity)/(T)
    print('avg_fr: ', avg_fr)
    return np.transpose(lattice_activity), np.transpose(filt_activity)
