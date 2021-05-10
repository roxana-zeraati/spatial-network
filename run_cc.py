"""
This script runs the simulation code for given parameters and computes cross-correlations between units on a certain distance from each other.
"""

import activity_generator as act
import network_setup as net_set
import numpy as np
import random

L = 100 # lattice size as L*L
m = 0.99 # branching parameter (ps + 8pr)
ps = 0.88 # self-excitation probability
pext = 0.0001 # external input
ptype = 'local' # connectivity structure
num_trials = 1 # number of trials
T = 100000 # duration of trials (in time steps)
timelag = 300 # maximum time-lag for computing cross-correlations

N = L*L # total number of units
num_nei_cc = 4 # number of units from each distance
min_dist = 1 # minimum distance
max_dist = 7 # maximum distance
jump_dist = 2 # separation between distances

# directories for saving the results
ac_save_path = '../cc/'
sim_save_path = '../sim/'


for tr in range(num_trials):
    print('tr:', tr)

    # simulate the network
    laT = act.act_gen(L,m,pext,ps,ptype,T)
#     np.save(sim_save_path +ptype+'_ps'+str(pa)+'_sig'+str(m)+'_pext'+str(pext)+'_T'+str(T), laT)


    # compute cross-correlations
    for neigh_dist in range(min_dist,max_dist+1,jump_dist):
        print('distance:', neigh_dist)
        crosscor_sum = 0
        nonzero_cells = 0
        neigh_all = net_set.find_allneigh(L,L,neigh_dist) # find neigbors up to dist without self-excitation
        total_neigh = len(neigh_all[0][neigh_dist-1])

        # cc with units on a certian distance
        for cell_num in range(0,N,1):
#             print(cell_num)
            s1 = laT[cell_num]
 
            if np.sum(s1)>0:        
                # picking a random neighbor
                idx = random.sample(range(total_neigh), num_nei_cc)
                for i in range(len(idx)):        
                    neigh_id = neigh_all[cell_num][neigh_dist-1][idx[i]]
                    sn = laT[neigh_id]
                    if np.sum(sn)>0:  
                        crosscor = []
                        for lag in range(timelag):
                            s2 = np.concatenate((np.zeros(lag),sn[:T-lag]))
                            cc = (np.dot(s1,s2)/(T-lag) - np.mean(s1[lag:])*np.mean(sn[:T-lag]))/(np.std(s1)*np.std(sn))
                            crosscor.append(cc)
                        crosscor = np.asarray(crosscor)
                        if np.sum(np.isnan(crosscor)) < 1:
                            crosscor_sum = crosscor_sum + crosscor
                            nonzero_cells = nonzero_cells +1

        crosscor_avg = crosscor_sum/nonzero_cells
    #     np.save(cc_save_path +'CC_'+ptype+'_ps'+str(ps)+'_bp'+str(m)+'_pext'+str(pext)+'_dist'+str(neigh_dist)+'_T'+str(T)+ '_lag'+ str(timelag),crosscor_avg)


