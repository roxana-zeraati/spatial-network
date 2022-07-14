"""
This script runs the simulation code for given parameters and computes autocorrelations.
"""

import activity_generator as act
import numpy as np

L = 100 # lattice size as L*L
m = 0.99 # branching parameter (ps + 8pr)
ps = 0.88 # self-excitation probability
pext = 0.0001 # external input
conn_type = 'local' # connectivity structure (other connectivity types now only work with L = 100 size)
num_trials = 1 # number of trials
T = 100000 # duration of trials (in time steps)
timelag = 300 # maximum time-lag for computing autocorrelations

N = L*L # total number of units

# directories for saving the results
ac_save_path = '../ac/'
sim_save_path = '../sim/'


ac_all = [] # autocorrelations 
psth_all = [] # trial firing rate over time
for tr in range(num_trials):
    print('tr:', tr)

    # simulate the network, default is R=1 (connectivity radius), otherwise also define R
    laT = act.act_gen(L,m,pext,ps,conn_type,T)
#     np.save(sim_save_path +conn_type+'_ps'+str(pa)+'_sig'+str(m)+'_pext'+str(pext)+'_T'+str(T), laT)


    # compute autocorrelations
    autocor_sum = 0
    nonzero_cells = 0
    for cell_num in range(0,N,10): # average autocorrelation over every 10th unit
        s1 = laT[cell_num]
        if np.sum(s1)<1:
            continue
        autocor = []
        for lag in range(timelag):
            s2 = np.concatenate((np.zeros(lag),s1[:T-lag]))
            ac = (np.dot(s1,s2)/(T-lag) - np.mean(s1[lag:])*np.mean(s1[:T-lag]))
            autocor.append(ac)
        autocor = np.asarray(autocor)
        if np.sum(np.isnan(autocor)) < 1:
            autocor_sum = autocor_sum + autocor
            nonzero_cells = nonzero_cells +1
    autocor_avg = autocor_sum/nonzero_cells
    ac_all.append(autocor_avg)    
    psth_all.append(np.mean(laT, axis =0))
# np.save(ac_save_path +'ac_' + conn_type +'_T' + str(T)+'_trial'+str(num_trials)+'_ps'+str(ps)+'_bp'+str(m) +'_pext'+ str(pext),[psth_all, ac_all])


