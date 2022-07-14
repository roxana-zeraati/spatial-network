""" 
Functions for defining connectivity and computing interaction probabilities.
"""

import numpy as np

def onedim_to_twodim(k,m,n):
    i = k // m + 1 - 1
    j = k % m 
    return i,j


def twodim_to_onedim (i,j,m,n):
    i = i + 1
    j = j +1
    k = (i-1) * n + j -1 
    return k


def find_nth_neigh_general (k,m,n,nth):
    # n and m are lattice dimensions
    # k is the unit id
    # nth is radius of connectivity
    i,j = onedim_to_twodim(k,m,n)
    
    i_up_all = []
    i_down_all = []
    for ct in (np.arange(nth)+1):
        i_up = int(i-ct >= 0) * (i-ct) + (m-(ct-i)) * int(i-ct < 0)
        i_down = int(i+ct <= m-1) * (i+ct) + (ct - ((m-1)-i)-1) * int(i+ct > m-1)
        i_up_all.append(i_up)
        i_down_all.append(i_down)
        
    j_left_all = []
    j_right_all = []
    for ct in (np.arange(nth)+1):
        j_left = int(j-ct >= 0) * (j-ct) + (n-(ct-j)) * int(j-ct < 0)
        j_right = int(j+ct <= n-1) * (j+ct) + (ct - ((n-1)-j)-1) * int(j+ct > n-1)
        j_left_all.append(j_left)
        j_right_all.append(j_right)
        
    x = [i_up_all[-1]]*(2*nth+1)
    y = [i_down_all[-1]]*(2*nth+1)
    z = i_up_all[:-1] + [i] + i_down_all[:-1] 
    NB_i = np.array(x + y + z +z)
    
    xx = [j_right_all[-1]]*(2*nth-1)
    yy = [j_left_all[-1]]*(2*nth-1)
    zz = j_left_all + [j] + j_right_all 
    NB_j = np.array(zz + zz + xx + yy)
    NB = twodim_to_onedim (NB_i,NB_j,m,n)
    return NB


def find_allneigh(n,m,num_neigh):
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(num_neigh):
            NB = find_nth_neigh_general(i,m,n,j+1)
            temp_NB.append(NB)
        neigh_all.append(temp_NB)
    return neigh_all


def find_allneigh_strongSelf(n,m,num_neigh): # STRONG SELF_CONNECTIVITY
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_neigh = num_neigh + 1 #adding each cell as a neighbor to itself (self-excitation)
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(num_neigh):
            if j == 0:
            	NB = np.array([i])
            else:
            	NB = find_nth_neigh_general(i,m,n,j)
            temp_NB.append(NB)
        neigh_all.append(temp_NB)
    return neigh_all


def compute_p_hyp(num_neigh,pext,p0,sigma,self_excite):
    # Finding normalization constant
    if self_excite:
        c = sigma/(num_neigh*8+1)
    else:
        c = sigma/(num_neigh*8)
    # Computing p vector
    p = [pext]
    for i in range(1,num_neigh+1):
        if num_neigh==1:
            if self_excite:
                pi = sigma/(i*8+1)
            else:
                pi = sigma/(i*8)
        else:
            pi = c * (1/i)
        p.append(pi)
    return p

def compute_p_uniform(num_neigh,pext,p0,sigma,self_excite):
    # Compute total number of neighbors
    total_num_nei = (8 * num_neigh * (num_neigh+1))/2
    
    # Computing p vector
    p = [pext]
    if self_excite:
        sigma = sigma - p0
        p.append(p0)
        
    for i in range(1,num_neigh+1):
        pi = sigma/total_num_nei
        p.append(pi)
    return p

def compute_p_uniform_randomConn(num_neigh,pext,p0,sigma,self_excite):
    # Compute total number of neighbors
    total_num_nei = num_neigh
    
    # Computing p vector
    p = [pext]
    if self_excite:
        sigma = sigma - p0
        p.append(p0)
        
    pi = sigma/total_num_nei
    p.append(pi)
    return p

        

def update_network_states(k,m,n,s,p,neigh_all):
    pActive = 0
    for i in range(len(p)):
        if i==0:
            pActive = pActive + p[0]
        else:             
            NB = neigh_all[k][i-1]  
            active_NB = sum(s[NB]>0)
            pActive = pActive + p[i]*active_NB
    
    if pActive > 1:
        pActive = 1

    d = np.random.rand()
    s_new =  int(d < pActive)
    return s_new


# functions for the 1D network
def find_upto_nth_neigh_general_1D(k,n,nth):
    i_left_all = []
    i_right_all = []
    for ct in (np.arange(nth)+1):
        i_left = int(k-ct >= 0) * (k-ct) + (n-(ct-k)) * int(k-ct < 0)
        i_right = int(k+ct <= n-1) * (k+ct) + (ct - ((n-1)-k)-1) * int(k+ct > n-1)
        i_left_all.append(i_left)
        i_right_all.append(i_right)
    NB = i_left_all + i_right_all
    return NB




def find_allneigh_strongSelf_1D_uniform(n,num_neigh): # STRONG SELF_CONNECTIVITY
    # n is number of units
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_cell = n
    neigh_all = []
    for i in range(num_cell):       
        neigh_all.append([[i], find_upto_nth_neigh_general_1D(i,n,num_neigh)])
    return neigh_all

def update_network_states_1D(k,n,s,p,neigh_all):
    pActive = 0
    for i in range(len(p)):
        if i==0:
            pActive = pActive + p[0]
        else:             
            NB = neigh_all[k][i-1]  
            active_NB = sum(s[NB]>0)
            pActive = pActive + p[i]*active_NB
    
    if pActive > 1:
        pActive = 1

    d = np.random.rand()
    s_new =  int(d < pActive)
    return s_new