
import numpy as np
import networkx as nx

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

def find_nth_neigh_vonNeu (k,m,n,nth):
    i,j = onedim_to_twodim(k,m,n) 
    
    i_up = int(i-nth >= 0) * (i-nth) + (m-(nth-i)) * int(i-nth < 0)
    i_down = int(i+nth <= m-1) * (i+nth) + (nth - ((m-1)-i)-1) * int(i+nth > m-1)
    
    j_left = int(j-nth >= 0) * (j-nth) + (n-(nth-j)) * int(j-nth < 0)
    j_right = int(j+nth <= n-1) * (j+nth) + (nth - ((n-1)-j)-1) * int(j+nth > n-1)
        
    NB_i = np.array([i,i,i_down,i_up])
    NB_j = np.array([j_right,j_left,j,j])
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
    num_neigh = num_neigh + 1 #adding each cell as the strongest neighbor to itself
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

def find_allneigh_strongSelf_vonNeu(n,m,num_neigh): # STRONG SELF_CONNECTIVITY
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_neigh = num_neigh + 1 #adding each cell as the strongest neighbor to itself
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        temp_NB = []
        for j in range(num_neigh):
            if j == 0:
            	NB = np.array([i])
            else:
            	NB = find_nth_neigh_vonNeu(i,m,n,j)
            temp_NB.append(NB)
        neigh_all.append(temp_NB)
    return neigh_all

def find_allneigh_randomRec(n,m,num_neigh):
    d = num_neigh # in degree
    N = n*m
    G = nx.random_regular_graph(d, N)
    edges = np.array(G.edges())

    # Find connectivities for BM (because of reciprocal connectivities we have n*d/2 egedes!!!)
    neigh_all = []
    for i in range(N):
        nei_list_1 = (edges[np.where(edges[:,0]==i)[0],1]).astype(int)
        nei_list_2 = (edges[np.where(edges[:,1]==i)[0],0]).astype(int)
        nei_list = np.concatenate((nei_list_1,nei_list_2),axis=0)
        neigh_all.append([[i],nei_list])    
    return neigh_all

def find_allneigh_fullConn(n,m):
    # n and m are lattice dimensions
    # num_neigh is the maximum number of nearest neighbors(first,second,etc.)
    num_cell = n*m
    neigh_all = []
    for i in range(num_cell):
        NB = np.arange(num_cell)
        neigh_all.append([[i],NB])
    return neigh_all

def rewire(neigh_all,num_rewire):
    A = np.zeros((len(neigh_all),len(neigh_all)))
    for i in range(len(neigh_all)):
        nei_indx = neigh_all[i][0]
        A[i][nei_indx] = 1

    G = nx.from_numpy_matrix(A)
    G_rand = nx.double_edge_swap(G, nswap=num_rewire, max_tries=num_rewire*3)
    A_rand = np.asarray(nx.to_numpy_matrix(G_rand))
    
    neigh_all_rand = []
    for i in range(len(A_rand)):
        temp_NB = []
        temp_NB.append(np.array([i]))
        temp_NB.append(np.where(A_rand[i]==1)[0])
        neigh_all_rand.append(temp_NB)
        
    return neigh_all_rand

def add_random(n,m,num_neigh,num_longrange,neigh_all):
    num_neighbors = (2*num_neigh+1)**2-1
    nodes_weChange = np.random.randint(0,n*m,num_longrange)
    for i in(nodes_weChange):
        longrange = np.random.randint(0,n*m)
        while longrange == i:
                  longrange = np.random.randint(0,n*m)
        if num_neigh==1:
            edge_weDelete = np.random.randint(0,num_neighbors)
            neigh_all[i][1][edge_weDelete] = longrange
        else:
            layer_weDeletefrom = np.random.randint(1,num_neigh+1)
            edge_weDelete = np.random.randint(0,layer_weDeletefrom*8)
            neigh_all[i][layer_weDeletefrom][edge_weDelete] = longrange
    return(neigh_all)


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

        
# def update_network_state(k,m,n,s,s_ref,p,neigh_all,ref):  
#     p_notActive = 1
#     for i in range(len(p)):
#         if i==0:
#             p_notActive = p_notActive * (1-p[0])
#         else:             
#             NB = neigh_all[k][i-1]
#             active_NB = sum(s[NB]>0)
#             p_notActive = p_notActive * (1-p[i])**active_NB

#     d = np.random.rand()
#     s_new = int(s_ref[k] == 0) * int(d < (1-p_notActive)) + int(s_ref[k]>0) * (s[k])
#     return s_new    


def update_network_states_ref0 (k,m,n,s,p,neigh_all):
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

