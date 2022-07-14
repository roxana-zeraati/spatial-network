## Overview
Codes for running the network model with spatial connectivity (in 2D or 1D) to study how network structure affects spatiotemporal correlations.

For details check:   
Zeraati, R., Shi, Y., Steinmetz, N. A., Gieselmann, M. A., Thiele, A., Moore, T., Levina, A. & Engel, T. A. Attentional modulation of intrinsic timescales in visual cortex and spatial networks. bioRxiv 2021.05.17.444537 (2021). https://www.biorxiv.org/content/10.1101/2021.05.17.444537v1.  
Please cite this reference when you use these codes for a scientific publication.


- act_gen: the function that simulates the network (import from activity_generator.py)
- run_ac.py: the script for simulating the network and computing autocorrelations.
- run_cc.py: the script for simulating the network and computing cross-correlations.


Available connectivity structures for the network (conn_type variable):
- 'local': 2D network with local connectivity (defined within the radius R in Chebyshev distances, R=1 is the Moore neighborhood)
- 'local_1D': 1D network with local connectivity (defined within the radius R in Chebyshev distances)
- 'random': 2D network with random connectivity
- 'random_spR' (e.g., 'random_sp2'): 2D network with 8 random connections within the radius R = 2,3,5,7


## Dependencies
- Python >= 3.7.1
- Numpy >= 1.15.4 
