## Overview
Codes for running the network model with spatial connectivity (in 2D or 1D) to study how network structure affects spatiotemporal correlations.

For details and neuroscietific implication of this model check:   
Zeraati, R., Shi, Y., Steinmetz, N. A., Gieselmann, M. A., Thiele, A., Moore, T., Levina, A. & Engel, T. A. (2021). Intrinsic timescales in the visual cortex change with selective attention and reflect spatial connectivity. bioRxiv 2021.05.17.444537. https://www.biorxiv.org/content/10.1101/2021.05.17.444537v2.  

Analytical derivations of spatiotemporal correlations in this model are provided in:  
Shi, Y. L., Zeraati, R., Levina, A., & Engel, T. A. (2022). Spatial and temporal correlations in neural networks with structured connectivity. arXiv preprint arXiv:2207.07930. 
https://arxiv.org/abs/2207.07930.

Please cite the above two references when you use these codes for a scientific publication.


Functions/scripts descriptions:
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
