[![DOI](https://zenodo.org/badge/366038263.svg)](https://zenodo.org/badge/latestdoi/366038263)


## Overview
Codes for running the network model with spatial connectivity (in 2D or 1D) to study how network structure affects spatiotemporal correlations.

For model details and the implications of these models to study the neural dynamics in visual cortex check:   
Zeraati, R., Shi, Y., Steinmetz, N. A., Gieselmann, M. A., Thiele, A., Moore, T., Levina, A. & Engel, T. A. (2023). Intrinsic timescales in the visual cortex change with selective attention and reflect spatial connectivity. Nature communications, 14(1), 1858. https://doi.org/10.1038/s41467-023-37613-7.  

Analytical derivations of spatiotemporal correlations in this model are provided in:  
Shi, Y.L., Zeraati, R., Levina, A. and Engel, T.A. (2023). Spatial and temporal correlations in neural networks with structured connectivity. Physical Review Research, 5(1), p.013005.
https://doi.org/10.1103/PhysRevResearch.5.013005.

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


Additional network models:
- Network with random connectivity (conn_type = 'random_hetro') and two different cell types: act_gen defined in activity_generator_hetro.py
- Network with synaptic timescales: act_gen defined in activity_generator_wSynFilter.py  


## Dependencies
- Python >= 3.7.1
- Numpy >= 1.15.4 
