## Overview
Codes for running the network model with spatial connectivity to study how network structure affects local and global temporal dynamics.

For details check: preprint link
Please cite this reference when you use this code for a scientific publication.

- act_gen: the function that simulates the network (import from activity_generator.py)
- run_ac.py: the script for simulating the network and computing autocorrelations.
- run_cc.py: the script for simulating the network and computing cross-correlations.

Available connectivity structures for the network (conn_type variable):
- 'local': Moore neighborhood
- 'random': random connectivity
- 'random_spR' (e.g., 'random_sp2'): 8 random connections within the radius R


## Dependencies
- Python >= 3.7.1
- Numpy >= 1.15.4 
