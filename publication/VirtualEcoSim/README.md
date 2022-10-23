# trophicSDM
This codes reproduce exactly the same results published in the paper.

The following steps should be followed:

1: Run simulations from theoretical models by running 1_Community_simulations.R .  The data produced by the script are available at : 10.5281/zenodo.7242258
In particular notice that this script runs one 'replication' (i.e. a set of 2550 communities). As specified inside the script, this should be run 100 times (with a different parameter 'job' specified) to obtain exactly the same results of the paper. 1_Tools.R :contains all the functions for runningn Community_simulations.R . 

2: Analyse the simulated datasets by running 2_tSDM_inference_sim. R. This script also compare the results and compute goodness of fit metrics for the realised and fundamental niche.
In particular notice that the script analyses one 'replication' only (i.e. a set of 2550 communities). It should be run 100 times, each time with a different parameter 'job' specified, so that all the 100 replications are analysed. 2_Tools_tSDM_functions.R  contains unctions for running 2_tSDM_inference_sim.R

This two steps can be time consuming, the second one in particular (around 1h for each of the 100 replications). 

3: Obtain summary figures by running 3_Compare_all_sim.R. The script compare the results of the 100 replications together and plots the results presented in the paper.

