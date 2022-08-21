# trophicSDM
# This codes reproduce exactly the same results published in the paper. Notice that this can be time consuming and we advice to run everything on a cluster, following the advices of every script. We remark that running the whole simulation settings (with 100 replications) could be time consumer. Otherwise, it is possible to run only 1 replications by only running script 1 and 2 once. 

The following steps should be followed:
# 1: Run simulations from theoretical models by running 1_Community_simulations.R
In particular notice that:
1_Community_simulations.R : runs one simulation 'replication'. As specified inside the script this sould be run 100 times (with a different parameter 'job' specified) to obtain exactly the same results of the paper.
1_Tolls.R : functions for running Community_simulations.R

# 2: Analyse datasets by running 2_tSDM_inference_sim. R. This script also compare the results and compute goodness of fit metrics for the realised and fundamental niche.
In particular
2_tSDM_inference_sim. R : Analyse one 'replication'. It should be run 100 times, each time with a different parameter 'job' specified. This parameter tells which of the 100 replications should be analysed. 
2_Tools_tSDM_functions.R : Functions for running 2_tSDM_inference_sim.R

# Obtain summary measures by running:
3_Compare_all_sim.R : The script compare the results of 100 runs of simulations together and plots the results presented in the paper.

