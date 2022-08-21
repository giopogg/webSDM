# trophicSDM
This codes reproduce exactly the same results published in the paper. 
The following steps should be followed:
# 1: Run simulations from theoretical models by running 1_Community_simulations.R
In particular notice that this script runs one simulation 'replication'. As specified inside the script this sould be run 100 times (with a different parameter 'job' specified) to obtain exactly the same results of the paper. 1_Tolls.R :contains all the functions for running Community_simulations.R

# 2: Analyse datasets by running 2_tSDM_inference_sim. R. This script also compare the results and compute goodness of fit metrics for the realised and fundamental niche.
In particular notice that the script analyses one 'replication' only. It should be run 100 times, each time with a different parameter 'job' specified, so that all the 100 replications are analysed. 2_Tools_tSDM_functions.R  contains unctions for running 2_tSDM_inference_sim.R

This two steps can be time consuming, the second one in particular (for each repetition it takes around 10 minutes to run the first script, and around 1h for the second one). We advice to run only one 'repetition', and to eventually decrease the number of iteration of the MCMC sampling algorithm, to have a first glance to the results. Then, in the aim to fully reproduce the study, we advice to use an external cluster. 

# 3: Obtain summary figures by running 3_Compare_all_sim.R. The script compare the results of 100 runs of simulations together and plots the results presented in the paper.

