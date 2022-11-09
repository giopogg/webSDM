# trophicSDM
This codes reproduce exactly the same results published in the paper. The results from the scripts are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7242258.svg)](https://doi.org/10.5281/zenodo.7242258)



The following steps should be followed:

1: Run simulations from theoretical models by running 1_Community_simulations.R .  The data produced by the script are available at : https://doi.org/10.5281/zenodo.7242258 . Each folder contains one of the 100 replications.
In particular notice that this script runs one 'replication' (i.e. a set of 2550 communities). As specified inside the script, this should be run 100 times (with a different parameter 'job' specified) to obtain exactly the same results of the paper. The script 1_Tools.R contains all the functions for running Community_simulations.R . 

2: Analyse the simulated datasets by running 2_tSDM_inference_sim. R. This script has to be applied to each of the 100 replication. The data produced by the application of this script are the "SIMList_job_*.RData" inside each of the folder provided in https://doi.org/10.5281/zenodo.7242258 . 
In order to reproduce the results, that script has to be applied to each of the 100 datasets created by the previous step, each time with a different parameter 'job' specified. The script 2_Tools_tSDM_functions.R  contains unctions for running 2_tSDM_inference_sim.R

This two steps can be time consuming, the second one in particular (around 1h for each of the 100 replications). 

3: Obtain summary figures by running 3_Compare_all_sim.R. The script compare the results of the 100 replications together and plots the results presented in the paper. Starting from the data  provided in https://doi.org/10.5281/zenodo.7242258, the script reproduces the same figures of the paper.
