# trophicSDM
This codes reproduce exactly the same results of the theoretical models 'Ricker' and 'GLV_Kbasal' published in the paper.

The results from the scripts are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7309691.svg)](https://doi.org/10.5281/zenodo.7309691) . The main Zenodo repository contains 100 folders, each corresponding to one replications. The folder contains the results from the theoretical model simulations corresponding to step 1 (the .csv files) and a SIMList_job_*.RData file that contains the results of their statistical analyses.

In order to reproduce the data using the following steps should be followed:

1: Run simulations from theoretical models by running 1_Community_simulations.R .
In particular notice that this script runs one 'replication' (i.e. a set of 2550 communities). As specified inside the script, this should be run 100 times (with a different parameter 'job' specified) to obtain exactly the same results of the paper. The script 1_Tools.R contains all the functions for running Community_simulations.R . 

2: Analyse the simulated datasets by running 2_tSDM_inference_sim. R. This script has to be applied to each of the 100 replication.
In order to reproduce the results, that script has to be applied to each of the 100 datasets created by the previous step, each time with a different parameter 'job' specified. The script 2_Tools_tSDM_functions.R  contains unctions for running 2_tSDM_inference_sim.R

This two steps can be time consuming, the second one in particular (around 1h for each of the 100 replications). 

3: Obtain summary figures by running 3_Compare_all_sim.R. The script compare the results of the 100 replications together and plots the results presented in the paper. Starting from the results of the statistical analyses step 2 (i.e. the SIMList_job_*.RData provided in https://doi.org/10.5281/zenodo.7242258), the script reproduces the same figures of the paper.
