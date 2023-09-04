# trophicSDM
These codes reproduce exactly the same results published in the paper.

Given the huge amount of data produced by these codes, we chose to store online only the summary tables containing SDMs and trophic SDMs perfomances in inferring the realised and fundamental niches, for all species, parameter settings and repetitions.
These tables are 'fund_list.rds', 'bin_realised_list.rds' and 'CI_width.rds', they are available in this GitHub directory, and correspond to the output of the step 3 ('3.Merge_all_sim.RData'). 

In order to reproduce the data, the following steps should be followed (step 0 to 3 can be skipped by loadingthe above-mentioned directories):

0: Use latin hypercube sampling to sample the 50 combinations of parameters of the Lotka-Volterra model.
This produces the file 'Parameter_set_full.txt' that is available in this GitHub directory.

1: Run simulations from theoretical models by running 1_Community_simulations.R
   In particular, this scripts creates the 5000 simulated species distribution datasets. The 'for' loop can be splitted and parallelized on a cluster to reduce computational times
   The script 1_Tools.R contains all the functions for running Community_simulations.R . 

2: Analyse the simulated datasets by running 2_tSDM_inference_sim.R
   In particular, this scripts analyses each of the 5000 simulated species distribution datasets created by the previous step.
   The 'for' loop needs be splitted and parallelized on a cluster to reduce computational times (as you should count 1h for each of the 5000 datasets)
   The script 2_Tools_tSDM_functions.R  contains unctions for running 2_tSDM_inference_sim.R

3: Merge all evaluation metrics across simulation settings and repetitions by running 3.Merge_all_sim.R
   This creates the files bin_realised_list.rds, fund_list.rds and CI_width.rds that are stored in this GitHub folder.

4: Analyse and plot the evaluation metrics for each simulation settings, across the 100 repetitions, by running 4.Figures_per_setting.R
   Creates the files summary_bin_realised.rds, summary_fund.rds that contain the summary of these evaluation metrics for each simulation setting

5: Analyse the evaluation metrics across all simulation settings and repetitions by running 5.Sensitivity_analysis.R
   This is the last necessary step of the workflow. It creates the final plots summ_bin_real_all.pdf, summ_fund_all.pdf and widthCI.pdf (corresponding to figures S2, S3 and S4) 
   and the tables all_summ_bin_real, regress_bin_real, all_summ_fund, regress_fund corresponding to tables S1, S2, S3, S4.

6: Create figures 4 and 5 by running 6.Final_figures.R

Finally, note that scripts and results for the 'Ricker' and 'GLV_Kbasal' theoretical models are available in the folder 'Other_theo_models'.