########################################################################################################################
########################################################################################################################
# Data analysis script from Giovanni Poggiato and Jérémy Andréoletti
# This scripts analysis all the 5000 species distribution datasets obtained using the script 1_Community_simulations.R

# !!! Only the working directory and the parameter job should be specified to run the script

# By running this script, it will enrich the folders created by 1_Community_simulations.R

# !!!! The for loop is extremely long to run. We advice to split the for loop and run each of the 5000 iterations 
# separately on a cluster

############################################################


# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

rm(list=ls())

library(igraph)
library(Matrix)
library(GGally)
library(intergraph)
library(gridExtra)
library(dismo)
library(coda)
library(transport)
library(dplyr)
library(rstan)
library(rstanarm)
library(arm)
library(glmnet)
library(formula.tools)
library(abind)
library(brms)
library(stringr)
library(bayesplot)
library(parallel)
library(cheddar)
library(ggpubr)
# ! If windows
# Parallel processing with windows
# library(devtools)
# install_github('nathanvan/parallelsugar')
# library(parallelsugar)


source("2.Tools_tSDM_inference.R") # Script for running functions. The documentation could be outdated, please refer to the documentation of the R package webSDM.


# !!!!!!! Specify model parameters

## Read parameter settings
set_param_full = read.table("Parameter_set_full.txt", header = TRUE)


# Analyse the simulated species distribution dataset. Each dataset corresponds to a line of the set_param_full_table.
# !! This is very long but can be fully parallelized
# To fully reproduce all the results, we advice to split the following loop, and run
# each iteration (or a given set of iteration) separately on a cluster
for(i in 1:nrow(set_param_full)){
  
  
  job = set_param_full$job[i]
  S = set_param_full$S[i]
  p = set_param_full$p[i]
  L = set_param_full$L[i]
  nicheBreadthRatio = set_param_full$niche_breadthGRratio[i]
  nEnv = set_param_full$nEnv[i]
  nRep = set_param_full$nRep[i]
  
  simPath = paste0(dir, "SimulationOutputs/Simulations_S_",S,"_L",L,"_p",p,
                   "_nicheBreadthRatio",nicheBreadthRatio,"_nEnv",nEnv,"_nRep",nRep,"_job",job,"/")
  
  if(!dir.exists(simPath)){ 
    
    cat("Folder does not exist, stop! \n")
    
  }else{
    
    
    # the function runTrophicSDMAnalysis can return an error (and prints its motivations), for example 
    # if at least one species was never or was always present in the dataset
    try(runTrophicSDMAnalysis(S = S, L = L, p = p, niche_breadthGR = nicheBreadthRatio/S,
                              nEnv = nEnv, nRep = nRep, simPath = simPath))
    
  }
  
  if(i %% 100 == 0) print(paste0(i,"\n"))
  
}
