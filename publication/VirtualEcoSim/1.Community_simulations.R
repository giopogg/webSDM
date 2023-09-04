########################################################################################################################
########################################################################################################################
# Community simulation script from Giovanni Poggiato and Jérémy Andréoletti
# By running this script once, you simulate the 5000 species distribution datasets 
# (i.e., 100 repetition for each simulation settings)

# !!! Only the working directory has be specified to run the script

# Running the script will create 5000 directory with all the simulated communities and model parameters inside

# The for loop can be fully splitted in several scripts and run on a cluster
############################################################

# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


library("igraph")
library("gtools")
library("cheddar")
library("devtools") 
library("reshape2")
#library("gplots")
library("magrittr")
library("purrr")
library("readr")
library("matrixcalc")
library("vegan")
library("parallel")
library("truncnorm")
# ! If windows
# Parallel processing with windows
# library(devtools)
# install_github('nathanvan/parallelsugar')
# library(parallelsugar)

rm(list=ls())

source("1.Tools_Community_simulations.R")


# Specify directory
dir = ""

dir.create(paste0(dir, "/SimulationOutputs/"))

## Read parameter settings
set_param_full = read.table("Parameter_set_full.txt", header = TRUE)
## Each line corresponds to a repetition for a given parameter setting.
# The column 'job' corresponds to the repetition index. All other columns are the parameter setting



# Simulate Lokta-volterra communities for each line of the set_param_full_table.
# !! This can be very long but can be fully parallelized
# To fully reproduce all the results, we advice to split the following loop, and run
# each iteration (or a given set of iteration) separately on a cluster
for(i in 1:nrow(set_param_full)){

# Name the directory in which to store the simulations
simPath = paste0(dir,"/SimulationOutputs/Simulations_S_",set_param_full$S[i],"_L",set_param_full$L[i],
                 "_p",set_param_full$p[i],
                 "_nicheBreadthRatio",set_param_full$niche_breadthGRratio[i],
                 "_nEnv",set_param_full$nEnv[i],"_nRep",set_param_full$nRep[i],"_job",
                 set_param_full$job[i],"/")

dir.create(simPath, showWarnings = FALSE)

set.seed(set_param_full$job[i])

S = set_param_full$S[i]

runGRSimulation(S = set_param_full$S[i], L = set_param_full$L[i], p = set_param_full$p[i],
                niche_breadthGR = set_param_full$niche_breadthGRratio[i]/S, nEnv = set_param_full$nEnv[i],
                           nRep = set_param_full$nRep[i], simPath = simPath)

if(i %% 100 == 0) print(paste0(i,"\n"))
}

