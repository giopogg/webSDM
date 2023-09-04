##############################################################################################################
########### Script for creating the model parameters for the sensitivity analyses of the Lotka Volterra model.

########################################################################################################################
########################################################################################################################
# Script for creating the model parameters for the sensitivity analyses of the Lotka Volterra model.

# !!! Only the working directory has be specified to run the script

# Running the script will the file parameter_set_full.txt containing all simulation parameters that will be used in the pipeline

############################################################

# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

library(lhs)
library(dplyr)

# Specify directory
dir = ""

# define all the possible values of model parameters 
S = c(10, 20, 50)
p = c(0.2, 0.35, 0.5)
L = c(3,4)
niche_breadthGRratio = c(2, 6, 10)
nEnv = c(30,50,70)
nRep = c(30, 50, 70)

params = list(S = S,p = p, L = L, niche_breadthGRratio = niche_breadthGRratio, nEnv = nEnv, nRep = nRep)

# Sample in the latin hypercube
set.seed(123)

X = improvedLHS(50, length(params))

# Assign to the corresponding categorical variables
parameter_set = matrix(NA, ncol = length(params), nrow = nrow(X), dimnames = list(NULL, names(params)))

for(i in 1:length(params)){
  
  param = params[[i]]
  intervals = seq(0, 1, length.out = length(param) +1)
  parameter_set[,i] = sapply(1:nrow(X),
                    function(n) param[which(sapply(1:length(S),
                                               function(s) between(X[n,i], intervals[s],
                                                                   intervals[s + 1])))]
                    )
  
}

# For each parameter setting add 100 repetitions
parameter_set_full = parameter_set[rep(1:nrow(parameter_set), each = 100),]
parameter_set_full = data.frame(parameter_set_full, job = rep(1:100, nrow(parameter_set)))


## Each line corresponds to a repetition for a given parameter setting.
# The column 'job' corresponds to the repetition index. All other columns are the parameter setting
write.table(parameter_set_full, paste0(dir,"Parameter_set_full.txt"), row.names = FALSE)

