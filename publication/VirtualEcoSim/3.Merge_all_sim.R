########################################################################################################################
########################################################################################################################
# Script for merging simulation results throughout all simulation settings from Giovanni Poggiato

# This scripts merges all the evaluation metrics from the analysis all the 5000 species distribution datasets obtained using the script 2_tSDM_inference.R, into few final lists

# !!! Only the working directory should be specified to run the script

# By running this script, it will create bin_realised_list, fund_list and CI_width that are available in the GitHub
############################################################


# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


library(dplyr)
library(tidyr)

rm(list = ls())

dir = ""

data_raw = paste0(dir,"/SimulationOutputs/")
fig = paste0(dir,"/SimulationOutputs/Fig_all")

###############################################################################################################################################################################
# Assemble simulation results in one list


realised_list = fund_list = CI_width = data.frame()

dirs = list.dirs(recursive = FALSE, full.names = FALSE)

for(i in dirs){ # load files
  
  # If SIMlist exists and can be loaded
  if(file.exists(paste0(data_raw, i,"/SIMlist.RData"))) {
    if( class(try(load(file = paste0(data_raw, i,"/SIMlist.RData")))) != "try-error"){
      
      niche_optima = read.csv2(file = paste0(data_raw,
                                             i,"/niche_optima.csv"))[,1]
      
      G = read.csv2(file = paste0(data_raw,
                                  i,"/InteractionMatrix.csv"))[,-1]
      
      rownames(G) = colnames(G)
      
      S = as.numeric(sub("_L.*","",sub(".*S_","", i)))
      L = as.numeric(sub("_p.*","",sub(".*_L","", i)))
      p = as.numeric(sub("_niche.*","",sub(".*_p","", i)))
      nbLink = length(which(G[upper.tri(G)]!=0))
      meanConnectivity = 2*nbLink/(S*(S-1)) #nb of link over all possible links
      meanStrength = mean(G[upper.tri(G)][which(G[upper.tri(G)]!=0)]) #mean non-zero values (decreases with S, increases with p)
      meanStrengthS = meanStrength/S #Should provide an info that combine both meanStrength and nbLink
      nicheBreadthRatio = as.numeric(sub("_nEnv.*","",sub(".*Ratio","", i)))
      nEnv = as.numeric(sub("_nRep.*","",sub(".*_nEnv","", i)))
      nRep = as.numeric(sub("_job.*","",sub(".*_nRep","", i)))
      job = as.numeric(sub(".*_job","", i))
      
      #Evaluation metrics for realised niche
      if(!is.null(SIMlist$eval.realised)){
        
        realised_list = rbind(realised_list,
                              data.frame(cbind(SIMlist$eval.realised,
                                               source = i,
                                               S = S,
                                               L = L,
                                               p = p,
                                               nbLink = nbLink,
                                               meanConnectivity = meanConnectivity,
                                               meanStrength = meanStrength,
                                               meanStrenghS = meanStrengthS,
                                               nicheBreadthRatio = nicheBreadthRatio,
                                               nEnv = nEnv,
                                               nRep = nRep,
                                               job = job
                              )))
        
      }
      
      #Evaluation metrics for fundamental niche
      if(!is.null(SIMlist$eval.fund)){
        
        fund_list = rbind(fund_list,
                          data.frame(cbind(SIMlist$eval.fund,
                                           source = i,
                                           S = S,
                                           L = L,
                                           p = p,
                                           nbLink = nbLink,
                                           meanConnectivity = meanConnectivity,
                                           meanStrength = meanStrength,
                                           meanStrenghS = meanStrengthS,
                                           nicheBreadthRatio = nicheBreadthRatio,
                                           nEnv = nEnv,
                                           nRep = nRep,
                                           job = job
                          )))
        
      }
      
      
      # Width CI
      CI_width = rbind(CI_width,
                       data.frame(cbind(SIMlist$eval.widthCI,
                                        source = i,
                                        S = S,
                                        L = L,
                                        p = p,
                                        nbLink = nbLink,
                                        meanConnectivity = meanConnectivity,
                                        meanStrength = meanStrength,
                                        meanStrenghS = meanStrengthS,
                                        nicheBreadthRatio = nicheBreadthRatio,
                                        nEnv = nEnv,
                                        nRep = nRep,
                                        job = job
                       )))
      
      rm(SIMlist)
      if(job == 100) cat("### ",i, "### \n")
    }
  }
}




saveRDS(realised_list,file=paste0(data_raw,"realised_list.rds"))

## Spread table and compute difference: tSDM - SDM
bin_realised_list = realised_list %>% filter(type != 'prob') %>%
                                      select(-type) %>%
                                      spread(key = model, value =value)
bin_realised_list$diff = bin_realised_list$tSDM - bin_realised_list$SDM


fund_list = fund_list  %>% spread(key = model, value =value)
fund_list$diff = fund_list$tSDM - fund_list$SDM


saveRDS(bin_realised_list,file=paste0(data_raw,"bin_realised_list.rds"))
saveRDS(fund_list,file=paste0(data_raw,"fund_list.rds"))
saveRDS(CI_width, file=paste0(data_raw,"CI_width.rds"))
