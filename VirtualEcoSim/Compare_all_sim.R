#### Script for analysing the simulations for 100 simulations
wd = "~/Documents/Phd/Futureweb/Code/"
setwd(wd)

data_raw = "~/Documents/Phd/Futureweb/Code/MY_SCRIPTS/VirtualEcoSim/SIMlists_all_cluster/tSDM/"

realised_list = fund_list = vector(mode = "list", length = 100)

methods = c("GLV_abioticGR", "GLV_abioticKbasal", "Ricker_abioticKbasal",
            "SOI_abioticER", "SOI_abioticERbasal", "VC_abiotic")

for(i in 1:100){ # load files
  
  load(file = paste0(data_raw,"Simulations_S20L3_nEnv51_nRep50_maxBI5_job_",i,"/SIMlist_job_",i,".RData"))
  print(i)
  for(j in methods){
    print(j)
    if(j %in% names(SIMlist)){
      if(!is.null(eval(parse(text=paste0("SIMlist$", j, "$eval.realised"))))){
      eval(parse(text=paste0("realised_list[[",i,"]]$",j," = SIMlist$", j, "$eval.realised"))) 
      }
      if(!is.null(eval(parse(text=paste0("SIMlist$", j, "$eval.fund"))))){
      eval(parse(text=paste0("fund_list[[",i,"]]$",j," = SIMlist$", j, "$eval.fund"))) 
      }
    }
  }
  
  rm(SIMlist)
}

# First check the most stricking LV results (order by difference depending on the metric) and find the plot

delta_wass = delta_loo = vector(length=100)

for(i in 1:100){
  tab = realised_list[[i]]$GLV_abioticGR
  delta_wass[i] = tab[which(model == "tSDM" & metric == "wasserstein" &
                              TL !1 & CV & type == "bin"),] -
                  tab[which(model == "SDM" & metric == "wasserstein" &
                                                          TL !1 & CV),]
    
  delta_loo[i] = tab[which(model == "tSDM" & metric == "loo" &
                             TL !1 & CV & type == "bin"),] -
    tab[which(model == "SDM" & metric == "wasserstein" &
                TL !1 & CV),]
}


# Plot for all the species and all the plots the pairs