#### Script for analysing the simulations for 100 simulations
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

rm(list = ls())

wd = "~/Documents/Phd/Futureweb/Code/"
setwd(wd)

data_raw = "~/Documents/Phd/Futureweb/Code/MY_SCRIPTS/VirtualEcoSim/SIMlists_all_cluster/tSDM/"
fig = "~/Documents/Phd/Futureweb/Code/MY_SCRIPTS/VirtualEcoSim/SIMlists_all_cluster/Fig_all"
load_list = T

methods = c("GLV_abioticGR", "GLV_abioticKbasal", "Ricker_abioticKbasal",
            "SOI_abioticER", "SOI_abioticERbasal", "VC_abiotic")


###############################################################################################################################################################################
# Assemble simulation results in one list

if(!load_list){
realised_list = fund_list = BioControl_wass = BioControl_cor =vector(mode = "list", length = 100)

sd_methods = c(0.3, 0.05, 0.05, 0.3, 0.3, 0.3)
names(sd_methods) = methods
S=20

for(i in 1:100){ # load files
  
  load(file = paste0(data_raw,"Simulations_S20L3_nEnv51_nRep50_maxBI5_job_",
                     i,"/SIMlist_job_",i,".RData"))
  
  niche_optima = read.csv2(file = paste0(data_raw,
                                         "Simulations_S20L3_nEnv51_nRep50_maxBI5_job_",
                                         i,"/niche_optima.csv"))[,1]
  
  print(i)
  for(j in methods){
    print(j)
    if(j %in% names(SIMlist)){
      
      if(!is.null(SIMlist[[j]]$eval.realised)){
        realised_list[[i]][[j]] = cbind(SIMlist[[j]]$eval.realised,
                                        source = i,
                                        method = j)
      }
      
      if(!is.null(SIMlist[[j]]$eval.fund)){
        fund_list[[i]][[j]] = cbind(SIMlist[[j]]$eval.fund,
                                        source = i,
                                        method = j)
      }
      
      # Compute interaction strength:

      #1) Wasserstein distance between realised niche and growth rate
      BioControl_wass[[i]][[j]] = data.frame( value = sapply(1:S,
                               function(s)transport::wasserstein1d(
                                 SIMlist[[j]]$survivalProba[,s],
                                 dnorm(seq(0,1,length.out=51),
                                       mean = niche_optima[s],
                                       sd = sd_methods[[j]]) )) ,
                               species = colnames(SIMlist[[j]]$survivalProba),
                               source = i,
                               methods = j)

      # 2) Correlation between realised niche and growth rate
      BioControl_cor[[i]][[j]] = data.frame(value = sapply(1:S,
                              function(s) cor(SIMlist[[j]]$survivalProba[,s], 
                                              dnorm(seq(0,1,length.out=51),
                                                    mean = niche_optima[s],
                                                    sd = sd_methods[[j]]) )),
                              species = colnames(SIMlist[[j]]$survivalProba),
                              source = i,
                              methods = j)
    }
  }

  rm(SIMlist)
}

saveRDS(realised_list,file=paste0(data_raw,"realised_list.rds"))
saveRDS(fund_list,file=paste0(data_raw,"fund_list.rds"))

saveRDS(BioControl_cor,file=paste0(data_raw,"BioControl_cor.rds"))
saveRDS(BioControl_wass,file=paste0(data_raw,"BioControl_wass.rds"))

}else{
  
realised_list = readRDS(file=paste0(data_raw,"realised_list.rds"))
fund_list = readRDS(file=paste0(data_raw,"fund_list.rds"))

BioControl_cor = readRDS(file=paste0(data_raw,"BioControl_cor.rds"))
BioControl_wass = readRDS(file=paste0(data_raw,"BioControl_wass.rds"))

}

###############################################################################################################################################################################
# First check the most stricking LV results (order by difference depending on the metric) 

delta_wass = delta_loo = vector(length = 100)

for(i in 1:100){
  tab = realised_list[[i]]$GLV_abioticGR
  
  delta_wass[i] =  colMeans( (tab %>% filter(model == "tSDM" & metric == "wasserstein" &
                   TL !=1 & CV=="CV" & type == "bin") %>% select(value) ) -
                   (tab %>% filter(model == "SDM" & metric == "wasserstein" &
                     TL !=1 & CV=="CV") %>% select(value) ) )
    
  delta_loo[i] =  colMeans( (tab %>% filter(model == "tSDM" & metric == "loo" &
                                                TL !=1)
                                    %>% select(value) ) -
                    (tab %>% filter(model == "SDM" & metric == "loo" &
                                      TL !=1) %>% select(value) ) )
  
  
}
names(delta_loo) = names(delta_wass) = as.character(1:100)

sort(delta_loo)
sort(delta_wass)



# We choose simulation number one (so modify plots to obtain the ones in the presentation)


###############################################################################################################################################################################
#  Plot for all the species and all the plots the pairs...for all methods!

realised_list = do.call(rbind,lapply(realised_list, function(x) do.call(rbind, x)))

fund_list = do.call(rbind, lapply(fund_list, function(x) do.call(rbind, x)))

BioControl_cor = do.call(rbind, lapply(BioControl_cor , function(x) do.call(rbind, x)))
BioControl_cor$TL = gsub(".*\\.TL", "", BioControl_cor$species)
BioControl_cor$id = paste0(BioControl_cor$species,BioControl_cor$source)

BioControl_wass = do.call(rbind, lapply(BioControl_wass , function(x) do.call(rbind, x)))
BioControl_wass$TL = gsub(".*\\.TL", "", BioControl_wass$species)
BioControl_wass$id = paste0(BioControl_wass$species,BioControl_wass$source)

#########################################################################################
######### Realised niche

metrics = c("wasserstein","loo","AUC", "calibration")

# for all methods
for(j in methods){
  
  BioControl_cor_temp = BioControl_cor %>% filter(methods == j, TL != 1)
  
  BioControl_wass_temp = BioControl_wass %>% filter(methods ==j, TL !=1)
  
  p = tab_temp = tab_temp_spread = BC_wass_temp = BC_cor_temp = BC_SDM_temp =  BC_SDM = BC_cor = list()
  anova_TL_temp = t_temp = vector()
  
  # put different metrics in the same plot
  for(k in 1:length(metrics)){
  
  # select the given metric and given method
  if(metrics[k] != "loo"){
  tab_temp[[k]] = realised_list %>% filter(method == j,
                                    type != "prob",
                                    CV== "CV",
                                    metric == metrics[k],
                                    TL != 1) %>% select(-type)
  }else{
    tab_temp[[k]] = realised_list %>% filter(method == j,
                                             metric == metrics[k],
                                             TL != 1) %>% select(-type)
    
  }
  
  # test wheter tSDM improvement is significant
  t_temp[k] = t.test(tab_temp[[k]][which(tab_temp[[k]]$model =="tSDM"),"value"],
                tab_temp[[k]][which(tab_temp[[k]]$model =="SDM"),"value"],
                alternative = ifelse(metrics[k] == "wasserstein","less","greater"),
                paired = TRUE)$p.value
  
  # create species-repetition id
  tab_temp[[k]]$id  = paste0(tab_temp[[k]]$species,tab_temp[[k]]$source)
  
  # compute differences between tSDM and SDM performances
  diff = data.frame(diff=NULL,id=NULL,TL=NULL)
  for(i in unique(tab_temp[[k]]$id)){
    diff = rbind(diff, 
               data.frame( diff = tab_temp[[k]][which(tab_temp[[k]]$id ==i &
                                              tab_temp[[k]]$model == 'tSDM'), "value"] - 
                        tab_temp[[k]][which(tab_temp[[k]]$id ==i &
                                         tab_temp[[k]]$model == 'SDM'),"value"] ,
                      id = i,
                      TL = tab_temp[[k]][which(tab_temp[[k]]$id == i &
                                            tab_temp[[k]]$model == 'tSDM'), "TL"]
               ))
  }

  #Test wheter the trophic level is significant: are predictions better improved for TL 2?
  #anova_TL_temp[k] = summary(lm(diff ~ as.factor(TL), diff))$coefficients[2,4]
  anova_TL_temp[k] = t.test(diff[which(diff$TL == 2), "diff"],
                            diff[which(diff$TL == 3), "diff"], 
                            alternative = ifelse(metrics[k] == "wasserstein",
                                                 "less", "greater")
                            )$p.value
  
  # Test if the improvement depends on the 'biotic control' previously computed
  
  # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences. So we expect negative coefficients
  BC_SDM_temp[[k]] = merge(diff,tab_temp[[k]][which(tab_temp[[k]]$model == 'SDM'),], by = "id")
  
  #if(k==1) BC_SDM_temp[[k]]$diff = - BC_SDM_temp[[k]]$diff
  
  BC_SDM[[k]] = summary(lm(scale(diff)~ scale(value),
                         BC_SDM_temp[[k]]))$coefficients[2,c(1,4)]
  
  # using correlation

  BC_cor_temp[[k]] = merge(BioControl_cor_temp, diff, by = "id")
  
  # If metric = wasserstein change the sense to be consistent: positive coeff means increase in performance with stronger biotic control
  if(k==1) BC_cor_temp[[k]]$diff = - BC_cor_temp[[k]]$diff
  
  # same here: higher correlation means less biotic control so we change sign
  BC_cor[[k]] = summary(lm(scale(diff)~ I(scale(-value)),
                         BC_cor_temp[[k]]))$coefficients[2,c(1,4)]
  
  #using wasserstein distance
  #BC_wass_temp[[k]] = merge(BioControl_wass_temp, diff, by = "id")
  
  # If metric = wasserstein change the sense to be consistent: positive coeff means increase in performance with stronger biotic control
  #if(k==1) BC_wass_temp[[k]]$diff = - BC_wass_temp[[k]]$diff
  
  # We use minus as 
  #summary(lm(diff~ value, BC_wass_temp[[k]]))$coefficients[2,c(1,4)]
  
  
  # Spread the table to do plot
  tab_temp_spread[[k]] = tab_temp[[k]] %>% spread(key = model, value = value)
  
  # Plot species-repetition pairs in the tSDM-SDM performance space
  p[[2*k-1]] = ggplot( data = tab_temp_spread[[k]] ) + 
    geom_point(aes(x = SDM, y = tSDM,
                   col = as.factor(source)),size = 0.8) +
    scale_colour_discrete(guide = "none") + geom_abline(slope=1, intercept=0) +
    annotate("text", x= min(tab_temp_spread[[k]]$SDM) + 2 *sd(tab_temp_spread[[k]]$SDM),
             y = max(tab_temp_spread[[k]]$tSDM) - 0.01,
             label = paste0("TL 2 > 3 pvalue ",
                            ifelse(anova_TL_temp[k]< 2.2*10^(-16), paste0("< ", 2.2*10^(-16)),
                            paste0("= ",round(anova_TL_temp[k], 
                                  digits = -floor(log10(anova_TL_temp[k]))))))) +
    annotate("text", x= max(tab_temp_spread[[k]]$SDM) - 2.2 *sd(tab_temp_spread[[k]]$SDM),
             y = min(tab_temp_spread[[k]]$tSDM) + 0.8*sd(tab_temp_spread[[k]]$tSDM),
             label = paste0("Bio control test: \n SDM: coef =",
                            round(BC_SDM[[k]][1],2)," pval ",
                            ifelse(BC_SDM[[k]][2]< 2.2*10^(-16),
                                   paste0("< ", 2.2*10^(-16)),
                                   paste0("= ",
                                          round(BC_SDM[[k]][2], 
                                                     digits = -floor(log10(BC_SDM[[k]][2]))))
                                   )#,
                            # " \n cor: coef = ", 
                            # round(BC_cor[[k]][1],2)," pval ",
                            # ifelse(BC_cor[[k]][2]< 2.2*10^(-16),
                            #        paste0("< ", 2.2*10^(-16)),
                            #        paste0("= ",
                            #               round(BC_cor[[k]][2], 
                            #                     digits = -floor(log10(BC_cor[[k]][2]))))
                            # )
                            
                            )) +
    theme_classic()
  
  #Boxplots summarising all repetitions
  p[[2*k]] = ggplot(data = tab_temp[[k]]) + geom_boxplot(aes(x = model,
                                                             y = value,
                                                             col = model)) +
            annotate("text", x= 1.5, y = max(tab_temp[[k]]$value)-0.01, 
                      label = paste0("paired t-test pvalue ",
                                      ifelse(t_temp[k]< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                             paste0("= ",round(t_temp[k],
                                                   digits = -floor(log10(t_temp[k]))))))) +
              theme_classic() + theme(legend.position="top") 
  
  if(k!=1) p[[2*k]] =  p[[2*k]] + scale_colour_discrete(guide = "none") 


  }
  
  p = grid.arrange( arrangeGrob(p[[1]], p[[2]], top=metrics[1], ncol =2),
                               arrangeGrob(p[[3]], p[[4]], top=metrics[2], ncol =2),
                               arrangeGrob(p[[5]], p[[6]], top=metrics[3], ncol =2),
                               arrangeGrob(p[[7]], p[[8]], top=metrics[4], ncol =2),
                               ncol=1)
  
  ggsave(p, file = paste0(fig, "/realised_all", j,".pdf"), width = 8, height = 13)
  
}

#########################################################################################
######### Fundamental niche

methods = unique(fund_list$method)

metrics = c("wasserstein", "calibration")

for(j in methods){
  
  p = tab_temp = tab_temp_spread = BC_cor_temp = BC_SDM_temp = BC_cor = BC_SDM = list()
  anova_TL_temp = t_temp = vector()
  
  for(k in 1:length(metrics)){
    
    
    tab_temp[[k]] = fund_list %>% filter(method == j,
                                               CV== "CV",
                                               metric == metrics[k],
                                               TL != 1) 
    
    t_temp[k] = t.test(tab_temp[[k]][which(tab_temp[[k]]$model =="tSDM"),"value"],
                       tab_temp[[k]][which(tab_temp[[k]]$model =="SDM"),"value"],
                       alternative = ifelse(metrics[k] == "wasserstein","less","greater"),
                       paired = TRUE)$p.value
    
    tab_temp[[k]]$id  = paste0(tab_temp[[k]]$species,tab_temp[[k]]$source)
    
    diff = data.frame(diff=NULL,id=NULL,TL=NULL)
    for(i in unique(tab_temp[[k]]$id)){
      diff = rbind(diff, 
                   data.frame( diff = tab_temp[[k]][which(tab_temp[[k]]$id ==i &
                                                            tab_temp[[k]]$model == 'tSDM'), "value"] - 
                                 tab_temp[[k]][which(tab_temp[[k]]$id ==i &
                                                       tab_temp[[k]]$model == 'SDM'),"value"] ,
                               id = i,
                               TL = tab_temp[[k]][which(tab_temp[[k]]$id == i &
                                                          tab_temp[[k]]$model == 'tSDM'), "TL"]
                   ))
    }
    
    #anova_TL_temp[k] = summary(lm(diff ~ as.factor(TL), diff))$coefficients[2,4]
    
    anova_TL_temp[k] = t.test(diff[which(diff$TL == 2), "diff"],
                              diff[which(diff$TL == 3), "diff"], 
                              alternative = ifelse(metrics[k] == "wasserstein","less", "greater")
                              )$p.value
    
    # Test if the improvement depends on the 'biotic control' previously computed
    
    # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences. So we expect negative coefficients
    BC_SDM_temp[[k]] = merge(diff,tab_temp[[k]][which(tab_temp[[k]]$model == 'SDM'),], by = "id")
    
    #if(k==1) BC_SDM_temp[[k]]$diff = - BC_SDM_temp[[k]]$diff
    
    BC_SDM[[k]] = summary(lm(scale(diff)~ scale(value),
                             BC_SDM_temp[[k]]))$coefficients[2,c(1,4)]
    
    # using correlation
    
    BC_cor_temp[[k]] = merge(BioControl_cor_temp, diff, by = "id")
    
    # If metric = wasserstein change the sense to be consistent: positive coeff means increase in performance with stronger biotic control
    if(k==1) BC_cor_temp[[k]]$diff = - BC_cor_temp[[k]]$diff
    
    # same here: higher correlation means less biotic control so we change sign
    BC_cor[[k]] = summary(lm(scale(diff)~ I(scale(-value)),
                             BC_cor_temp[[k]]))$coefficients[2,c(1,4)]
    
      
    tab_temp_spread[[k]] = tab_temp[[k]] %>% spread(key = model, value = value)
    
    p[[2*k-1]] = ggplot( data = tab_temp_spread[[k]] ) + 
      geom_point(aes(x = SDM, y = tSDM,
                     col = as.factor(source)),size = 0.8) +
      scale_colour_discrete(guide = "none") + geom_abline(slope=1, intercept=0) +
      annotate("text", x= min(tab_temp_spread[[k]]$SDM) + 2 *sd(tab_temp_spread[[k]]$SDM),
               y = max(tab_temp_spread[[k]]$tSDM) - 0.01,
               label = paste0("TL 2 > 3 pvalue = ",
                              ifelse(anova_TL_temp< 2.2*10^(-16), paste0(" < ",2.2*10^(-16)),
                                     paste0(" = ",round(anova_TL_temp[k], 
                                           digits = -floor(log10(anova_TL_temp[k]))))))) +
      annotate("text", x= max(tab_temp_spread[[k]]$SDM) - 2.2 *sd(tab_temp_spread[[k]]$SDM),
               y = min(tab_temp_spread[[k]]$tSDM) + 0.4*sd(tab_temp_spread[[k]]$tSDM),
               label = paste0("Bio control test: \n SDM: coef =",
                              round(BC_SDM[[k]][1],2)," pval ",
                              ifelse(BC_SDM[[k]][2]< 2.2*10^(-16),
                                     paste0("< ", 2.2*10^(-16)),
                                     paste0("= ",
                                            round(BC_SDM[[k]][2], 
                                                  digits = -floor(log10(BC_SDM[[k]][2]))))
                              )#,
                              # " \n cor: coef = ", 
                              # round(BC_cor[[k]][1],2)," pval ",
                              # ifelse(BC_cor[[k]][2]< 2.2*10^(-16),
                              #        paste0("< ", 2.2*10^(-16)),
                              #        paste0("= ",
                              #               round(BC_cor[[k]][2], 
                              #                     digits = -floor(log10(BC_cor[[k]][2]))))
                              # )
                              
               )) +
      theme_classic()
    
    
    
    
    p[[2*k]] = ggplot(data = tab_temp[[k]]) + geom_boxplot(aes(x = model,
                                                               y = value,
                                                               col = model)) +
      annotate("text", x= 1.5, y = max(tab_temp[[k]]$value)-0.01, 
               label = paste0("paired t-test pvalue ",
                              ifelse(t_temp[k]< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                     paste0("= ",round(t_temp[k],
                                           digits = -floor(log10(t_temp[k]))))))) +
      theme_classic() + theme(legend.position="top") 
    
    if(k!=1) p[[2*k]] =  p[[2*k]] + scale_colour_discrete(guide = "none") 
    
    
  }
  
  p = grid.arrange( arrangeGrob(p[[1]], p[[2]], top=metrics[1], ncol =2),
                    arrangeGrob(p[[3]], p[[4]], top=metrics[2], ncol =2),
                    ncol=1)
  
  ggsave(p, file = paste0(fig, "/fund_all", j,".pdf"), width = 8, height = 10)
  
}


