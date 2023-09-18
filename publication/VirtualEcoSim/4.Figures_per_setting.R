########################################################################################################################
########################################################################################################################
# Script for obtaining plots for each simulation setting (across the different repetitions) from Giovanni Poggiato

# This scripts merges all the evaluation metrics across the 100 repetitions for each simulation settings

# !!! Only the working directory should be specified to run the script

# By running this script, it will create subfolders in the "/SimulationOutputs/Fig_all", each containing the summary figures for a given parameter setting
# and it creates the files summary_bin_realised, summary_fund, fund_list and CI_width in "/SimulationOutputs", containing sumary metrics for every parameter settings (across the 100 repetitions)

############################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

rm(list = ls())

# ! Set working directory
dir = ""

data_raw = paste0(dir,"/SimulationOutputs/")
fig = paste0(dir,"/SimulationOutputs/Fig_all")


# Load overall lists
bin_realised_list = readRDS(file=paste0(data_raw,"bin_realised_list.rds"))
fund_list =readRDS(file=paste0(data_raw,"fund_list.rds"))


###############################################################################################################################################################################
#  Plot for all the species and all the plots the pairs...for all methods!

all_sett = unique(sub("_job.*","",bin_realised_list$source))

######### Realised niche

metrics = c("wasserstein","loo","AUC", "calibration")

summary_bin_realised = data.frame()

# for all methods
for(j in all_sett){
  
  if(!dir.exists(paste0(fig,
                        "/",
                        sub("Simulations_","",j)))) dir.create(paste0(fig,
                                                                      "/",
                                                                      sub("Simulations_","",j)))
  
  p = tab_temp = tab_temp_gather= BC_SDM = list()
  anova_TL_p_val =  anova_TL_coef = t_temp = vector()
  
  # put different metrics in the same plot
  for(k in 1:length(metrics)){
    
    # select the given metric and given method
    if(metrics[k] != "loo"){
      tab_temp[[k]] = bin_realised_list %>% filter(source %in% paste0(j, "_job",1:100),
                                                   CV == "CV",
                                                   metric == metrics[k],
                                                   TL != 1)
      
    }else{
      tab_temp[[k]] = bin_realised_list %>% filter(source %in% paste0(j, "_job",1:100),
                                                   metric == metrics[k],
                                                   TL != 1)
      
    }
    
    # test wheter tSDM improvement is significant
    t_temp[k] = t.test(tab_temp[[k]][,"tSDM"],
                       tab_temp[[k]][,"SDM"],
                       alternative = ifelse(metrics[k] == "wasserstein","less","greater"),
                       paired = TRUE)$p.value
    
    rel_diff = tab_temp[[k]]$diff/abs(tab_temp[[k]][,"SDM"])
    
    impr_spp = ifelse(metrics[k] != "wasserstein",
                      length(which(tab_temp[[k]]$diff >= 0))/length(tab_temp[[k]]$diff[which(!is.na(tab_temp[[k]]$diff))]),
                      length(which(tab_temp[[k]]$diff <= 0))/length(tab_temp[[k]]$diff[which(!is.na(tab_temp[[k]]$diff))]))
    
    max_impr = ifelse(metrics[k] != "wasserstein", max(rel_diff[which(rel_diff != Inf)]),min(rel_diff))
    
    neg_spp20 = ifelse(metrics[k] != "wasserstein",
                       length(which(rel_diff < -0.2))/length(rel_diff[which(!is.na(rel_diff))]),
                       length(which(rel_diff > 0.2))/length(rel_diff[which(!is.na(rel_diff))]))
    
    
    which_neg = ifelse(metrics[k] != "wasserstein", which(tab_temp[[k]]$diff < 0), which(tab_temp[[k]]$diff > 0))
    
    mean_neg = median(tab_temp[[k]]$diff[which_neg]/abs(tab_temp[[k]][which_neg,"SDM"]), na.rm = TRUE)
    
    ## Test wheter the trophic level is significant: are predictions better improved for higher TL ?
    # If so, we expect positive coefficients
    
    lm_TL = lm(diff ~ TL, data = tab_temp[[k]])
    
    anova_TL_coef[k] = ifelse(metrics[k] != "wasserstein", coef(lm_TL)[2], -coef(lm_TL)[2])
    
    anova_TL_p_val[k] = summary(lm_TL)$coefficients[2,4]
    
    
    ## Test if the improvement depends on the 'biotic control' previously computed
    
    # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences.
    # So we expect negative coefficients, also for wasserstein (as it's a distance)
    BC_SDM[[k]] = summary(lm(scale(diff)~ scale(SDM),
                             tab_temp[[k]]))$coefficients[2,c(1,4)]
    
    # Summary table
    summary_bin_realised = rbind(summary_bin_realised,
                                 data.frame(source = sub("Simulations_","",j),
                                            metric = metrics[k],
                                            S = unique(tab_temp[[k]]$S),
                                            L = unique(tab_temp[[k]]$L),
                                            p = unique(tab_temp[[k]]$p),
                                            nbLink = mean(unique(tab_temp[[k]]$nbLink)),
                                            meanConnectivity = mean(unique(tab_temp[[k]]$meanConnectivity)),
                                            meanStrength = mean(unique(tab_temp[[k]]$meanStrength)),
                                            meanStrengthS = mean(unique(tab_temp[[k]]$meanStrenghS)),
                                            nicheBreadthRatio = mean(unique(tab_temp[[k]]$nicheBreadthRatio)),
                                            nEnv = unique(tab_temp[[k]]$nEnv),
                                            nRep = unique(tab_temp[[k]]$nRep),
                                            mean_value = mean(tab_temp[[k]]$tSDM),
                                            rel_diff = mean(rel_diff[which(rel_diff != Inf)], na.rm = TRUE),
                                            impr_spp = impr_spp,
                                            max_impr = max_impr,
                                            neg_spp20 = neg_spp20,
                                            mean_neg = mean_neg,
                                            BC_coef = as.numeric(BC_SDM[[k]][1]),
                                            BC_pval = as.numeric(BC_SDM[[k]][2]),
                                            TL_coef = anova_TL_coef[k],
                                            TL_pval = anova_TL_p_val[k]
                                 ))
    
    # Plot species-repetition pairs in the tSDM-SDM performance space
    p[[2*k-1]] = ggplot( data = tab_temp[[k]] ) + 
      geom_point(aes(x = SDM, y = tSDM,
                     col = as.factor(source)),size = 0.8) +
      scale_colour_discrete(guide = "none") + geom_abline(slope=1, intercept=0) +
      annotate("text", x= min(tab_temp[[k]]$SDM) + 2 *sd(tab_temp[[k]]$SDM),
               y = max(tab_temp[[k]]$tSDM) - 0.01,
               label = paste0("TL pvalue ",
                              ifelse(anova_TL_p_val[k]< 2.2*10^(-16), paste0("< ", 2.2*10^(-16)),
                                     paste0("= ",round(anova_TL_p_val[k], 
                                                       digits = -floor(log10(anova_TL_p_val[k]))))),
                              "\n TL coef ", round(anova_TL_coef[k], digits = 2))) +
      annotate("text", x= max(tab_temp[[k]]$SDM) - 2.2 *sd(tab_temp[[k]]$SDM),
               y = min(tab_temp[[k]]$tSDM) + 0.8*sd(tab_temp[[k]]$tSDM),
               label = paste0("Bio control test: \n coef =",
                              round(BC_SDM[[k]][1],2)," \n pval ",
                              ifelse(BC_SDM[[k]][2]< 2.2*10^(-16),
                                     paste0("< ", 2.2*10^(-16)),
                                     paste0("= ",
                                            round(BC_SDM[[k]][2], 
                                                  digits = -floor(log10(BC_SDM[[k]][2]))))
                              ))) +
      theme_classic()
    
    tab_temp_gather[[k]] = gather(tab_temp[[k]], "model", "value", c(SDM,tSDM))
    
    #Boxplots summarising all repetitions
    p[[2*k]] = ggplot(data = tab_temp_gather[[k]]) + geom_boxplot(aes(x = model,
                                                                      y = value,
                                                                      col = model)) +
      annotate("text", x= 1.5, y = max(tab_temp_gather[[k]]$value)-0.01, 
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
  
  ggsave(p, file = paste0(fig, "/",sub("Simulations_","",j),"/bin_realised_all.pdf"), width = 8, height = 13)
  
}


#########################################################################################
######### Potential niche

metrics = c("wasserstein", "calibration", "wass_growth_rate", "dist_niche_opt")


summary_fund = data.frame()

# for all methods
for(j in all_sett){
  
  if(!dir.exists(paste0(fig,
                        "/",
                        sub("Simulations_","",j)))) dir.create(paste0(fig,
                                                                      "/",
                                                                      sub("Simulations_","",j)))
  
  p = tab_temp = tab_temp_gather= BC_SDM = list()
  anova_TL_p_val =  anova_TL_coef = t_temp = vector()
  
  # put different metrics in the same plot
  for(k in 1:length(metrics)){
    
    # select the given metric
    
    tab_temp[[k]] = fund_list %>% filter(source %in% paste0(j, "_job",1:100),
                                         metric == metrics[k],
                                         CV == "CV",
                                         TL != 1)
    
    
    # test wheter tSDM improvement is significant
    t_temp[k] = t.test(tab_temp[[k]][,"tSDM"],
                       tab_temp[[k]][,"SDM"],
                       alternative = ifelse(metrics[k] != "calibration","less","greater"),
                       paired = TRUE)$p.value
    
    rel_diff = tab_temp[[k]]$diff/abs(tab_temp[[k]][,"SDM"])
    
    impr_spp = ifelse(metrics[k] == "calibration",
                      length(which(tab_temp[[k]]$diff >= 0))/length(tab_temp[[k]]$diff[which(!is.na(tab_temp[[k]]$diff))]),
                      length(which(tab_temp[[k]]$diff <= 0))/length(tab_temp[[k]]$diff[which(!is.na(tab_temp[[k]]$diff))]))
    
    max_impr = ifelse(metrics[k] == "calibration", max(rel_diff[which(rel_diff != Inf)]),min(rel_diff))
    
    neg_spp20 = ifelse(metrics[k] == "calibration",
                       length(which(rel_diff < -0.2))/length(rel_diff[which(!is.na(rel_diff))]),
                       length(which(rel_diff > 0.2))/length(rel_diff[which(!is.na(rel_diff))]))
    
    
    which_neg = ifelse(metrics[k] == "calibration", which(tab_temp[[k]]$diff < 0), which(tab_temp[[k]]$diff > 0))
    
    mean_neg = median(tab_temp[[k]]$diff[which_neg]/abs(tab_temp[[k]][which_neg,"SDM"]), na.rm = TRUE)
    
    ## Test wheter the trophic level is significant: are predictions better improved for higher TL ?
    # If so, we expect positive coefficients
    
    lm_TL = lm(diff ~ TL, data = tab_temp[[k]])
    
    anova_TL_coef[k] = ifelse(metrics[k] == "calibration", coef(lm_TL)[2], -coef(lm_TL)[2])
    
    anova_TL_p_val[k] = summary(lm_TL)$coefficients[2,4]
    
    
    ## Test if the improvement depends on the 'biotic control' previously computed
    
    # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences.
    # So we expect negative coefficients, also for wasserstein (as it's a distance)
    BC_SDM[[k]] = summary(lm(scale(diff)~ scale(SDM),
                             tab_temp[[k]]))$coefficients[2,c(1,4)]
    
    # Summary table
    summary_fund = rbind(summary_fund,
                         data.frame(source = sub("Simulations_","",j),
                                    metric = metrics[k],
                                    S = unique(tab_temp[[k]]$S),
                                    L = unique(tab_temp[[k]]$L),
                                    p = unique(tab_temp[[k]]$p),
                                    nbLink = mean(unique(tab_temp[[k]]$nbLink)),
                                    meanConnectivity = mean(unique(tab_temp[[k]]$meanConnectivity)),
                                    meanStrength = mean(unique(tab_temp[[k]]$meanStrength)),
                                    meanStrengthS = mean(unique(tab_temp[[k]]$meanStrenghS)),
                                    nicheBreadthRatio = mean(unique(tab_temp[[k]]$nicheBreadthRatio)),
                                    nEnv = unique(tab_temp[[k]]$nEnv),
                                    nRep = unique(tab_temp[[k]]$nRep),
                                    mean_value = mean(tab_temp[[k]]$tSDM),
                                    rel_diff = mean(rel_diff[which(rel_diff != Inf)],  na.rm = TRUE),
                                    impr_spp = impr_spp,
                                    max_impr = max_impr,
                                    neg_spp20 = neg_spp20,
                                    mean_neg = mean_neg,
                                    BC_coef = as.numeric(BC_SDM[[k]][1]),
                                    BC_pval = as.numeric(BC_SDM[[k]][2]),
                                    TL_coef = anova_TL_coef[k],
                                    TL_pval = anova_TL_p_val[k]
                         ))
    
    # Plot species-repetition pairs in the tSDM-SDM performance space
    p[[2*k-1]] = ggplot( data = tab_temp[[k]] ) + 
      geom_point(aes(x = SDM, y = tSDM,
                     col = as.factor(source)),size = 0.8) +
      scale_colour_discrete(guide = "none") + geom_abline(slope=1, intercept=0) +
      annotate("text", x= min(tab_temp[[k]]$SDM) + 2 *sd(tab_temp[[k]]$SDM),
               y = max(tab_temp[[k]]$tSDM) - 0.01,
               label = paste0("TL pvalue ",
                              ifelse(anova_TL_p_val[k]< 2.2*10^(-16), paste0("< ", 2.2*10^(-16)),
                                     paste0("= ",round(anova_TL_p_val[k], 
                                                       digits = -floor(log10(anova_TL_p_val[k]))))),
                              "\n TL coef ", round(anova_TL_coef[k], digits = 2))) +
      annotate("text", x= max(tab_temp[[k]]$SDM) - 2.2 *sd(tab_temp[[k]]$SDM),
               y = min(tab_temp[[k]]$tSDM) + 0.8*sd(tab_temp[[k]]$tSDM),
               label = paste0("Bio control test: \n coef =",
                              round(BC_SDM[[k]][1],2)," \n pval ",
                              ifelse(BC_SDM[[k]][2]< 2.2*10^(-16),
                                     paste0("< ", 2.2*10^(-16)),
                                     paste0("= ",
                                            round(BC_SDM[[k]][2], 
                                                  digits = -floor(log10(BC_SDM[[k]][2]))))
                              ))) +
      theme_classic()
    
    tab_temp_gather[[k]] = gather(tab_temp[[k]], "model", "value", c(SDM,tSDM))
    
    #Boxplots summarising all repetitions
    p[[2*k]] = ggplot(data = tab_temp_gather[[k]]) + geom_boxplot(aes(x = model,
                                                                      y = value,
                                                                      col = model)) +
      annotate("text", x= 1.5, y = max(tab_temp_gather[[k]]$value)-0.01, 
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
  
  ggsave(p, file = paste0(fig, "/",sub("Simulations_","",j),"/fund_all.pdf"), width = 8, height = 13)
  
}


saveRDS(summary_bin_realised, file = paste0(data_raw, "summary_bin_realised.rds"))
saveRDS(summary_fund, file = paste0(data_raw, "summary_fund.rds"))

