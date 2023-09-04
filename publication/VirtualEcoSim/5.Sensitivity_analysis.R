########################################################################################################################
########################################################################################################################
# Script to analyse the evaluation metrics across simulation settings and repetitions from Giovanni Poggiato
# This is the last step of the whole pipeline!


# This scripts merges all the evaluation metrics across all the simulation settings and repetitions

# !!! Only the working directory should be specified to run the script

# By running this script, it creates the final plots summ_bin_real_all.pdf, summ_fund_all.pdf and widthCI.pdf (corresponding to figures S2, S3 and S4) 
# And the tables all_summ_bin_real, regress_bin_real, all_summ_fund, regress_fund corresponding to tables S1, S2, S3, S4.

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


# Load summaries for each parameter setting
summary_bin_realised = readRDS(file = paste0(data_raw, "summary_bin_realised.rds"))
summary_fund = readRDS(file = paste0(data_raw, "summary_fund.rds"))


######### Overall improvements across different settings

### Realised niche predictions

metrics = c("wasserstein","loo","AUC", "calibration")

all_summ_bin_real = data.frame()

t_temp = vector()

p = list()

for(k in 1:length(metrics)){
  
  if(metrics[k] != "loo"){
    tab_temp = bin_realised_list %>% filter(CV == "CV",
                                            metric == metrics[k],
                                            TL != 1)
  }else{
    tab_temp = bin_realised_list %>% filter(metric == metrics[k],
                                            TL != 1)
  }
  
  t_temp[k] = t.test(tab_temp[,"tSDM"],
                     tab_temp[,"SDM"],
                     alternative = ifelse(metrics[k] == "wasserstein","less","greater"),
                     paired = TRUE)$p.value
  
  mean_value = mean(tab_temp[,"tSDM"], na.rm = TRUE)
  
  ## Test wheter the trophic level is significant: are predictions better improved for higher TL ?
  # If so, we expect positive coefficients
  
  lm_TL = lm(diff ~ TL, data = tab_temp)
  
  TL_coef = ifelse(metrics[k] != "wasserstein", coef(lm_TL)[2], -coef(lm_TL)[2])
  
  TL_p_val = summary(lm_TL)$coefficients[2,4]
  
  
  ## Test if the improvement depends on the 'biotic control' previously computed
  
  # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences.
  # So we expect negative coefficients, also for wasserstein (as it's a distance)
  BioControl = summary(lm(scale(diff)~ scale(SDM),
                          tab_temp))$coefficients[2,c(1,4)]
  BC_coef = BioControl[1]
  BC_pval = BioControl[2]
  
  tab_temp_gather = gather(tab_temp, "model", "value", c(SDM,tSDM))
  
  #Boxplots summarising all repetitions
  p[[k]] = ggplot(data = tab_temp_gather) + geom_boxplot(aes(x = model,
                                                             y = value,
                                                             col = model), outlier.shape = NA) +
    ggtitle(metrics[k]) +
    annotate("text", x= 1.5, y = max(tab_temp_gather$value)-0.01, 
             label = paste0("paired t-test pvalue ",
                            ifelse(t_temp[k]< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                   paste0("= ",round(t_temp[k],
                                                     digits = -floor(log10(t_temp[k]))))))) +
    theme_classic() + theme(legend.position="top") 
  
  
  # Mean improvements across simulations (so that results do not depend on the number of species!)
  
  tab_temp = summary_bin_realised  %>% filter(metric == metrics[k])
  
  meanRelDiff = mean(tab_temp$rel_diff[which(tab_temp$rel_diff != Inf)], na.rm = TRUE)
  meanImprSpp = mean(tab_temp$impr_spp, na.rm = TRUE)
  meanMaxImpr = ifelse(metrics[[k]] != "wasserstein",
                       max(tab_temp$max_impr),
                       min(tab_temp$max_impr))
  meanNegSpp20 = mean(tab_temp$neg_spp20, na.rm = TRUE)
  
  
  all_summ_bin_real = rbind(all_summ_bin_real,
                            data.frame(metric = metrics[k],
                                       mean_value = mean_value,
                                       meanRelDiff = meanRelDiff,
                                       t_test_pval = t_temp[k],
                                       meanImprSpp = meanImprSpp,
                                       meanMaxImpr = meanMaxImpr,
                                       meanNegSpp20 = meanNegSpp20,
                                       TL_coef = TL_coef,
                                       TL_p_val = TL_p_val,
                                       BC_coef = BC_coef,
                                       BC_pval = BC_pval
                            ))
  
  
}


p = grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],
                 ncol=2)

ggsave(p, file = paste0(fig, "/summ_bin_real_all.pdf"), width = 8, height = 13)

all_summ_bin_real

### Fundamental niche

metrics = c("wasserstein", "calibration")

all_summ_fund = data.frame()

t_temp = vector()

p = list()

for(k in 1:length(metrics)){
  
  tab_temp = fund_list %>% filter(metric == metrics[k],
                                  CV == "CV",
                                  TL != 1)
  
  
  t_temp[k] = t.test(tab_temp[,"tSDM"],
                     tab_temp[,"SDM"],
                     alternative = ifelse(metrics[k] != "calibration","less","greater"),
                     paired = TRUE)$p.value
  
  mean_value = mean(tab_temp[,"tSDM"], na.rm = TRUE)
  
  ## Test wheter the trophic level is significant: are predictions better improved for higher TL ?
  # If so, we expect positive coefficients
  
  lm_TL = lm(diff ~ TL, data = tab_temp)
  
  TL_coef = ifelse(metrics[k] == "calibration", coef(lm_TL)[2], -coef(lm_TL)[2])
  
  TL_p_val = summary(lm_TL)$coefficients[2,4]
  
  
  ## Test if the improvement depends on the 'biotic control' previously computed
  
  # diff ~ SDM performances: we expect that the less SDM are good, the biggest the differences.
  # So we expect negative coefficients, also for wasserstein (as it's a distance)
  BioControl = summary(lm(scale(diff)~ scale(SDM),
                          tab_temp))$coefficients[2,c(1,4)]
  BC_coef = BioControl[1]
  BC_pval = BioControl[2]
  
  
  tab_temp_gather = gather(tab_temp, "model", "value", c(SDM,tSDM))
  
  #Boxplots summarising all repetitions
  p[[k]] = ggplot(data = tab_temp_gather) + geom_boxplot(aes(x = model,
                                                             y = value,
                                                             col = model), outlier.shape = NA) +
    ggtitle(metrics[k]) +
    annotate("text", x= 1.5, y = max(tab_temp_gather$value)-0.01, 
             label = paste0("paired t-test pvalue ",
                            ifelse(t_temp[k]< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                   paste0("= ",round(t_temp[k],
                                                     digits = -floor(log10(t_temp[k]))))))) +
    theme_classic() + theme(legend.position="top") 
  
  
  
  # Mean improvements across simulations (so that results do not depend on the number of species!)
  
  tab_temp = summary_fund %>% filter(metric == metrics[k])
  
  meanRelDiff = mean(tab_temp$rel_diff[which(tab_temp$rel_diff != Inf)], na.rm = TRUE)
  meanImprSpp = mean(tab_temp$impr_spp, na.rm = TRUE)
  meanMaxImpr = ifelse(metrics[[k]] == "calibration",
                       max(tab_temp$max_impr[which(tab_temp$max_impr != Inf)]),
                       min(tab_temp$max_impr))
  meanNegSpp20 = mean(tab_temp$neg_spp20, na.rm = TRUE)
  
  
  all_summ_fund = rbind(all_summ_fund,
                            data.frame(metric = metrics[k],
                                       mean_value = mean_value,
                                       meanRelDiff = meanRelDiff,
                                       t_test_pval = t_temp[k],
                                       meanImprSpp = meanImprSpp,
                                       meanMaxImpr = meanMaxImpr,
                                       meanNegSpp20 = meanNegSpp20,
                                       TL_coef = TL_coef,
                                       TL_p_val = TL_p_val,
                                       BC_coef = BC_coef,
                                       BC_pval = BC_pval
                            ))
  
  
}


p = grid.arrange(p[[1]],p[[2]],
                 ncol=2)

ggsave(p, file = paste0(fig, "/summ_fund_all.pdf"), width = 8, height = 13)




####### Sensitivity analysis

## Realised niche binary predictions

metrics = c("wasserstein","AUC")

regress_bin_real = data.frame()

for(k in 1:length(metrics)){
  
  tab_temp = summary_bin_realised %>% filter(metric == metrics[k]) %>% select(-source,-metric) 
  tab_temp = as.data.frame(apply(tab_temp, 2, scale))
  
  if(metrics[k] == "wasserstein") tab_temp$rel_diff = -tab_temp$rel_diff
  
  model = summary(lm(rel_diff ~ S + L + p + nicheBreadthRatio + nEnv + nRep,
                     data = tab_temp))
  
  regress_bin_real = rbind(regress_bin_real,
                           cbind(metric = metrics[k],
                                 variable = rownames(model$coefficients),
                                 model$coefficients[,c(1,4)]))
  
}


## Fundamental niche predictions

regress_fund = data.frame()


tab_temp = summary_fund %>% filter(metric =="wasserstein") %>% select(-source,-metric) 
tab_temp = as.data.frame(apply(tab_temp, 2, scale))
  
tab_temp$rel_diff = -tab_temp$rel_diff
  
model = summary(lm(rel_diff ~ S + L + p + nicheBreadthRatio + nEnv + nRep,
                     data = tab_temp))
  
regress_fund = model$coefficients[,c(1,4)]
  


# Effect of CI
CI_width = readRDS(file = "CI_width.rds")

CI_width = filter(CI_width, type == "bin")

mod = summary(lm(scale(value) ~ scale(TL) + scale(PreyNb), CI_width))

p_1 = ggplot(data=CI_width) + geom_boxplot(aes(x=TL, y = value, group= TL),
                                           color = "#F8766D", outlier.shape = NA) + theme_classic() +
  xlab("Trophic level") + ylab("Width of credible intervals")

p_2 = ggplot(data=CI_width) + geom_boxplot(aes(x=PreyNb, y = value, group= PreyNb),
                                           color = "#F8766D", outlier.shape = NA) + theme_classic() +
  xlab("Number of preys") + ylab("Width of credible intervals")

p = grid.arrange(p_1, p_2, ncol = 2)

ggsave(p, file = paste0(fig, "/widthCI.pdf"), width = 15, height = 8)
