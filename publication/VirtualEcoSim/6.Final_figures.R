########################################################################################################################
########################################################################################################################
# Script for reproducing figures 4 and 5 from Giovanni Poggiato
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
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cheddar)


rm(list = ls())

dir = ""

data_raw = paste0(dir,"/SimulationOutputs/")
fig = paste0(dir,"/SimulationOutputs/Fig_all")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


###### Figure 4

### Realised niche

bin_realised_list = readRDS(file=paste0(data_raw,"bin_realised_list.rds"))


tab_temp = bin_realised_list %>% filter(metric == "AUC",
                                        CV == "CV",
                                        TL != 1)


t_temp = t.test(tab_temp[,"tSDM"],
                tab_temp[,"SDM"],
                alternative = "greater",
                paired = TRUE)$p.value

tab_temp_gather = gather(tab_temp, "model", "value", c(SDM,tSDM))

#Boxplots summarising all repetitions
p = ggplot(data = tab_temp_gather) + geom_boxplot(aes(x = model,
                                                      y = value,
                                                      col = model), outlier.shape = NA) +
  annotate("text", x= 1.5, y = max(tab_temp_gather$value)-0.01, 
           label = paste0("paired t-test pvalue ",
                          ifelse(t_temp< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                 paste0("= ",round(t_temp,
                                                   digits = -floor(log10(t_temp))))))) +
  theme_classic() + theme(legend.position="top") 

ggsave(p, file = paste0(fig, "/AUC_realised.pdf"), height = 15, width = 10)


### Potential niche

fund_list = readRDS(file=paste0(data_raw,"fund_list.rds"))


tab_temp = fund_list %>% filter(metric == "wasserstein",
                                CV == "CV",
                                TL != 1)

t_temp = t.test(tab_temp[,"tSDM"],
                tab_temp[,"SDM"],
                alternative = "less",
                paired = TRUE)$p.value

tab_temp_gather = gather(tab_temp, "model", "value", c(SDM,tSDM))

#Boxplots summarising all repetitions
p = ggplot(data = tab_temp_gather) + geom_boxplot(aes(x = model,
                                                      y = value,
                                                      col = model), outlier.shape = NA) +
  annotate("text", x= 1.5, y = max(tab_temp_gather$value)-0.01, 
           label = paste0("paired t-test pvalue ",
                          ifelse(t_temp< 2.2*10^(-16), paste0("< ",2.2*10^(-16)),
                                 paste0("= ",round(t_temp,
                                                   digits = -floor(log10(t_temp))))))) +
  theme_classic() + theme(legend.position="top") 

ggsave(p, file = paste0(fig, "/Wass_fund.pdf"), height = 15, width = 10)


###### Figure 5


#### Low normalised niche breadth (here niche breadth ratio)

source_dir = "Simulations_S_20_L3_p0.2_nicheBreadthRatio2_nEnv50_nRep50_job25"

# AUC values
filter(bin_realised_list, source == source_dir, metric == "AUC", CV == "CV", species == "Sp13.TL2")


SIM = get(load(paste0(data_raw, source_dir, "/SIMlist.RData")))

nRep = 50
nEnv = 50
nbMerge = 1
S = 20
intercept = TRUE
nEnv = nrow(SIM$X)/(nRep*nbMerge)
envs = SIM$X[1:nEnv, ifelse(intercept, -1, -(K+1))]


IntMat <- as.matrix(read.csv2(paste0(data_raw, source_dir,"/InteractionMatrix.csv"), row.names = 1))  # matrix of species interactions

PredMat <- sign(IntMat) ; PredMat[lower.tri(PredMat, diag=TRUE)] <- 0                    # matrix of trophic links
spNames = colnames(PredMat)        
strengthBI = 5
# species names
trophL = as.numeric(sub(pattern = "Sp.*TL", replacement = "", x = spNames))         # species trophic levels
#Stroph = table(trophL)   
community <- Community(nodes=data.frame(node=spNames),
                       trophic.links=PredationMatrixToLinks(PredMat),
                       properties=list(title='Random community'))

load(paste0(data_raw, source_dir,"/spNewNames.R"))

nameOrder = cumsum(table(spNewNames)[spNewNames])[paste0("Y", 1:S)]  # table to alternate between name ordering

proba = list()

proba$pred = SIM$p.mean.stan_bin
proba$est.02 = SIM$p.qinf.stan_bin
proba$est.97 = SIM$p.qsup.stan_bin

proba$SDM.pred = SIM$SDM.p.mean.stan
proba$SDM.est.02 = SIM$SDM.p.qinf.stan
proba$SDM.est.97 = SIM$SDM.p.qsup.stan

realizedNiche = SIM$prob

# build niche table
table=data.frame(env=envs,rbind(data.frame(stack(realizedNiche),type=rep("True",nrow(stack(realizedNiche)))),
                                data.frame(stack(as.data.frame(proba$pred)),type=rep("tSDM",nrow(stack(realizedNiche)))),
                                data.frame(stack(as.data.frame(proba$SDM.pred)),type=rep("SDM",nrow(stack(realizedNiche))))))


Preys=data.frame()
for(s in 1:S){
  for (prey in community$trophic.links[community$trophic.links$consumer == spNames[nameOrder[s]],]$resource){
    intStrength = abs(IntMat[prey,spNames[nameOrder[s]]])/(strengthBI*2)  # strength of the biotic interaction
    Preys = rbind(Preys, data.frame(focal=paste0("Y", gsub("Sp|\\..*", "", spNames[nameOrder[s]])),env=envs,
                                    value=realizedNiche[,paste0("Y", gsub("Sp|\\..*", "", prey))]*intStrength,
                                    predprey=paste0("Y", gsub("Sp|\\..*", "", prey))))
    
  }
}

s = 13

table_plot = table[which(table$ind==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),]
tablePrey = Preys[which(Preys$focal==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),]


cols = c("grey75",  gg_color_hue(2)[2], gg_color_hue(2)[1])

p1 = ggplot() +
  geom_line(data=table_plot[table_plot$type=="True",],
            aes(x=env,y=values),lwd=1, col=cols[1])+
  geom_line(data=table_plot[table_plot$type=="tSDM",],
            aes(x=env,y=values),lwd=1, lty=1, col=cols[2])+
  geom_line(data=table_plot[table_plot$type=="SDM",],
            aes(x=env,y=values),lwd=1, lty=2, col=cols[3]) +
  theme_classic()+ggtitle(spNames[s]) +
  xlab("Env") + ylab("Probability of presence")

p1 = p1 + geom_line(data=tablePrey,
                  aes(x=env,y=value, group=predprey),
                  color="blue",lwd =0.5 ,linetype=2,alpha=0.3)






#### High normalised niche breadth (here niche breadth ratio)

source_dir = "Simulations_S_20_L3_p0.2_nicheBreadthRatio10_nEnv50_nRep70_job56"

# AUC values
filter(bin_realised_list, source == source_dir, metric == "AUC", CV == "CV", species == "Sp15.TL2")

SIM = get(load(paste0(data_raw, source_dir, "/SIMlist.RData")))

IntMat <- as.matrix(read.csv2(paste0(data_raw, source_dir,"/InteractionMatrix.csv"),
                              row.names = 1))  # matrix of species interactions
PredMat <- sign(IntMat) ; PredMat[lower.tri(PredMat, diag=TRUE)] <- 0                    # matrix of trophic links
spNames = colnames(PredMat)       

# species names
trophL = as.numeric(sub(pattern = "Sp.*TL", replacement = "", x = spNames))         # species trophic levels
#Stroph = table(trophL)   
community <- Community(nodes=data.frame(node=spNames),
                       trophic.links=PredationMatrixToLinks(PredMat),
                       properties=list(title='Random community'))

load(paste0(data_raw,source_dir, "/spNewNames.R"))

nameOrder = cumsum(table(spNewNames)[spNewNames])[paste0("Y", 1:S)]  # table to alternate between name ordering


proba = list()


proba$pred = SIM$p.mean.stan_bin
proba$est.02 = SIM$p.qinf.stan_bin
proba$est.97 = SIM$p.qsup.stan_bin

proba$SDM.pred = SIM$SDM.p.mean.stan
proba$SDM.est.02 = SIM$SDM.p.qinf.stan
proba$SDM.est.97 = SIM$SDM.p.qsup.stan

realizedNiche = SIM$prob

# build niche table
table=data.frame(env=envs,rbind(data.frame(stack(realizedNiche),type=rep("True",nrow(stack(realizedNiche)))),
                                data.frame(stack(as.data.frame(proba$pred)),type=rep("tSDM",nrow(stack(realizedNiche)))),
                                data.frame(stack(as.data.frame(proba$SDM.pred)),type=rep("SDM",nrow(stack(realizedNiche))))))


Preys=data.frame()
for(s in 1:S){
  for (prey in community$trophic.links[community$trophic.links$consumer == spNames[nameOrder[s]],]$resource){
    intStrength = abs(IntMat[prey,spNames[nameOrder[s]]])/(strengthBI*2)  # strength of the biotic interaction
    Preys = rbind(Preys, data.frame(focal=paste0("Y", gsub("Sp|\\..*", "", spNames[nameOrder[s]])),env=envs,
                                    value=realizedNiche[,paste0("Y", gsub("Sp|\\..*", "", prey))]*intStrength,
                                    predprey=paste0("Y", gsub("Sp|\\..*", "", prey))))
    
  }
}

s = 15

table_plot = table[which(table$ind==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),]
tablePrey = Preys[which(Preys$focal==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),]


cols = c("grey75",  gg_color_hue(2)[2], gg_color_hue(2)[1])

p2 = ggplot() +
  geom_line(data=table_plot[table_plot$type=="True",],
            aes(x=env,y=values),lwd=1, col=cols[1])+
  geom_line(data=table_plot[table_plot$type=="tSDM",],
            aes(x=env,y=values),lwd=1, lty=1, col=cols[2])+
  geom_line(data=table_plot[table_plot$type=="SDM",],
            aes(x=env,y=values),lwd=1, lty=2, col=cols[3]) +
  theme_classic()+ggtitle(spNames[s]) +
  xlab("Env") + ylab("Probability of presence")

p2 = p2 + geom_line(data=tablePrey,
                    aes(x=env,y=value, group=predprey),
                    color="blue",lwd =0.5, linetype=2,alpha=0.3)

p = grid.arrange(p1,p2, ncol=2)

# save
ggsave(p, file = paste0(fig, "/Fig5_raw.pdf"), width = 13, height = 8)

