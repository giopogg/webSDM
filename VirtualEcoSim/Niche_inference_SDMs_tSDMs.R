####### Jeremy: tSDMs vs. SDMs inference of realized and fundamental niches

rm(list=ls())

library(ggplot2)
library(brms)
library(igraph)
library(cheddar)
library(gridExtra)


root= "/Users/poggiatg/Documents/Phd/Futureweb/Code/SIM_JEREMY"

setwd(root)


if (grepl("A","CA")) print(1)

fitPreds=F
horsh=F
richness=F
simPath = "Simulations_S20L3_nEnv51_nRep50_strengthBI5_asymBI2/"    # directory from which to load the simulations

if(!fitPreds){if(!horsh)figPath=paste0(simPath,"Preys_Gio/") else figPath=paste0(simPath,"Preys_Horshoe/")}else{figPath=simPath}
#simPath = "Simulations_S6L3_nEnv51_nRep50_strengthBI5_asymBI2/"
load(file=paste0(figPath, "SIMlist.Rdata"))


S = as.numeric(gsub(".*S|L.*", "", simPath))                       # number of species
L = as.numeric(gsub(".*L|_nEnv.*", "", simPath))                   # number of trophic levels
simMethods = gsub("_finalStates_abiotic|.csv", "", grep("finalStates", list.files(simPath), value = TRUE))
nEnv = as.numeric(gsub(".*nEnv|_nRep.*", "", simPath))             # number of environments
nRep = as.numeric(gsub(".*nRep|_strength.*", "", simPath))         # number of replicates
strengthBI = as.numeric(gsub(".*strengthBI|_asym.*", "", simPath)) # strength of biotic interactions
mergeReplicates = TRUE                   # wether to merge the simulation replicates to train the models
trainFraction = 20/nRep                  # fraction of the simulations used to train the models
testFraction = 1/nRep                    # fraction of the simulations used to test the models
#testFraction = 1-trainFraction
intercept = TRUE                         # wether to infer intercepts
linear = F                               # wether to include a quadratic term for X
fitPreds = F                             # wether to fit predators
if(!linear) poly=F                       # if yes, whether to take a raw polynomial (poly=T) or an orthogonal one (poly=F)

load(paste0(simPath,"spNewNames.R"))



# Load niche optima used in the simulations
niche_optima <- read.csv2(paste0(simPath, "niche_optima.csv"))[[1]]
names(niche_optima) = paste0("Y", 1:S)
niche_optima = niche_optima[spNewNames]
K = ifelse(is.vector(niche_optima), 1, ncol(niche_optima))    # number of environmental variables



IntMat <- as.matrix(read.csv2(paste0(simPath, "InteractionMatrix.csv"), row.names = 1))  # matrix of species interactions
PredMat <- sign(IntMat) ; PredMat[lower.tri(PredMat, diag=TRUE)] <- 0                    # matrix of trophic links
#PredMat2 <- PredMat
#while(is.dag(graph_from_adjacency_matrix(t(PredMat2))) | !any(PredMat2 + t(PredMat2) == 2)){
#    # add some interactions to test with a non-DAG
#    PredMat2 = PredMat + (1-PredMat)*rbinom(length(PredMat), 1, 0.05)
#    diag(PredMat2) <- 0
#}
#IntMat <- IntMat + (PredMat2-PredMat) - t(PredMat2-PredMat)/2
#PredMat <- PredMat2
spNames = colnames(PredMat)                                                         # species names
trophL = as.numeric(sub(pattern = "Sp.*TL", replacement = "", x = spNames))         # species trophic levels
Stroph = table(trophL)                                                              # number of species by trophic level
community <- Community(nodes=data.frame(node=spNames),
                       trophic.links=PredationMatrixToLinks(PredMat),
                       properties=list(title='Random community'))                   # cheddar community
dimnames(PredMat) <- list(paste0("Y", 1:S),paste0("Y", 1:S))

# Build the graph of trophic interactions from PredMat
G = graph_from_adjacency_matrix(t(PredMat))


####### Comparisons
####Distributions


# Load simulated fundamental niches, corresponding to probabilities of presence when all preys and no predator are present.

if ("glvGR" %in% simMethods & file.exists(paste0(simPath, "glv.fundNicheTh.abioticGR.csv"))){
  SIMlist$GLV_abioticGR$fundNiche <- read.csv2(paste0(simPath, "glv.fundNicheTh.abioticGR.csv"), row.names=1)
  colnames(SIMlist$GLV_abioticGR$fundNiche) <- paste0("Y",1:S)
  SIMlist$GLV_abioticGR$fundNiche <- SIMlist$GLV_abioticGR$fundNiche[,spNewNames]
}

if ("glvKbasal" %in% simMethods & file.exists(paste0(simPath, "glv.fundNicheTh.abioticKbasal.csv"))){
  SIMlist$GLV_abioticKbasal$fundNiche <- read.csv2(paste0(simPath, "glv.fundNicheTh.abioticKbasal.csv"), row.names=1)
  colnames(SIMlist$GLV_abioticKbasal$fundNiche) <- paste0("Y",1:S)
  SIMlist$GLV_abioticKbasal$fundNiche <- SIMlist$GLV_abioticKbasal$fundNiche[,spNewNames]
}

if ("rickerKbasal" %in% simMethods & file.exists(paste0(simPath, "ricker.fundNicheTh.abioticKbasal.csv"))){
  SIMlist$Ricker_abioticKbasal$fundNiche <- read.csv2(paste0(simPath, "ricker.fundNicheTh.abioticKbasal.csv"), row.names=1)
  colnames(SIMlist$Ricker_abioticKbasal$fundNiche) <- paste0("Y",1:S)
  SIMlist$Ricker_abioticKbasal$fundNiche <- SIMlist$Ricker_abioticKbasal$fundNiche[,spNewNames]
}


nameOrder = cumsum(table(spNewNames)[spNewNames])[paste0("Y", 1:S)]  # table to alternate between name ordering



plotDistributions(SIMlist$GLV_abioticGR, plotprey=T, RN=T, main="Predicted and observed species distribution: GLV Growth Rates",filename=paste0(figPath,"RealizedNicheLV_GR.pdf"))
ggsave(filename=paste0(figPath,"RealizedNicheLV_GR.pdf"),width=15,height=10)

plotDistributions(SIMlist$GLV_abioticGR, plotprey=F, FN=T, main="Predicted and observed species distribution: GLV Growth Rates",filename=paste0(figPath,"FundNicheLV_GR.pdf"))



if ("glvKbasal" %in% simMethods){
  plotDistributions(SIMlist$GLV_abioticKbasal, plotprey=T, RN=T,main="Predicted and observed species distribution: GLV Basal")
  ggsave(filename=paste0(figPath,"RealizedNicheLV_basalK.pdf"))
  
  plotDistributions(SIMlist$GLV_abioticKbasal,plotprey=F, FN=T, main="Predicted and observed species distribution: GLV Basal")
  ggsave(filename=paste0(figPath,"FundNicheLV_basalK.pdf"))
  
}


if ("rickerKbasal" %in% simMethods){
  plotDistributions(SIMlist$Ricker_abioticKbasal,  plotprey=T, RN=T, main="Predicted and observed species distribution: Ricker Basal")
  ggsave(filename=paste0(figPath,"RealizedNicheRicker_Kbasal.pdf"))
  
  plotDistributions(SIMlist$Ricker_abioticKbasal,  plotprey=T, FN=T, main="Predicted and observed species distribution: Ricker Basal")
  ggsave(filename=paste0(figPath,"FundamentalNicheRicker_Kbasal.pdf"))
  
}


if ("soiER" %in% simMethods){
  plotDistributions(SIMlist$SOI_abioticER, plotprey=T, RN=T,main="Predicted and observed species distribution: SOI Extinction Rates")
  ggsave(filename=paste0(figPath,"FundamentalNicheSOIER.pdf"))
  
}


if ("soiERbasal" %in% simMethods){
  plotDistributions(SIMlist$SOI_abioticERbasal, plotprey=T, RN=T,main="Predicted and observed species distribution: SOI Basal")
  ggsave(filename=paste0(figPath,"FundamentalNicheSOIERbasal.pdf"))
  
}



if ("vc" %in% simMethods){
  plotDistributions(SIMlist$VC_abiotic, plotprey=T, RN=T, main="Predicted and observed species distribution: VirtualCom")
  ggsave(filename=paste0(figPath,"FundamentalNicheVirtualCom.pdf"))
  
}


############# Recovery of the realized niche

for (i in 1:length(SIMlist)){
  SIM = SIMlist[[i]]
  envs = SIM$X[1:nEnv,ifelse(intercept,-1,-(K+1))]
  
  # Realized niche
  realizedNiche = SIM$prob
  
  # Predicted realized niche: trophic SDM
  proba = list()
  if (!is.null(SIM$table_no_err_pr_prob_bin)){
    proba$pred = matrix(SIM$table_no_err_pr_prob_bin$pred, ncol=S)[1:51,]
    
    #SIMlist[[i]]$realizedNicheError = colMeans(abs(realizedNiche - proba$pred))
    SIMlist[[i]]$realizedNicheError = sapply(1:S, function(s)transport::wasserstein1d(realizedNiche[,s], proba$pred[,s]))
    SIMlist[[i]]$realizedNicheModeError = abs(envs[apply(realizedNiche, 2, which.max)]-envs[apply(proba$pred, 2, which.max)])
  }
  
  # Predicted realized niche: non-trophic SDM
  if (!is.null(SIM$SDM.table_no_err_pr_prob_bin)){
    proba$SDM.pred = matrix(SIM$SDM.table_no_err_pr_prob_bin$pred, ncol=S)[1:51,]
    
    #SIMlist[[i]]$SDM.realizedNicheError = colMeans(abs(realizedNiche - proba$SDM.pred))
    SIMlist[[i]]$SDM.realizedNicheError = sapply(1:S, function(s)transport::wasserstein1d(realizedNiche[,s], proba$SDM.pred[,s]))
    SIMlist[[i]]$SDM.realizedNicheModeError = abs(envs[apply(realizedNiche, 2, which.max)]-envs[apply(proba$SDM.pred, 2, which.max)])
  }
}

for (i in which(sapply(SIMlist, function(SIM)!is.null(SIM$realizedNicheError)))){
  
  assign(paste0("t",i), data.frame(niche=c(SIMlist[[i]]$realizedNicheError,SIMlist[[i]]$SDM.realizedNicheError,SIMlist[[i]]$R2["R2_bin",],SIMlist[[i]]$SDM.R2),
                                   name=names(SIMlist)[i],type=c(rep("Wasserstein d.",2*S),rep("R2",2*S)),
                                   focal=spNewNames,
                                   model=rep(c(rep("tSDM",S),rep("SDM",S)),2),TrophLev=paste0("TL",trophL[nameOrder])))
  assign(paste0("p_all",i), ggplot(data=get(paste0("t",i))[which(get(paste0("t",i))$TrophLev!="TL1"),]) + geom_boxplot(aes(y=niche,x=model,col=model)) + ggtitle(get(paste0("t",i))$name) + facet_wrap(.~type,scale="free")  + ylab("Distance"))
  assign(paste0("p_diff",i), ggplot(data=get(paste0("t",i))) + geom_boxplot(aes(y=niche,x=model,col=model)) + ggtitle(get(paste0("t",i))$name) + facet_wrap(TrophLev~type,scale="free",nrow=3)  + ylab("Distance"))
  
}

main="Boxplot Niche error"
p_all=eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p_all",1:length(SIMlist)),collapse=","),"),nrow=2,top=)")))
ggsave(filename = paste0(figPath,"BoxplotRealizes_ALL.pdf"))

main="Boxplot Niche error per trophic level"
p_diff=eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p_diff",1:length(SIMlist)),collapse=","),"),nrow=2,top=)")))
ggsave(filename = paste0(figPath,"BoxplotRealizes_ALL.pdf"))

p_all1
ggsave(filename = paste0(figPath,"BoxplotRealizes_LV_GRall.pdf"))

### Fig de Jeremy

# Wasserstein

#t1$focal=as.numeric(sub("Y","",t1$focal))
t1$focal <- factor(t1$focal, levels=paste0("Y",(c(1:S))))
p1=ggplot(data=t1[which(t1$type=="Wasserstein d."),]) + geom_point(aes(x=focal,y=niche,col=model,shape=model),size=4)+
    geom_line(aes(x=focal,y=niche),arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "closed"))+
    ylab("Wass. distance from the niche") + theme_classic() +theme(legend.position="top")
ggsave(filename = paste(figPath,"Wass_all1.pdf"),p1,width=10)
p2= ggplot(data=t1[which(t1$type=="Wasserstein d." & t1$TrophLev!="TL1"),]) + geom_boxplot(aes(y=niche,col=model))  +
    theme_classic() +theme(legend.position="top")
p2
ggsave(filename = paste(figPath,"Wass_all2.pdf"))

# R2

#t1$focal=as.numeric(sub("Y","",t1$focal))
t1$focal <- factor(t1$focal, levels=paste0("Y",(c(1:S))))
p1=ggplot(data=t1[which(t1$type=="R2"),]) + geom_point(aes(x=focal,y=niche,col=model,shape=model))+
  geom_line(aes(x=focal,y=niche),arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "closed")) + ylab("R2")
ggsave(filename = paste(figPath,"R2_all1.pdf"),p1)
p2= ggplot(data=t1[which(t1$type=="R2"),]) + geom_boxplot(aes(y=niche,col=model)) + facet_wrap(.~TrophLev)
p2
ggsave(filename = paste(figPath,"R2_all2.pdf"))






for (i in which(sapply(SIMlist, function(SIM)!is.null(SIM$realizedNicheError)))){
  SIM = SIMlist[[i]]
  Ymax = ceiling(max(SIM$SDM.realizedNicheError)*20)/20
  plot(SIM$realizedNicheError, pch=16, cex=2.5, col="lightcoral", ylim=c(0,1.2*Ymax), cex.main=2.5, cex.lab=2, cex.axis=2, 
       xaxt='n', xlab="", ylab="Wasserstein distance", main=paste("RN error:", names(SIMlist)[i]))
  abline(v=cumsum(Stroph)[-L]+0.5, lty=2)
  points(SIM$SDM.realizedNicheError, col="cornflowerblue", pch=17, cex=2.5)
  abs.diff = abs((SIM$realizedNicheError) - SIM$SDM.realizedNicheError)
  rel.diff = abs.diff/SIM$SDM.realizedNicheError
  for (s in which(abs.diff/Ymax>0.03 & is.finite(rel.diff)))arrows(x0=s,x1=s,y0=SIM$SDM.realizedNicheError[s],y1=(SIM$realizedNicheError)[s], length=rel.diff[s]**0.3/7)
  legend(x=2, y=1.15*Ymax, c("SDM","tSDM"), pch=c(17,16), pt.cex=3, col=c("cornflowerblue","lightcoral"), x.intersp=0.5, y.intersp=-1, cex=2, text.width=1, box.lty=0, ncol=2)
  text(x=cumsum(c(1.5,Stroph[-max(trophL)])), y=1.15*Ymax, labels=paste0("TL",1:max(trophL)), cex=2, col="gray30")
  title(xlab="Species", line=1.5, cex.lab=2)
  
  boxplot(list(Basal.species=SIM$SDM.realizedNicheError[trophL==1], Higher.trophic.levels=SIM$SDM.realizedNicheError[trophL>1]), boxwex=0.7, col="cornflowerblue", ylim=c(0,1.2*Ymax), cex.main=2.5, cex.lab=2, cex.axis=2, 
          ylab="Wasserstein distance", main=paste("RN error:", names(SIMlist)[i]))
  boxplot(list(Basal.species=(SIM$realizedNicheError)[trophL==1], Higher.trophic.levels=(SIM$realizedNicheError)[trophL>1]), boxwex=0.2, cex.axis=2, col="lightcoral", add=T)
  legend(x=0.7, y=1.15*Ymax, c("SDM","tSDM"), lty=1, col=c("cornflowerblue","lightcoral"), x.intersp=0.5, y.intersp=-1, cex=2, text.width=0.2, seg.len=1, lwd=10, box.lty=0, ncol=2)
}
invisible(dev.copy(png, file=paste0(simPath, "Wasserstein_distance_RealizedNiche_nTrain", round(nRep*trainFraction), ".png"), width=12, height=4*LEN, unit="in", res=150))
invisible(dev.off())











options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 150)
predAccuracy = do.call(rbind, lapply(SIMlist, function(SIM)reshape::melt(rbind(SIM$AUC, SIM$TSS), varnames=c("Prediction.Quality", "Species"))))
predAccuracy$simMethod <- gsub("\\..*", "", rownames(predAccuracy))
predAccuracy$Data.type <- gsub(".*_", "", predAccuracy$Prediction.Quality)
predAccuracy$Metrics <- gsub("_.*", "", predAccuracy$Prediction.Quality)
predAccuracy$trophL <- as.factor(trophL[as.numeric(gsub("Y", "", predAccuracy$Species))])
predAccuracy$relative.value <- predAccuracy$value/rep(predAccuracy[predAccuracy$Data.type=="simul",]$value, each=3)

predAccuracy.SDM = do.call(rbind, lapply(SIMlist, function(SIM)reshape::melt(rbind(SIM$SDM.AUC, SIM$SDM.TSS), varnames=c("Prediction.Quality", "Species"))))
predAccuracy.SDM$simMethod <- gsub("\\..*", "", rownames(predAccuracy.SDM))
predAccuracy.SDM$Data.type <- gsub(".*_", "", predAccuracy.SDM$Prediction.Quality)
predAccuracy.SDM$Metrics <- gsub("_.*", "", predAccuracy.SDM$Prediction.Quality)
predAccuracy.SDM$trophL <- as.factor(trophL[as.numeric(gsub("Y", "", predAccuracy.SDM$Species))])
predAccuracy.SDM$relative.value <- predAccuracy.SDM$value/rep(predAccuracy.SDM[predAccuracy.SDM$Data.type=="simul",]$value, each=2)



options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 150)

predAccuracy.combined = rbind(cbind(predAccuracy, Model="tSDM"), cbind(predAccuracy.SDM, Model="SDM"))
p <- ggplot(predAccuracy.combined, aes(x=trophL, y=sapply(relative.value, min, 1), fill=Model, shape=Data.type)) +
  #geom_point() +
  #geom_violin() +
  geom_boxplot() + scale_fill_brewer(palette="Pastel1", direction=-1) +
  xlab("Trophic levels") + ylab("Relative values compared to ideal simulation data") +
  ggtitle("AUC and TSS comparison between prediction methods, relative to actual presence probabilities") +
  facet_wrap(Metrics ~ simMethod, nrow=2, scales="free") + theme(legend.position="bottom")
p
ggsave(paste0(figPath,"AUC_TSS_Comparison.pdf"))

# 
# 
# options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 150)
# 
# predAccuracy.combined = rbind(cbind(predAccuracy, Model="tSDM"), cbind(predAccuracy.SDM, Model="SDM"))
# p <- ggplot(predAccuracy.combined[predAccuracy.combined$Data.type=="bin",],
#             aes(x=trophL, y=sapply(relative.value, min, 1), fill=Model)) +
#   #geom_point() +
#   #geom_violin() +
#   geom_boxplot() + scale_fill_brewer(palette="Pastel1", direction=-1) +
#   xlab("Trophic levels") + ylab("Relative values compared to ideal simulation data") + 
#   ggtitle("AUC and TSS comparison between prediction methods, relative to actual presence probabilities") +
#   facet_wrap(Metrics ~ simMethod, nrow=2, scales="free") + theme(legend.position="bottom")
# p
# ggsave(filename = paste0(simPath, "Prediction_accuracies_nTrain", round(nRep*trainFraction), ".png"), p, height=5, width=10)


#########################################################################################################
########## fundamental niche

knownFundNiche = sapply(SIMlist, function(X)"fundNiche" %in% names(X))
for (i in which(knownFundNiche)){
  
  SIM = SIMlist[[i]]
  envs = SIM$X[1:nEnv,ifelse(intercept,-1,-(K+1))]
  
  # Fundamental niche
  # SIM$fundNiche
  
  # Predicted fundamental niche: trophic SDM
  coeff = list()
  if (length(grep("brmsfit",class(SIM$m_stan$model[[1]])))>0){  # brms package
    coeff$abiotic = sapply(SIM$m_stan$model, function(sp)posterior_summary(sp)[1:(intercept+K*(2-linear)),"Estimate",drop=FALSE])
    coeff$biotic = sapply(SIM$m_stan$model, function(sp)posterior_summary(sp)[-c(1:(intercept+K*(2-linear)),nrow(posterior_summary(sp))),"Estimate",drop=FALSE])
  }else {  # rstan package
    coeff$abiotic = sapply(SIM$m_stan$model, function(sp)coef(sp)[1:(intercept+K*(2-linear))])
    coeff$biotic = sapply(SIM$m_stan$model, function(sp)coef(sp)[-(1:(intercept+K*(2-linear)))])}
  
  coeff$infered.preys.sum = sapply(1:S, function(s)sum(coeff$biotic[[paste0("Y",s)]][which(as.numeric(gsub(".*Y","",names(coeff$biotic[[paste0("Y",s)]])))<s)], na.rm=T))
  
  a.tSDM = coeff$abiotic[3,]
  b.tSDM = coeff$abiotic[2,]
  c.tSDM = coeff$abiotic[1,]

  predFundNiche = plogis(envs**2%*%t(a.tSDM) + envs%*%t(b.tSDM) + rep(1,nEnv)%*%t(c.tSDM) + rep(1,nEnv)%*%t(coeff$infered.preys.sum))
 
  SIMlist[[i]]$fundNicheError = sapply(1:S, function(s)transport::wasserstein1d(SIM$fundNiche[,s],predFundNiche[,s]))
  SIMlist[[i]]$fundNicheModeError = abs(niche_optima-(0:100/100)[apply((0:100/100)**2%*%t(a.tSDM) + (0:100/100)%*%t(b.tSDM) + rep(1,101)%*%t(c.tSDM), 2, which.max)])
  SIMlist[[i]]$fundNicheR2 = sapply(1:S, function(s)cor(SIM$fundNiche[,s],predFundNiche[,s])^2)
  # Predicted fundamental niche: non-trophic SDM
  if (!is.null(SIM$SDM_stan)) {
    coeff$SDM.abiotic = sapply(SIM$SDM_stan$model, function(sp)coef(sp)[1:(intercept+K*(2-linear))])
  }
  
  a.SDM = coeff$SDM.abiotic[3,]
  b.SDM = coeff$SDM.abiotic[2,]
  c.SDM = coeff$SDM.abiotic[1,]
  #SIMlist[[i]]$SDM.fundNicheError = colMeans(abs(SIM$fundNiche -
  #     plogis(envs**2%*%t(a.SDM) + envs%*%t(b.SDM) + rep(1,nEnv)%*%t(c.SDM))))
  predFundNiche.SDM = plogis(envs**2%*%t(a.SDM) + envs%*%t(b.SDM) + rep(1,nEnv)%*%t(c.SDM))
  SIMlist[[i]]$SDM.fundNicheError = sapply(1:S, function(s)transport::wasserstein1d(SIM$fundNiche[,s],predFundNiche.SDM[,s]))
  SIMlist[[i]]$SDM.fundNicheModeError = abs(niche_optima-(0:100/100)[apply((0:100/100)**2%*%t(a.SDM) + (0:100/100)%*%t(b.SDM) + rep(1,101)%*%t(c.SDM), 2, which.max)])
  SIMlist[[i]]$SDM.fundNicheR2 = sapply(1:S, function(s)cor(SIM$fundNiche[,s],predFundNiche.SDM[,s])^2)
  
}



for (i in which(sapply(SIMlist, function(SIM)!is.null(SIM$fundNicheError)))){
  
  # assign(paste0("t",i), data.frame(niche=c(SIMlist[[i]]$fundNicheError,SIMlist[[i]]$SDM.fundNicheError,SIMlist[[i]]$fundNicheModeError,SIMlist[[i]]$SDM.fundNicheModeError),
  #                                  name=names(SIMlist)[i],type=c(rep("Wasserstein distance",2*S),rep("Niche optimum error",2*S)),
  #                                  model=rep(c(rep("tSDM",S),rep("SDM",S)),2),TrophLev=paste0("TL",trophL[nameOrder])))
  assign(paste0("t",i), data.frame(niche=c(SIMlist[[i]]$fundNicheError,SIMlist[[i]]$SDM.fundNicheError,SIMlist[[i]]$fundNicheModeError,SIMlist[[i]]$SDM.fundNicheModeError,SIMlist[[i]]$fundNicheR2,SIMlist[[i]]$SDM.fundNicheR2),
                                   focal=spNewNames,
                                   name=names(SIMlist)[i],type=c(rep("Wasserstein distance",2*S),rep("Niche optimum error",2*S),rep("R2",2*S)),
                                   model=rep(c(rep("tSDM",S),rep("SDM",S)),3),TrophLev=paste0("TL",trophL[nameOrder])))
  
  
  assign(paste0("p_all",i), ggplot(data=get(paste0("t",i))[which(get(paste0("t",i))$TrophLev!="TL1"),]) + geom_boxplot(aes(y=niche,x=model,col=model)) + ggtitle(get(paste0("t",i))$name) + facet_wrap(.~type,scale="free")  + ylab("Distance"))
  assign(paste0("p_diff",i), ggplot(data=get(paste0("t",i))) + geom_boxplot(aes(y=niche,x=model,col=model)) + ggtitle(get(paste0("t",i))$name) + facet_wrap(TrophLev~type,scale="free",nrow=3)  + ylab("Distance"))
  
}

main="Boxplot Niche error"
p_all=eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p_all",1:length(SIMlist)),collapse=","),"),nrow=2,top=)")))
ggsave(filename = paste0(figPath,"BoxplotFund_ALL.pdf"))

main="Boxplot Niche error per trophic level"
p_diff=eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p_diff",1:length(SIMlist)),collapse=","),"),nrow=2,top=)")))
ggsave(filename = paste0(figPath,"BoxplotFund_ALL.pdf"))

p_all1
ggsave(filename = paste0(figPath,"BoxplotFund_LV_GRall.pdf"))


###### Figure Jeremy

#t1$focal=as.numeric(sub("Y","",t1$focal))
t1$focal <- factor(t1$focal, levels=paste0("Y",(c(1:S))))
p1=ggplot(data=t1[which(t1$type=="Wasserstein distance"),]) + geom_point(aes(x=focal,y=niche,col=model,shape=model),size=4)+
  geom_line(aes(x=focal,y=niche),arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "closed"))+
  ylab("Wass. distance from the niche") + theme_classic()+theme(legend.position="top")
ggsave(filename = paste(figPath,"Wass_fund_all1.pdf"),p1,width=10)
p2= ggplot(data=t1[which(t1$type=="Wasserstein distance" & t1$TrophLev!="TL1"),]) + geom_boxplot(aes(y=niche,col=model)) +theme_classic()+theme(legend.position = "top")
p2
ggsave(filename = paste(figPath,"Wass_fund2.pdf"))

# R2

#t1$focal=as.numeric(sub("Y","",t1$focal))
t1$focal <- factor(t1$focal, levels=paste0("Y",(c(1:S))))
p1=ggplot(data=t1[which(t1$type=="R2"),]) + geom_point(aes(x=focal,y=niche,col=model,shape=model),size=4)+
  geom_line(aes(x=focal,y=niche),arrow = arrow(length=unit(0.1,"cm"), ends="first", type = "closed"))
ggsave(filename = paste(figPath,"R2_all1",p1))
p2= ggplot(data=t1[which(t1$type=="R2"),]) + geom_boxplot(aes(y=niche,col=model)) + facet_wrap(.~TrophLev)
p2
ggsave(filename = paste(figPath,"R2_all2.pdf"))















# Relationship with the strength of the biotic control (for LV GR only)
SIM=SIMlist[[1]]
delta_R2=SIM$R2["R2_bin",which(trophL[nameOrder]!=1)]-SIM$SDM.R2[trophL[nameOrder]!=1]
delta_WRN=SIM$realizedNicheError[trophL[nameOrder]!=1]-SIM$SDM.realizedNicheError[trophL[nameOrder]!=1]
delta_WFN=SIM$fundNicheError[trophL[nameOrder]!=1]-SIM$SDM.fundNicheError[trophL[nameOrder]!=1]
delta_ModeFN=SIM$fundNicheModeError[trophL[nameOrder]!=1]-SIM$SDM.fundNicheModeError[trophL[nameOrder]!=1]


IntMat2=IntMat;IntMat2[lower.tri(PredMat, diag=TRUE)] <- 0

BioticPress=colSums(IntMat2)[which(trophL[nameOrder]!=1)]

lmR2=lm(delta_R2~BioticPress)
summary(lmR2)

lmWRN=lm(delta_WRN~BioticPress)
summary(lmWRN)

lmWFN=lm(delta_WFN~BioticPress)
summary(lmWFN)

lmMFN=lm(delta_ModeFN~BioticPress)
summary(lmMFN)

