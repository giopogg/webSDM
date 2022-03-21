####### Jeremy: test of tSDM and SDMs on simulation data
rm(list=ls())
library(igraph)
library(Matrix)
library(cheddar)
library(igraph)
library(GGally)
library(intergraph)
library(gridExtra)
library(dismo)
library(coda)
library(cheddar)
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

root= "/Users/poggiatg/Documents/GitHub/trophicSDM/VirtualEcoSim/"

setwd(root)

source("../tSDM_functions.R")
  
#################################################################################################
#### Test on simulated data
################################################################################################

#job=args[1]
job=NA
simPath = "Simulations_S20L3_nEnv51_nRep50_strengthBI5_asy/"   # directory from which to load the simulations
S = as.numeric(gsub(".*S|L.*", "", simPath))                       # number of species
L = as.numeric(gsub(".*L|_nEnv.*", "", simPath))                   # number of trophic levels
nEnv = as.numeric(gsub(".*nEnv|_nRep.*", "", simPath))             # number of environments
nRep = as.numeric(gsub(".*nRep|_strengthBI.*", "", simPath))         # number of replicates
strengthBI = as.numeric(gsub(".*strengthBI|_asy.*", "", simPath)) # strength of biotic interactions
mergeReplicates = TRUE   # whether to merge the simulation replicates to train the models
nbMerge = 1                              # between 0 and 1, ratio of replicate to keep for fitting
linear = F                               # whether to include a quadratic term for X
fitPreds = F                             # whether to fit predators
if(!linear) poly=F                       # if yes, whether to take a raw polynomial (poly=T) or an orthogonal one (poly=F)

horsh=F

iter=200
pred_samples=200

figPath=paste0(simPath,"Fig/")
# Gather outputs in a single list
SIMlist = list()

# Check the simulation methods available (among "glvGR", "glvKbasal", "rickerKbasal", "soiER", "soiERbasal", "vc")
simMethods = gsub("_finalStates_abiotic|.csv", "", grep("finalStates", list.files(simPath), value = TRUE))
#simMethods = simMethods[1]




### loadData : the function that loads the simulated community datasets.

loadData <- function(filePath){
  df.merged <- read.csv2(filePath, check.names = FALSE)
  df.list <- lapply(unique(df.merged$datasets), function(i) {
    df <- subset(df.merged, datasets == i, select = -c(datasets))
    rownames(df) <- df$rowN
    return(df[-(1:2)])
  })
  return(df.list)
}

### finalStatesToXY : the function that converts the simulated datasets into environment (X) and 
#                     presence/absence (Y) matrices, for training and testing.


finalStatesToXY <- function(finalStates, nbMerge=1){
  # Training dataset
  finalStates.combineRep = Reduce(cbind, finalStates[1:round(length(finalStates)*nbMerge)])
  Y = as.data.frame(t(finalStates.combineRep))
  colnames(Y) = paste("Y", 1:ncol(Y), sep="")
  X = cbind(1, X1=as.numeric(colnames(finalStates.combineRep)))
  XY.train = list(X=X,Y=Y)
  
  return(list(XY.train=XY.train))
}

checkCV = function(out = out, K = 5, nEnv){
  
  X = out$X[1:nEnv,"X1"]
  partition = rep(1:K,each = round(nEnv/K))
  S = length(ncol(out$Y))
  
  if(length(partition)<nEnv) partition = c(partition, rep(K,nEnv-length(partition)))
  
  index_all = vector(length=nrow(out$X))
  
  for(i in 1:K){
    
    index_all[ which(out$X[,"X1"] %in% X[partition==i])] = i
    
  }
  
  # check for all training datasets
  temp_check = vector(length = K)
  for(i in 1:K){
    # check if all species have at least one zero and one 1
    temp_check[i] = any(apply(out$Y[index_all != i,],2,function(x) var(x) ==0))
    
  }
  return(!any(temp_check))
}



### Build the known graph

IntMat <- as.matrix(read.csv2(paste0(simPath, "InteractionMatrix.csv"), row.names = 1))  # matrix of species interactions
PredMat <- sign(IntMat) ; PredMat[lower.tri(PredMat, diag=TRUE)] <- 0                    # matrix of trophic links

spNames = colnames(PredMat)                                                         # species names
trophL = as.numeric(sub(pattern = "Sp.*TL", replacement = "", x = spNames))         # species trophic levels
Stroph = table(trophL)                                                              # number of species by trophic level
community <- Community(nodes=data.frame(node=spNames),
                       trophic.links=PredationMatrixToLinks(PredMat),
                       properties=list(title='Random community'))                   # cheddar community
dimnames(PredMat) <- list(paste0("Y", 1:S),paste0("Y", 1:S))

# Build the graph of trophic interactions from PredMat
G = graph_from_adjacency_matrix(t(PredMat))
#spNewNames = names(topo_sort(G,mode="in"))
spNewNames = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=T)]$name
save(spNewNames, file = paste0(simPath,"spNewNames.R"))


#options(repr.plot.width = 12, repr.plot.height = 3, repr.plot.res = 150)
#ggnet2(intergraph::asNetwork(G),arrow.gap = 0.05,arrow.size = 10,label=TRUE)+ggtitle("Known graph")


# Load niche optima used in the simulations
niche_optima <- read.csv2(paste0(simPath, "niche_optima.csv"))[[1]]
names(niche_optima) = paste0("Y", 1:S)
niche_optima = niche_optima[spNewNames]
niche_optima


K = ifelse(is.vector(niche_optima), 1, ncol(niche_optima))    # number of environmental variables

# Create formulas for the environmental part of each species
if(linear){
  env.form=as.list(rep("~X1",S))
}else{
  if(!poly){
    env.form=as.list(rep("~X1+I(X1^2)",S))
  }else{
    env.form=as.list(rep("~poly(X1,2)",S))
  }
}
env.form=lapply(env.form,FUN=as.formula)
names(env.form) = paste("Y", 1:S, sep="")



#Composite variables
#Instead of looking individual species interactions, integrate several species in a single composite variables, for example richness or diversity.
#sp.formula = "richness+richness^2"
#sp.formula = "richness"
#sp.partition = list(c("Y1","Y2"),c("Y3","Y4","Y5","Y6"))


#### GLV GR

if ("glvGR" %in% simMethods){
  glv.finalStates.abioticGR <- loadData(paste0(simPath, "glv_finalStates_abioticGR.csv"))
  head(glv.finalStates.abioticGR[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(glv.finalStates.abioticGR, nbMerge = nbMerge)
  }else{out = finalStatesToXY(glv.finalStates.abioticGR[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$GLV_abioticGR = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$GLV_abioticGR$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                   if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$GLV_abioticGR$survivalProba <- data.frame(t(Reduce("+", glv.finalStates.abioticGR)/nRep))
  }else{
    simMethods = simMethods[-which(simMethods=="glvGR")]
  }
}

#### GLV GR K basal

  if("glvKbasal" %in% simMethods){
  glv.finalStates.abioticKbasal <- loadData(paste0(simPath, "glv_finalStates_abioticKbasal.csv"))
  if(nRep != length(glv.finalStates.abioticKbasal)) warning("Wrong number of replicates: ", length(glv.finalStates.abioticKbasal), " instead of ", nRep)
  head(glv.finalStates.abioticKbasal[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(glv.finalStates.abioticKbasal, nbMerge = nbMerge)
  }else{out = finalStatesToXY(glv.finalStates.abioticKbasal[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$GLV_abioticKbasal = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$GLV_abioticKbasal$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                       if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$GLV_abioticKbasal$survivalProba <- data.frame(t(Reduce("+", glv.finalStates.abioticKbasal)/nRep))
  }else{
    simMethods = simMethods[-which(simMethods=="glvKbasal")]
  }
}


### Ricker model

if ("rickerKbasal" %in% simMethods){
  ricker.finalStates.abioticKbasal <- loadData(paste0(simPath, "ricker_finalStates_abioticKbasal.csv"))
  if(nRep != length(ricker.finalStates.abioticKbasal)) warning("Wrong number of replicates: ", length(ricker.finalStates.abioticKbasal), " instead of ", nRep)
  head(ricker.finalStates.abioticKbasal[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(ricker.finalStates.abioticKbasal, nbMerge = nbMerge)
  }else{out = finalStatesToXY(ricker.finalStates.abioticKbasal[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$Ricker_abioticKbasal = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$Ricker_abioticKbasal$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                          if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$Ricker_abioticKbasal$survivalProba <- data.frame(t(Reduce("+", ricker.finalStates.abioticKbasal)/nRep))
  }else{
    simMethods = simMethods[-which(simMethods=="rickerKbasal")]
  }
}


### SOI ER
if ("soiER" %in% simMethods){
  soi.finalStates.abioticER <- loadData(paste0(simPath, "soi_finalStates_abioticER.csv"))
  if(nRep != length(soi.finalStates.abioticER)) warning("Wrong number of replicates: ", length(soi.finalStates.abioticER), " instead of ", nRep)
  head(soi.finalStates.abioticER[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(soi.finalStates.abioticER, nbMerge = nbMerge)
  }else{out = finalStatesToXY(soi.finalStates.abioticER[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$SOI_abioticER = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$SOI_abioticER$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                   if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$SOI_abioticER$survivalProba <- data.frame(t(Reduce("+", soi.finalStates.abioticER)/nRep))
  }else{
    simMethods = simMethods[-which(simMethods=="soiER")]
  }
}

### SOR ER basal
if ("soiERbasal" %in% simMethods){
  soi.finalStates.abioticERbasal <- loadData(paste0(simPath, "soi_finalStates_abioticERbasal.csv"))
  if(nRep != length(soi.finalStates.abioticERbasal)) warning("Wrong number of replicates: ", length(soi.finalStates.abioticERbasal), " instead of ", nRep)
  head(soi.finalStates.abioticERbasal[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(soi.finalStates.abioticERbasal, nbMerge = nbMerge)
  }else{out = finalStatesToXY(soi.finalStates.abioticERbasal[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$SOI_abioticERbasal = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$SOI_abioticERbasal$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                        if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$SOI_abioticERbasal$survivalProba <- data.frame(t(Reduce("+", soi.finalStates.abioticERbasal)/nRep))
  }else{
    simMethods = simMethods[-which(simMethods=="soiERbasal")]
  }
}



### Virtual COM


if ("vc" %in% simMethods){
  vc.finalStates.abiotic <- loadData(paste0(simPath, "vc_finalStates_abiotic.csv"))
  if(nRep != length(vc.finalStates.abiotic)) warning("Wrong number of replicates: ", length(vc.finalStates.abiotic), " instead of ", nRep)
  head(vc.finalStates.abiotic[[1]])
  
  # Format the abiotic conditions (X) and biotic presence/absence (Y) matrices
  if (mergeReplicates){out = finalStatesToXY(vc.finalStates.abiotic, nbMerge = nbMerge)
  }else{out = finalStatesToXY(vc.finalStates.abiotic[1])}
  
  check_data = checkCV(out$XY.train, K = 5, nEnv = nEnv)
  
  if(check_data){
  SIMlist$VC_abiotic = out$XY.train
  
  # Construct an analogue of the coefficients matrix
  SIMlist$VC_abiotic$B <- cbind(rep(c(1,-1),c(Stroph[1],S-Stroph[1])),
                                if (linear) 10*(niche_optima-0.5) else niche_optima)
  
  # Estimate survival probabilities
  SIMlist$VC_abiotic$survivalProba <- data.frame(t(Reduce("+", vc.finalStates.abiotic)/nRep))
  colnames(SIMlist$VC_abiotic$survivalProba) <- spNames
  }else{
    simMethods = simMethods[-which(simMethods=="vc")]
  }
}


#save(SIMlist,file = paste0(simPath,"SIMlist.Rdata"))

#####################################################################################################
########## Fit tSDM


# Choose the inference algorithms
#algos = c("glm", "stan", "bayes")
algos = c( "stan")

### GLM
if ("glm" %in% algos){
  for (i in 1:length(SIMlist)){
    print(names(SIMlist[i]))
    SIM = SIMlist[[i]]
    
    # Basic GLM
    #SIMlist[[i]]$m_glm = trophicSDM(Y=SIM$Y, X=SIM$X, G=G, formulas=env.form, sp.formula=sp.formula, 
    #                                sp.partition=sp.partition, penal=NULL, method="glm", 
    #                                family=binomial(link="logit"), fitPreds=fitPreds, iter=iter)
    
    # Basic GLM + constraint on coefficients signs (+ for preys->predators and - for predators->preys)
    SIMlist[[i]]$m_glm = trophicSDM(Y=SIM$Y, X=SIM$X, G=G, formulas=env.form, sp.formula=sp.formula, 
                                    sp.partition=sp.partition, penal="coeff.signs", method="glm", 
                                    family=binomial(link="logit"), fitPreds=fitPreds, iter=iter)
  }
}


## Stan

if ("stan" %in% algos){
  rstan_options(auto_write = TRUE)
  for (i in 1:length(SIMlist)){
    print(names(SIMlist[i]))
    SIM = SIMlist[[i]]
    
    # STAN GLM
    if(!horsh){
    SIMlist[[i]]$m_stan = trophicSDM(Y=SIM$Y,X=SIM$X,G=G,formulas=env.form,penal=NULL,method="stan_glm",
                                     family=binomial(link = "logit"),fitPreds=fitPreds,iter=iter,run.parallel = F)
    }else{
      SIMlist[[i]]$m_stan = trophicSDM(Y=SIM$Y,X=SIM$X,G=G,formulas=env.form,penal="horshoe",method="stan_glm",
                                       family=binomial(link = "logit"),fitPreds=fitPreds,iter=iter,run.parallel = TRUE)
      
    }
    # STAN GLM + constraint on coefficients signs (+ for preys->predators and - for predators->preys)
    #SIMlist[[i]]$m_stan = trophicSDM(Y=SIM$Y, X=SIM$X, G=G, formulas=env.form, sp.formula=sp.formula,
                                     #sp.partition=sp.partition, penal="coeff.signs", method="stan_glm",
                                     #family=bernoulli(link = "logit"), fitPreds=fitPreds, iter=iter)
  }
}

#save(SIMlist,file = paste0(figPath,"model.RData"))

##################################################################################################
### Check model convergence (here done only with one species from the stan model and GLV GR)

for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))
  SIM = SIMlist[[i]]
  Rhats = unlist(lapply(1:S, function(s){
                        rhat(SIM$m_stan$model[[s]]) 
                }))
  Rhats = data.frame(names=names(Rhats),value=Rhats, Bioabio = as.numeric((1:length(Rhats)) %in% grep("Y",names(Rhats))))
  
  ness = unlist(lapply(1:S, function(s){
    neff_ratio(SIM$m_stan$model[[s]]) 
  }))
  
  ness = data.frame(names=names(ness),value=ness, Bioabio = as.numeric((1:length(ness)) %in% grep("Y",names(ness))))
  
  table=rbind(data.frame(Rhats,type="Rhats"),data.frame(ness,type="ness"))
  
  ggplot(data= table,aes(x=value,color=as.factor(Bioabio))) + geom_histogram(alpha=0.8)+
    scale_color_discrete(name="Biotic or abiotic",labels=c("Abiotic","Biotic")) +
    theme_classic() + ggtitle(names(SIMlist[i])) + facet_wrap(.~type,scale="free")
  
  ggsave(filename = paste0(figPath, names(SIMlist[i]),"_Convergence.pdf"))
  
}

##################################################################################################
# Analyse results
#Obtain posterior means and 95% credible regions for the niche optimum and interaction parameters of each species, and compare against the true ones


# p estimate finds the intercept, the niche optima and the biotic variables
for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))
  SIM = SIMlist[[i]]
  
  estimates_stan=estimates_glm=estimates_bayes=data.frame( p.est=double(),
                                                           est.97=double(),
                                                           est.02=double(),
                                                           true=double(),
                                                           type=factor(levels = c("biotic","abiotic")),
                                                           sp.id=factor(),
                                                           sp.trophL=factor())
  
  for (j in 1:S){
    if(fitPreds){
    neigh.sp = which((IntMat-diag(diag(IntMat)))[,paste0("Sp",j,".TL",trophL[j])]!=0)   # neighbour species in the trophic web
    }else{
      IntMat_temp = IntMat
      IntMat_temp [lower.tri(IntMat_temp , diag=TRUE)] = 0  
      neigh.sp = which((IntMat_temp-diag(diag(IntMat_temp)))[,paste0("Sp",j,".TL",trophL[j])]!=0)   # neighbour species in the trophic web
      
    }
    if ("glm" %in% algos){
      print("glm")
      model_j = SIM$m_glm$model[[paste0("Y",j)]]
      post_j = coef(model_j)
      print(list(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", row.names(post_j))]),
                 est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], warning=function(cond){message(paste0("Warning GLM Y",j));return(NA)}, error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
                 est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], warning=function(cond)return(NA), error=function(cond)return(NA)) else NA,
                 true=c(SIM$B[paste0("Y",j),],IntMat[neigh.sp,j]),
                 type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(neigh.sp)))),
                 sp.id=as.factor(j),
                 sp.trophL=paste("Trophic level", trophL[j])))
      estimates_glm=rbind(estimates_glm,data.frame(#p.est=if(linear)post_j else c(post_j[1], -(post_j[2])/(2*post_j[3]), post_j[-(1:3)]),
        #p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[which(suppressWarnings(as.numeric(gsub("Y","",rownames(post_j)), drop = T))<j)]),
        p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", row.names(post_j))]),
        est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], warning=function(cond){message(paste0("Warning GLM Y",j));return(NA)}, error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
        est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], warning=function(cond)return(NA), error=function(cond)return(NA)) else NA,
        true=c(SIM$B[paste0("Y",j),],IntMat[neigh.sp,j]),
        type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(neigh.sp)))),
        sp.id=as.factor(j),
        sp.trophL=paste("Trophic level", trophL[j])
      ))}
    
    if ("stan" %in% algos){
      #print("stan")
      model_j = SIM$m_stan$model[[paste0("Y",j)]]
      post_j = as.data.frame(posterior_summary(model_j))
      #post_j = post_j[-nrow(post_j),] #only for gaussian data!!!
      # print(list(p.est=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[which(gsub(".*_", "", rownames(post_j)) %in% names(neighbors(G,paste0("Y", j),mode=c("out"))))]),
      #            p.est2=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[grepl("Y", row.names(post_j))]),
      #            est.02=if(linear)post_j$Q2.5 else NA,
      #            est.97=if(linear)post_j$Q97.5 else NA,
      #            true=c(SIM$B[paste0("Y",j),],IntMat[neigh.sp,j]),
      #            type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(neigh.sp)))),
      #            sp.id=as.factor(j),
      #            sp.trophL=paste("Trophic level", trophL[j])))
      #estimates_stan=rbind(estimates_stan,data.frame(p.est=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[which(gsub(".*_", "", rownames(post_j)) %in% names(neighbors(G,paste0("Y", j),mode=c("out"))))]),
      estimates_stan=rbind(estimates_stan,data.frame(p.est=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[grepl("Y", row.names(post_j))]),
                                                     est.02=if(linear)post_j$Q2.5 else NA,
                                                     est.97=if(linear)post_j$Q97.5 else NA,
                                                     true=c(SIM$B[paste0("Y",j),],IntMat[neigh.sp,j]),
                                                     type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(neigh.sp)))),
                                                     sp.id=as.factor(j),
                                                     sp.trophL=paste("Trophic level", trophL[j])
      ))}
    
    
    if ("bayes" %in% algos){
      print("bayes")
      prey.sp = which((IntMat-diag(diag(IntMat)))[,paste0("Sp",j,".TL",trophL[j])] > 0)   # keep only positive interactions
      #print(prey.sp)
      model_j = SIM$m_bayes$model[[paste0("Y",j)]]
      post_j = coef(model_j)
      #print(post_j)
      #print(post_j[grepl("Y", names(post_j))])
      # print(list(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", names(post_j))]),
      #            est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
      #            est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], error=function(cond)return(NA)) else NA,
      #            true=c(SIM$B[paste0("Y",j),],IntMat[prey.sp,j]),
      #            type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(prey.sp)))),
      #            sp.id=as.factor(j),
      #            sp.trophL=paste("Trophic level", trophL[j])))
      estimates_bayes=rbind(estimates_bayes,data.frame(#p.est=if(linear)post_j else c(post_j[1], -(post_j[2])/(2*post_j[3]), post_j[-(1:3)]),
        #p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[which(suppressWarnings(as.numeric(gsub("Y","",names(post_j)), drop = T))<j)]),
        p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", names(post_j))]),
        est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
        est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], error=function(cond)return(NA)) else NA,
        true=c(SIM$B[paste0("Y",j),],IntMat[prey.sp,j]),
        type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(prey.sp)))),
        sp.id=as.factor(j),
        sp.trophL=paste("Trophic level", trophL[j])
      ))}
  }
  SIMlist[[i]]$estimates_stan = estimates_stan
  SIMlist[[i]]$estimates_glm = estimates_glm
  SIMlist[[i]]$estimates_bayes = estimates_bayes
}


ggplotEstimates <- function(SIM, algos=c("stan", "glm", "bayes"), shapeLegend=FALSE, colorLegend=FALSE, wrap.type=TRUE, wrap.trophL=TRUE, onlyBio=FALSE, onlyAbio=FALSE){
  # Remove intercepts, for which there is no true value to compare to
  intercepts=abs(SIM[[paste0("estimates_",algos[1])]]$true)==1
  if (onlyBio){cond=!intercepts & SIM[[paste0("estimates_",algos[1])]]$type=="biotic"
  }else{if (onlyAbio){cond=!intercepts & SIM[[paste0("estimates_",algos[1])]]$type=="abiotic"
  }else cond=!intercepts}
  
  if ("stan" %in% algos){
    # Plot pred vs true for stan_glm
    p_stan=ggplot(data=SIM$estimates_stan[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) + 
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("stan_glm()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_stan <- p_stan + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_stan <- p_stan + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_stan <- p_stan + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  if ("glm" %in% algos){
    # Plot pred vs true for glm
    p_glm=ggplot(data=SIM$estimates_glm[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) +
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type)) + 
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("glm()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_glm <- p_glm + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_glm <- p_glm + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_glm <- p_glm + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  if ("bayes" %in% algos){
    # Plot pred vs true for bayes
    p_bayes=ggplot(data=SIM$estimates_bayes[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) +
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("bayes()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_bayes <- p_bayes + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_bayes <- p_bayes + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_bayes <- p_bayes + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  return(list(p_glm=if("glm" %in% algos){p_glm}else NA,
              p_stan=if("stan" %in% algos){p_stan}else NA, 
              p_bayes=if("bayes" %in% algos){p_bayes}else NA))
}


# plot together with the trophic network
for(i in 1:length(SIMlist)){
  ggplotEstimates(SIMlist[[i]], algos=algos)$p_stan + ggtitle(label=paste("True vs inferred regression coefficients:", simMethods[i]))
  ggsave(filename = paste0(figPath, "Biotic_abiotic_coefficients_",simMethods[i],"nTrain", round(nRep*nbMerge), ".pdf"))
}



# Compute the correlation coefficients
if(S>=20){
  R2_bin.abiotic = sapply(SIMlist,function(SIM){
    cond = SIM$estimates_stan$type=="abiotic" & abs(SIM$estimates_stan$true)!=1 #only abiotic
    TL = as.numeric(sub("Trophic level ", "", SIM$estimates_stan$sp.trophL))
    return(list(TL1 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==1])!=0, cor(SIM$estimates_stan$p.est[cond & TL==1],SIM$estimates_stan$true[cond & TL==1], method="spearman"),NA),
                TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman"),NA),
                TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman"),NA)))
  })

  
  R2_bin.biotic = sapply(SIMlist,function(SIM){
    cond = SIM$estimates_stan$type=="biotic" & abs(SIM$estimates_stan$true)!=1
    TL = as.numeric(sub("Trophic level ", "", SIM$estimates_stan$sp.trophL))
    if(fitPreds){
    return(list(TL1 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==1])!=0, cor(SIM$estimates_stan$p.est[cond & TL==1],SIM$estimates_stan$true[cond & TL==1], method="spearman"),NA),
                TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman"),NA),
                TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman"),NA)))
    }else{
      return(list(TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman"),NA),
                  TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman"),NA)))
    }
    
  })
  
  R2_bin_pVal.abiotic = sapply(SIMlist,function(SIM){
    cond = SIM$estimates_stan$type=="abiotic" & abs(SIM$estimates_stan$true)!=1
    TL = as.numeric(sub("Trophic level ", "", SIM$estimates_stan$sp.trophL))
    return(list(TL1 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==1])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==1],SIM$estimates_stan$true[cond & TL==1], method="spearman",exact=FALSE)$p.value,NA),
                TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman",exact=FALSE)$p.value,NA),
                TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman",exact=FALSE)$p.value,NA)))
  })
  R2_bin_pVal.biotic = sapply(SIMlist,function(SIM){
    cond = SIM$estimates_stan$type=="biotic" & abs(SIM$estimates_stan$true)!=1
    TL = as.numeric(sub("Trophic level ", "", SIM$estimates_stan$sp.trophL))
    if(fitPreds){
    return(list(TL1 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==1])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==1],SIM$estimates_stan$true[cond & TL==1], method="spearman",exact=FALSE)$p.value,NA),
           TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman",exact=FALSE)$p.value,NA),
           TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman",exact=FALSE)$p.value,NA)))
    }else{
      return(list(TL2 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==2])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==2],SIM$estimates_stan$true[cond & TL==2], method="spearman",exact=FALSE)$p.value,NA),
                  TL3 = ifelse(var(SIM$estimates_stan$p.est[cond & TL==3])!=0, cor.test(SIM$estimates_stan$p.est[cond & TL==3],SIM$estimates_stan$true[cond & TL==3], method="spearman",exact=FALSE)$p.value,NA)))
      
    }
    })
  
  options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 150)
  df.abiotic <- cbind(type=rownames(R2_bin.abiotic),stack(as.data.frame(R2_bin.abiotic)),stack(as.data.frame(R2_bin_pVal.abiotic))$values)
  df.biotic <- cbind(type=rownames(R2_bin.biotic),stack(as.data.frame(R2_bin.biotic)),stack(as.data.frame(R2_bin_pVal.biotic))$values)
  names(df.abiotic) <- names(df.biotic) <- c("Trophic.level", "Correlation.coef", "Simulation.model", "p.value")
  p.max = 0.2
  p1 <- ggplot(df.abiotic, aes(Simulation.model, Trophic.level)) + 
    geom_tile(data = subset(df.abiotic, p.value<=p.max), aes(fill=Correlation.coef), na.rm=T) +
    geom_tile(data = subset(df.abiotic, p.value>p.max), fill=NA, na.rm=T) + 
    geom_text(aes(label = ifelse(p.value<=p.max, paste("rho =", round(Correlation.coef, 2), "\np =", round(p.value,2)), 
                                 paste("p > ", p.max))), na.rm=T) +
    scale_fill_gradient2(low="red", mid="white", high="green4", limits=c(-1,1)) +
    theme(panel.background = element_rect(fill=alpha("blue", 0.1),colour=alpha("blue", 0.1),size = 0.5), panel.grid = element_blank()) +
    ggtitle("Correlations between inferred and 'true' ABIOTIC parameters, tSDM")
  p2 <- ggplot(df.biotic, aes(Simulation.model, Trophic.level)) + 
    geom_tile(data = subset(df.biotic, p.value<=p.max), aes(fill=Correlation.coef), na.rm=T) +
    geom_tile(data = subset(df.biotic, p.value>p.max), fill=NA, na.rm=T) + 
    geom_text(aes(label = ifelse(p.value<=p.max, paste("rho =", round(Correlation.coef, 2), "\np =", round(p.value,2)), 
                                 paste("p > ", p.max))), na.rm=T) +
    scale_fill_gradient2(low="red", mid="white", high="green4", limits=c(-1,1)) +
    theme(panel.background = element_rect(fill = alpha("red",0.1),colour = alpha("red",0.1),size = 0.5), panel.grid = element_blank()) +
    ggtitle("Correlations between inferred and 'true' BIOTIC parameters, tSDM")
  p <- grid.arrange(arrangeGrob(grobs = list(p1,p2), layout_matrix=matrix(rep(1:2,c(nrow(R2_bin.abiotic),nrow(R2_bin.biotic))))), ncol=1)
  ggsave(filename = paste0(figPath, "Correlation_coefficients_nTrain", round(nRep*nbMerge), ".png"), p, height=5, width=10)
  invisible(p)
}


##################################################################################################################
#### Predictions
#Let's predict probabilities of presence in the environments of the testing dataset in order 
#to evaluate the prediction accuracy.

################################################
### With probcov=F

# Compute mean and 95%CI predictions, R2, calibration (% of samples inside the CI)
for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))

  Xnew=SIMlist[[i]]$X[1:nEnv,] #Xnew are the environmental variables that we used to fit the model (actually if Xnew=NULL, this is the default choice)
  error_prop_sample=10 #error_prop_sample is an important parameter. The greater it is, the largest the error propagation (and therefore the predictive credible regions). see tSDM_functions for more details
  pred_samples=pred_samples

  # report the probabilities of presence (realized niches) obtained from the observed dataset
  SIMlist[[i]]$prob = SIMlist[[i]]$survivalProba
  colnames(SIMlist[[i]]$prob) = paste0("Y", 1:S)

  # run predictions
  SIMlist[[i]]$pred_stan_bin=trophicSDM_predict(m=SIMlist[[i]]$m_stan,Xnew=Xnew,binary.resp=F,prob.cov=F,mode=ifelse(fitPreds,"all","out"),pred_samples=pred_samples,error_prop_sample=error_prop_sample)

  # reorder the columns of prob to be sure that we are comparing the same things
  SIMlist[[i]]$prob = SIMlist[[i]]$prob[,spNewNames]
  SIMlist[[i]]$Y=SIMlist[[i]]$Y[,spNewNames]

  p.mean.stan.temp=lapply(SIMlist[[i]]$pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,mean, MARGIN = 1))
  SIMlist[[i]]$p.mean.stan_bin=do.call(cbind, p.mean.stan.temp)
  SIMlist[[i]]$p.qinf.stan_bin = do.call(cbind,lapply(SIMlist[[i]]$pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.025)))
  SIMlist[[i]]$p.qsup.stan_bin = do.call(cbind,lapply(SIMlist[[i]]$pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.975)))
  
  
  # Compute R2 and calibration
  R2_bin=sapply(spNewNames,function(x){cor(SIMlist[[i]]$p.mean.stan_bin[,x],SIMlist[[i]]$prob[,x])^2})
  SIMlist[[i]]$R2_bin = R2_bin
  
  SIMlist[[i]]$calibration_bin = apply(sapply(spNewNames, function(name){
    sapply(1:nEnv,function(x) ifelse(SIMlist[[i]]$prob[x,name] < SIMlist[[i]]$p.qsup.stan_bin[x,name] & SIMlist[[i]]$prob[x,name]>SIMlist[[i]]$p.qinf.stan_bin[x,name], 1, 0))
  }
  ),FUN=mean,MARGIN=2)
  
}



for (i in 1:length(SIMlist)){
  SIM = SIMlist[[i]]
   
  SIMlist[[i]]$table_pred_bin=data.frame(obs=as.vector(as.matrix(SIM$prob)), pred=as.vector(as.matrix(SIM$p.mean.stan_bin)),
                                                       est.02=as.vector(as.matrix(SIMlist[[i]]$p.qinf.stan_bin)), est.97=as.vector(as.matrix(SIMlist[[i]]$p.qsup.stan_bin)),
                                                       sp.name=rep(colnames(SIM$prob),each=nrow(SIM$prob)),
                                                       trophL=trophL[as.numeric(sub("Y", "", rep(colnames(SIM$prob),each=nrow(SIM$prob))))])
  
  #plot predicted vs true probabilities of presence                                             
  SIMlist[[i]]$p.predictions_bin <- ggplot(data=SIMlist[[i]]$table_pred_bin,mapping=aes(x=obs,y=pred,col=factor(sp.name, levels=unique(sp.name)))) + 
    geom_abline(slope=1,intercept = 0) + guides(col=guide_legend(title=NULL, nrow=10)) +
    geom_point(alpha=0.5) + geom_linerange(mapping=aes(ymin=est.02, ymax=est.97), alpha=0.5) +
    ggtitle(paste0("Stan-GLM predictions for ", names(SIMlist)[i], " simulations"))+
    xlab("Observed presence probability") + ylab("Predicted presence probability") + 
    facet_wrap(trophL~., labeller="label_both") + theme_minimal() + theme(legend.position="none")
}

#plot predicted vs true probabilities of presence 
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 150)
p <- grid.arrange(arrangeGrob(grobs=lapply(SIMlist, function(SIM)SIM$p.predictions_bin)))
ggsave(filename = paste0(figPath, "Presence_probability_predictions_Binary", round(nRep*nbMerge), ".pdf"), p, height=7, width=12)




####### probcov=T
#What happens when we use probabilities of presence for predictions? -> set prob.cov = T
for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))
  
  SIMlist[[i]]$pred_stan_prob=trophicSDM_predict(m=SIMlist[[i]]$m_stan,Xnew=Xnew,binary.resp=F,prob.cov=T,mode=ifelse(fitPreds,"all","out"),pred_samples=pred_samples,error_prop_sample=error_prop_sample)
  
  p.mean.stan.temp=lapply(SIMlist[[i]]$pred_stan_prob$sp.prediction,FUN=function(x) apply(x$predictions.prob,mean, MARGIN = 1))
  SIMlist[[i]]$p.mean.stan_prob=do.call(cbind, p.mean.stan.temp)
  SIMlist[[i]]$p.qinf.stan_prob = do.call(cbind,lapply(SIMlist[[i]]$pred_stan_prob$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.025)))
  SIMlist[[i]]$p.qsup.stan_prob = do.call(cbind,lapply(SIMlist[[i]]$pred_stan_prob$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.975)))
  
  
  R2_prob=sapply(spNewNames,function(x){cor( SIMlist[[i]]$p.mean.stan_prob[,x], SIMlist[[i]]$prob[,x])^2})
  SIMlist[[i]]$R2_prob = R2_prob

  SIMlist[[i]]$calibration_prob = apply(sapply(spNewNames, function(name){
    sapply(1:nEnv,function(x) ifelse(SIMlist[[i]]$prob[x,name] < SIMlist[[i]]$p.qsup.stan_prob[x,name] & SIMlist[[i]]$prob[x,name]>SIMlist[[i]]$p.qinf.stan_prob[x,name], 1, 0))
   }
  ),FUN=mean,MARGIN=2)
  
  
  SIMlist[[i]]$corr_prob_bin = sapply(spNewNames,function(x){cor( SIMlist[[i]]$p.mean.stan_prob[,x], SIMlist[[i]]$p.mean.stan_bin[,x])^2})

}




for (i in 1:length(SIMlist)){
  SIM = SIMlist[[i]]
  
  SIMlist[[i]]$table_pred_prob=data.frame(obs=as.vector(as.matrix(SIM$prob)), pred=as.vector(as.matrix(SIM$p.mean.stan_prob)),
                                          est.02=as.vector(as.matrix(SIMlist[[i]]$p.qinf.stan_prob)), est.97=as.vector(as.matrix(SIMlist[[i]]$p.qsup.stan_prob)),
                                          sp.name=rep(colnames(SIM$prob),each=nrow(SIM$prob)),
                                          trophL=trophL[as.numeric(sub("Y", "", rep(colnames(SIM$prob),each=nrow(SIM$prob))))])
  
  
  #plot predicted vs true probabilities of presence                                             
  SIMlist[[i]]$p.predictions_prob <- ggplot(data=SIMlist[[i]]$table_pred_prob,mapping=aes(x=obs,y=pred,col=factor(sp.name, levels=unique(sp.name)))) + 
    geom_abline(slope=1,intercept = 0) + guides(col=guide_legend(title=NULL, nrow=10)) +
    geom_point(alpha=0.5) + geom_linerange(mapping=aes(ymin=est.02, ymax=est.97), alpha=0.5) +
    ggtitle(paste0("Stan-GLM predictions for ", names(SIMlist)[i], " simulations"))+
    xlab("Observed presence probability") + ylab("Predicted presence probability") + 
    facet_wrap(trophL~., labeller="label_both") + theme_minimal() + theme(legend.position="none")
}

#plot predicted vs true probabilities of presence 
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 150)
p <- grid.arrange(arrangeGrob(grobs=lapply(SIMlist, function(SIM)SIM$p.predictions_prob)))
ggsave(filename = paste0(figPath, "Presence_probability_predictions_Prob", round(nRep*nbMerge), ".pdf"), p, height=7, width=12)




#### Compare the results``


options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 150)
predAccuracy = do.call(rbind, lapply(SIMlist, function(SIM)reshape::melt(rbind(bin=SIM$R2_bin,prob=SIM$R2_prob), varnames=c("Prediction.Quality", "Species"))))


predAccuracy$simMethod <- gsub("\\..*", "", rownames(predAccuracy))
predAccuracy$Data.type <- gsub(".*_", "", predAccuracy$Prediction.Quality)
predAccuracy$Metrics <- gsub("_.*", "", predAccuracy$Prediction.Quality)
predAccuracy$trophL <- as.factor(trophL[as.numeric(gsub("Y", "", predAccuracy$Species))])

p <- ggplot(predAccuracy, aes(x=trophL, y=value, fill=Data.type)) +
  #geom_point() +
  #geom_violin() +
  geom_boxplot() +
  facet_wrap(.~simMethod, nrow=2) + theme(legend.position="bottom") +
  ggtitle("R2 comparison between prediction methods, relative to actual presence probabilities")
p
ggsave(paste0(figPath,"R2_tSDM.pdf"))





################################################################################################
##### Fit Classical SDM
################################################################################################

# Define a graph without biotic interaction to infer based only on environmental variables
G_null = graph_from_adjacency_matrix(matrix(0, nrow=S, ncol=S, dimnames=c(list(spNewNames),list(spNewNames))))
#options(repr.plot.width = 12, repr.plot.height = 3, repr.plot.res = 150)
#ggnet2(intergraph::asNetwork(G_null),arrow.gap = 0.05,arrow.size = 10,label=TRUE)

# glm
if ("glm" %in% algos){
  for (i in 1:length(SIMlist)){
    print(names(SIMlist[i]))
    SIM = SIMlist[[i]]
    SIMlist[[i]]$SDM_glm = trophicSDM(Y=SIM$Y,X=SIM$X,G=G_null,formulas=env.form,penal=NULL,method="glm",family=binomial(link = "logit"),iter=iter)
  }
}


if ("stan" %in% algos){
  for (i in 1:length(SIMlist)){
    print(names(SIMlist[i]))
    SIM = SIMlist[[i]]
    SIMlist[[i]]$SDM_stan = trophicSDM(Y=SIM$Y,X=SIM$X,G=G_null,formulas=env.form,penal=NULL,method="stan_glm",
                                       family=binomial(link = "logit"),iter=iter,run.parallel = F)
  }
}


if ("bayes" %in% algos){
  for (i in 1:length(SIMlist)){
    print(names(SIMlist[i]))
    SIM = SIMlist[[i]]
    SIMlist[[i]]$SDM_bayes = trophicSDM(Y=SIM$Y,X=SIM$X,G=G_null,formulas=env.form,penal=NULL,method="bayesglm",
                                        family=binomial(link = "logit"),iter=iter)
  }
}

###########################################
## Analyse results



for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))
  SIM = SIMlist[[i]]
  
  SDM.estimates_stan=SDM.estimates_glm=SDM.estimates_bayes=data.frame( p.est=double(),
                                                                       est.97=double(),
                                                                       est.02=double(),
                                                                       true=double(),
                                                                       type=factor(levels = c("biotic","abiotic")),
                                                                       sp.id=factor(),
                                                                       sp.trophL=factor())
  
  for (j in 1:S){
    if ("glm" %in% algos){
      print("glm")
      model_j = SIM$SDM_glm$model[[paste0("Y",j)]]
      post_j = coef(model_j)
      print(list(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", row.names(post_j))]),
                 est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], warning=function(cond){message(paste0("Warning GLM Y",j));return(NA)}, error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
                 est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], warning=function(cond)return(NA), error=function(cond)return(NA)) else NA,
                 true=SIM$B[paste0("Y",j),],
                 type=as.factor(rep("abiotic",K)),
                 sp.id=as.factor(j),
                 sp.trophL=paste("Trophic level", trophL[j])))
      SDM.estimates_glm=rbind(SDM.estimates_glm,data.frame(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])], post_j[grepl("Y", row.names(post_j))]),
                                                           est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], warning=function(cond){message(paste0("Warning GLM Y",j));return(NA)}, error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
                                                           est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], warning=function(cond)return(NA), error=function(cond)return(NA)) else NA,
                                                           true=SIM$B[paste0("Y",j),],
                                                           type=as.factor(rep("abiotic",K)),
                                                           sp.id=as.factor(j),
                                                           sp.trophL=paste("Trophic level", trophL[j])
      ))}
    
    if ("stan" %in% algos){
      print("stan")
      model_j = SIM$SDM_stan$model[[paste0("Y",j)]]
      post_j = as.data.frame(posterior_summary(model_j))
      #post_j = post_j[-nrow(post_j),]
      print(list(p.est=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[grepl("Y", row.names(post_j))]),
                 est.02=if(linear)post_j$Q2.5 else NA,
                 est.97=if(linear)post_j$Q97.5 else NA,
                 true=SIM$B[paste0("Y",j),],
                 type=as.factor(rep("abiotic",K)),
                 sp.id=as.factor(j),
                 sp.trophL=paste("Trophic level", trophL[j])))
      SDM.estimates_stan=rbind(SDM.estimates_stan,data.frame(p.est=if(linear)post_j$Estimate else c(post_j$Estimate[1], (0:100/100)[which.max(post_j$Estimate[3]*(0:100/100)**2 + post_j$Estimate[2]*0:100/100 + post_j$Estimate[1])], post_j$Estimate[grepl("Y", row.names(post_j))]),
                                                             est.02=if(linear)post_j$Q2.5 else NA,
                                                             est.97=if(linear)post_j$Q97.5 else NA,
                                                             true=SIM$B[paste0("Y",j),],
                                                             type=as.factor(rep("abiotic",K)),
                                                             sp.id=as.factor(j),
                                                             sp.trophL=paste("Trophic level", trophL[j])
      ))}
    
    if ("bayes" %in% algos){
      print("bayes")
      model_j = SIM$SDM_bayes$model[[paste0("Y",j)]]
      post_j = coef(model_j)
      print(list(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])]),
                 est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
                 est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], error=function(cond)return(NA)) else NA,
                 true=SIM$B[paste0("Y",j),],
                 type=as.factor(rep("abiotic",K)),
                 sp.id=as.factor(j),
                 sp.trophL=paste("Trophic level", trophL[j])))
      SDM.estimates_bayes=rbind(SDM.estimates_bayes,data.frame(p.est=if(linear)post_j else c(post_j[1], (0:100/100)[which.max(post_j[3]*(0:100/100)**2 + post_j[2]*0:100/100 + post_j[1])]),
                                                               est.02=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,1], error=function(cond){message(paste0("Error bayes Y",j));return(NA)}) else NA,
                                                               est.97=if(linear)tryCatch(confint(profile(model_j),prob = 0.95)[,2], error=function(cond)return(NA)) else NA,
                                                               true=SIM$B[paste0("Y",j),],
                                                               type=as.factor(rep("abiotic",K)),
                                                               sp.id=as.factor(j),
                                                               sp.trophL=paste("Trophic level", trophL[j])
      ))}
  }
  SIMlist[[i]]$SDM.estimates_stan = SDM.estimates_stan
  SIMlist[[i]]$SDM.estimates_glm = SDM.estimates_glm
  SIMlist[[i]]$SDM.estimates_bayes = SDM.estimates_bayes
}




ggplotEstimates <- function(SIM, algos=c("stan", "glm", "bayes"), shapeLegend=FALSE, colorLegend=FALSE, wrap.type=TRUE, wrap.trophL=TRUE, onlyBio=FALSE, onlyAbio=FALSE){
  # Remove intercepts, for which there is no true value to compare to
  intercepts=abs(SIM[[paste0("estimates_",algos[1])]]$true)==1
  if (onlyBio){cond=!intercepts & SIM[[paste0("estimates_",algos[1])]]$type=="biotic"
  }else{if (onlyAbio){cond=!intercepts & SIM[[paste0("estimates_",algos[1])]]$type=="abiotic"
  }else cond=!intercepts}
  
  if ("stan" %in% algos){
    # Plot pred vs true for stan_glm
    p_stan=ggplot(data=SIM$estimates_stan[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) + 
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("stan_glm()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_stan <- p_stan + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_stan <- p_stan + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_stan <- p_stan + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  if ("glm" %in% algos){
    # Plot pred vs true for glm
    p_glm=ggplot(data=SIM$estimates_glm[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) +
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type)) + 
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("glm()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_glm <- p_glm + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_glm <- p_glm + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_glm <- p_glm + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  if ("bayes" %in% algos){
    # Plot pred vs true for bayes
    p_bayes=ggplot(data=SIM$estimates_bayes[cond,],aes(x=true,y=p.est,colour=type,shape=sp.id))+ 
      geom_point(cex=1.5) +
      #scale_shape_manual(values=c(1:2,4:(S+1))) +
      scale_shape_manual(values=rep(19,S)) + scale_colour_brewer(palette="Set1", direction=-1) +
      #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
      geom_abline(intercept=0, slope=1, colour=alpha("black",0.5), linetype="dashed")+
      #coord_fixed(ratio=1, xlim=NULL, ylim=NULL, expand=TRUE, clip="on")+ 
      ggtitle("bayes()") + xlab("Truth proxy") + ylab("Inferred") + 
      theme_classic() + theme(panel.grid = element_blank()) + 
      guides(shape=if(shapeLegend){guide_legend(nrow=5)}else{FALSE}, 
             color=if(colorLegend){guide_legend()}else{FALSE})
    
    if (wrap.type){
      if (wrap.trophL){p_bayes <- p_bayes + facet_wrap(type~sp.trophL, scales="free", nrow=ifelse(onlyBio|onlyAbio,L,2))
      }else{p_bayes <- p_bayes + facet_wrap(.~type, scales="free", nrow=L)}
    }else{if (wrap.trophL){p_bayes <- p_bayes + facet_wrap(sp.trophL~., scales="free", nrow=L)}}
  }
  
  return(list(p_glm=if("glm" %in% algos){p_glm}else NA,
              p_stan=if("stan" %in% algos){p_stan}else NA, 
              p_bayes=if("bayes" %in% algos){p_bayes}else NA))
}



# Compute the correlation coefficients
if (S>=20){
  SDM.R2_bin = sapply(SIMlist,function(SIM){
    cond = SIM$SDM.estimates_stan$type=="abiotic" & abs(SIM$SDM.estimates_stan$true)!=1
    TL = as.numeric(sub("Trophic level ", "", SIM$SDM.estimates_stan$sp.trophL))
    return(list(TL1 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==1])!=0, cor(SIM$SDM.estimates_stan$p.est[cond & TL==1],SIM$SDM.estimates_stan$true[cond & TL==1], method="spearman"),NA),
                TL2 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==2])!=0, cor(SIM$SDM.estimates_stan$p.est[cond & TL==2],SIM$SDM.estimates_stan$true[cond & TL==2], method="spearman"),NA),
                TL3 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==3])!=0, cor(SIM$SDM.estimates_stan$p.est[cond & TL==3],SIM$SDM.estimates_stan$true[cond & TL==3], method="spearman"),NA)))
    
  })
  SDM.R2_bin_pVal = sapply(SIMlist,function(SIM){
    cond = SIM$SDM.estimates_stan$type=="abiotic" & abs(SIM$SDM.estimates_stan$true)!=1
    TL = as.numeric(sub("Trophic level ", "", SIM$SDM.estimates_stan$sp.trophL))
    return(list(TL1 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==1])!=0, cor.test(SIM$SDM.estimates_stan$p.est[cond & TL==1],SIM$SDM.estimates_stan$true[cond & TL==1], method="spearman",exact=FALSE)$p.value,NA),
                TL2 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==2])!=0, cor.test(SIM$SDM.estimates_stan$p.est[cond & TL==2],SIM$SDM.estimates_stan$true[cond & TL==2], method="spearman",exact=FALSE)$p.value,NA),
                TL3 = ifelse(var(SIM$SDM.estimates_stan$p.est[cond & TL==3])!=0, cor.test(SIM$SDM.estimates_stan$p.est[cond & TL==3],SIM$SDM.estimates_stan$true[cond & TL==3], method="spearman",exact=FALSE)$p.value,NA)))
    
     })
  
  options(repr.plot.width = 12, repr.plot.height = 3, repr.plot.res = 150)
  df <- cbind(type=rownames(SDM.R2_bin),stack(as.data.frame(SDM.R2_bin)),stack(as.data.frame(SDM.R2_bin_pVal))$values)
  names(df) <- c("Trophic.level", "Correlation.coef", "Simulation.model", "p.value")
  p.max = 0.2
  p <- ggplot(df, aes(Simulation.model, Trophic.level)) + 
    geom_tile(data = subset(df, p.value<=p.max), aes(fill=Correlation.coef), na.rm=T) +
    geom_tile(data = subset(df, p.value>p.max), fill=NA, na.rm=T) + 
    geom_text(aes(label = ifelse(p.value<=p.max, paste("rho =", round(Correlation.coef, 2), "\np =", round(p.value,2)), 
                                 paste("p > ", p.max))), na.rm=T) +
    scale_fill_gradient2(low="red", mid="white", high="green4", limits=c(-1,1)) +
    theme(panel.background = element_rect(fill=alpha("blue", 0.1),colour=alpha("blue", 0.1),size = 0.5), panel.grid = element_blank()) +
    ggtitle("Correlations between inferred and 'true' ABIOTIC parameters, SDM")
  ggsave(filename = paste0(figPath, "Correlation_coefficients_SDM_nTrain", round(nRep*nbMerge), ".png"), p, height=5, width=10)
  #p
}


######################################################################################################
################# Prediction

### With SDM pred.cov=F or T is the same

for (i in 1:length(SIMlist)){
  print(names(SIMlist[i]))
  
  Xnew=SIMlist[[i]]$X[1:nEnv,] #Xnew are the environmental variables that we used to fit the model (actually if Xnew=NULL, this is the default choice) 
  error_prop_sample=10 #error_prop_sample is an important parameter. The greater it is, the largest the error propagation (and therefore the predictive credible regions). see tSDM_functions for more details
  pred_samples=pred_samples
  
  #run predictions
  SIMlist[[i]]$SDM.pred_stan_bin=trophicSDM_predict(m=SIMlist[[i]]$SDM_stan,Xnew=Xnew,binary.resp=F,prob.cov=F,pred_samples=pred_samples,error_prop_sample=error_prop_sample)
  
  p.mean.stan.temp=lapply(SIMlist[[i]]$SDM.pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,mean, MARGIN = 1))
  SIMlist[[i]]$SDM.p.mean.stan_bin=do.call(cbind, p.mean.stan.temp)
  
  # Compute CI
  SIMlist[[i]]$SDM.p.qinf.stan_bin = do.call(cbind,lapply(SIMlist[[i]]$SDM.pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.025)))
  SIMlist[[i]]$SDM.p.qsup.stan_bin = do.call(cbind,lapply(SIMlist[[i]]$SDM.pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.975)))
  
  # Compute R2
  R2_bin=sapply(spNewNames,function(x){cor(SIMlist[[i]]$SDM.p.mean.stan_bin[,x],SIMlist[[i]]$prob[,x])^2})
  SIMlist[[i]]$SDM.R2_bin = R2_bin
  
  SIMlist[[i]]$SDM.calibration_bin = apply(sapply(spNewNames, function(name){
    sapply(1:nEnv,function(x) ifelse(SIMlist[[i]]$prob[x,name] < SIMlist[[i]]$SDM.p.qsup.stan_bin[x,name] & SIMlist[[i]]$prob[x,name]>SIMlist[[i]]$SDM.p.qinf.stan_bin[x,name], 1, 0))
  }
  ),FUN=mean,MARGIN=2)
  
}


for (i in 1:length(SIMlist)){
  
  SIM = SIMlist[[i]]
 
  SIMlist[[i]]$SDM.table_pred_bin=data.frame(obs=as.vector(as.matrix(SIM$prob)), pred=as.vector(as.matrix(SIM$SDM.p.mean.stan_bin)),
                                                       est.02=as.vector(as.matrix(SIMlist[[i]]$SDM.p.qinf.stan_bin)), est.97=as.vector(as.matrix(SIMlist[[i]]$SDM.p.qsup.stan_bin)),
                                                       sp.name=rep(colnames(SIM$prob),each=nrow(SIM$prob)),
                                                       trophL=trophL[as.numeric(sub("Y", "", rep(colnames(SIM$prob),each=nrow(SIM$prob))))])
  
  #plot predicted vs true probabilities of presence                                             
  SIMlist[[i]]$SDM.p.predictions_bin <- ggplot(data=SIMlist[[i]]$SDM.table_pred_bin,mapping=aes(x=obs,y=pred,col=factor(sp.name, levels=unique(sp.name)))) + 
    geom_abline(slope=1,intercept = 0) + guides(col=guide_legend(title=NULL, nrow=10)) +
    geom_point(alpha=0.5) + geom_linerange(mapping=aes(ymin=est.02, ymax=est.97), alpha=0.5) +
    ggtitle(paste0("Stan-GLM predictions for ", names(SIMlist)[i], " simulations"))+
    xlab("Observed presence probability") + ylab("Predicted presence probability") + 
    facet_wrap(trophL~., labeller="label_both") + theme_minimal() + theme(legend.position="none")

}

#plot predicted vs true probabilities of presence 
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 150)
p <- grid.arrange(arrangeGrob(grobs=lapply(SIMlist, function(SIM)SIM$SDM.p.predictions_bin)))
ggsave(filename = paste0(figPath, "Presence_probability_predictions_SDM_nTrain", round(nRep*nbMerge), ".png"), p, height=7, width=12)


## Check AUS TSS
#We compute AUC and TSS in predictions the AUC and TSS of the simulated data (are probabilities of presence well separated from presence/absences)
figPath=simPath

SIMlist = lapply(SIMlist, function(SIM){SIM[!grepl("m_glm|m_bayes|SDM_glm|SDM_bayes|pEstim", names(SIM))]})
#lapply(SIMlist$GLV_abioticGR, function(SIM)pryr::object_size(SIM))
#save(SIMlist, file=paste0(figPath, "SIMlist.Rdata"))




###################################################################################################
### Load fundamental niches

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


# Compute the fundamental niche for tSDM (this does not depend on if prob.cov=T or F)
# The fundamental niche for SDM = realised niche!
for(i in 1:length(SIMlist)){
  
  SIM=SIMlist[[i]]
  intercept=T
  
  SIMlist[[i]]$pFund.mean.stan = SIMlist[[i]]$pFund.qinf.stan = SIMlist[[i]]$pFund.qsup.stan = matrix(NA, nrow=nEnv,ncol=S,
                                                                                       dimnames = list(NULL,names(SIM$m_stan$model)))
  
  for(s in 1:S){
    
    sp_mod = SIM$m_stan$model[[s]]

    newdata =  cbind(X1=SIM$X[1:nEnv,"X1"],
                     data.frame(matrix(1,nrow=nEnv,ncol=length(coef(sp_mod)[-(1:(intercept+K*(2-linear)))]),
                                       dimnames= list(NULL,names(coef(sp_mod)[-(1:(intercept+K*(2-linear)))]))))
                     )
    pred_temp = SDMpredict(m=SIM$m_stan,focal=names(SIM$m_stan$model)[s],newdata = newdata,pred_samples = SIM$m_stan$iter, binary.resp = F,prob.cov = T)
    
    
    SIMlist[[i]]$pFund.mean.stan[,s] = apply(pred_temp$predictions.prob,1, mean) 
    SIMlist[[i]]$pFund.qinf.stan[,s] = apply(pred_temp$predictions.prob,1, quantile,0.025) 
    SIMlist[[i]]$pFund.qsup.stan[,s] = apply(pred_temp$predictions.prob,1, quantile,0.975) 
    
  }
  

}
  

### Compute CV fundamental & realized for both prob and bin
for(i in 1:length(SIMlist)){
  
  SIM=SIMlist[[i]]
  
  probCV = tSDM_CV_SIMUL(mod = SIM$m_stan, K = 5, fundNiche = T, prob.cov = T, iter = SIM$m_stan$iter,
                        pred_samples = pred_samples, error_prop_sample = error_prop_sample,
                        fitPreds=F,run.parallel=F, nEnv = nEnv)
  
  SIMlist[[i]]$pCV.mean.stan_prob = probCV$meanPred
  SIMlist[[i]]$pCV.qinf.stan_prob = probCV$Pred975
  SIMlist[[i]]$pCV.qsup.stan_prob = probCV$Pred025
  
  SIMlist[[i]]$pFundCV.mean.stan = probCV$fundNiche$FN.meanPred
  SIMlist[[i]]$pFundCV.qinf.stan = probCV$fundNiche$FN.meanPred
  SIMlist[[i]]$pFundCV.qsup.stan = probCV$fundNiche$FN.meanPred
  
  binCV = tSDM_CV_SIMUL(mod = SIM$m_stan, K=5, fundNiche=T, prob.cov=F,iter=SIM$m_stan$iter,
                         pred_samples = pred_samples, error_prop_sample=error_prop_sample, fitPreds=F,
                        run.parallel=F, nEnv = nEnv)
  
  SIMlist[[i]]$pCV.mean.stan_bin = binCV$meanPred
  SIMlist[[i]]$pCV.qinf.stan_bin = binCV$Pred975
  SIMlist[[i]]$pCV.qsup.stan_bin = binCV$Pred025
  
}

    

    
plot(SIMlist[[i]]$pCV.mean.stan_bin,  SIMlist[[i]]$pCV.mean.stan_prob)
    
    
    
nameOrder = cumsum(table(spNewNames)[spNewNames])[paste0("Y", 1:S)]  # table to alternate between name ordering

plotDistributions(SIMlist$GLV_abioticGR, plotprey=T, RN=T, main="Predicted and observed species distribution: GLV Growth Rates",filename=paste0(figPath,"RealizedNicheLV_GR.pdf"))

