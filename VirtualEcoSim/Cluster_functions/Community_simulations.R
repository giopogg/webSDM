rm(list=ls())
.libPaths( c( .libPaths(), "~/my_r_libraries/") )   

# Choose wether to simulate or simply load already simulated datasets
loadSim = FALSE
setwd("~/Documents/GitHub/trophicSDM/VirtualEcoSim")

args= commandArgs(trailingOnly = TRUE)
job=args[1]
# Parameters
if (!loadSim){
    # Choose the main simulation parameter values
    S = 5       # number of species
    L = 3         # number of trophic levels
    nEnv = 51     # number of environmental samples
    nRep = 50     # number of replicates
    maxBI = 5     # maximum weight of biotic interaction
    niche_breadthGR = 0.3
    niche_breadth_Kbasal=0.05
    niche_breadthVC=0.3 
    # Name the directory in which to store the simulations
    simPath = paste0("Simulations_S",S,"L",L,"_nEnv",nEnv,"_nRep",nRep,"_maxBI",maxBI,"_",job,"/")
    dir.create(simPath, showWarnings = FALSE)
}else{
    # Choose the main simulation parameter values
    S = 5       # number of species
    L = 3         # number of trophic levels
    nEnv = 51     # number of environmental samples
    nRep = 50     # number of replicates
    maxBI = 5     # maximum weight of biotic interaction
    niche_breadthGR = 0.3
    niche_breadth_Kbasal=0.05
    niche_breadthVC=0.3 
    # Name the directory in which to store the simulations
    simPath = paste0("Simulations_S",S,"L",L,"_nEnv",nEnv,"_nRep",nRep,"_maxBI",maxBI,"_job",job,"/")
    dir.create(simPath, showWarnings = FALSE)
    # Choose the directory from which to load the simulations
    simPath = paste0("Simulations_S",S,"L",L,"_nEnv",nEnv,"_nRep",nRep,"_maxBI",maxBI,"_job",job,"/")

}


library("igraph")
library("gtools")
library("cheddar")
library("devtools") 
library("seqtime") 
library("reshape2")
library("magrittr")
library("purrr")
library("readr")
library("matrixcalc")
library("vegan")
library("parallel")
source("Sim_functions.R")
library("truncnorm")

if (!loadSim){
    # Stochastic version
    #trophL = NA
    #while(length(unique(trophL))!=L) trophL = sample(1:L, size = S, replace = TRUE, prob = 2**(L:1))  # Species trophic levels (predator:prey = 2)
    #Stroph = table(trophL)

    # Deterministic version
    Stroph = round(2**(L:1)/sum(2**(L:1))*S)   # species trophic levels (predator:prey = 2)
    names(Stroph) <- 1:L
    trophL = rep(1:L, Stroph)                  # number of species in each trophic level

    d = 1     # death rate due to intraspecific competition
    p = 0.2   # interaction probability between predators and available preys (will be decreasing logarithmically with the trophic level)
    spNames <- paste("Sp", 1:S, ".TL", trophL, sep='')  # species names
}

##################################################################################################  
####### Build interaction matrix

if (loadSim){
    IntMat <- as.matrix(read.csv2(paste0(simPath,"InteractionMatrix.csv"), row.names = 1))
    PredMat_w <- IntMat
    PredMat_w[lower.tri(PredMat_w, diag=TRUE)] <- 0
    PredMat <- sign(PredMat_w)
    S = ncol(PredMat)
    spNames = colnames(PredMat)
    trophL = as.numeric(sub(pattern = "Sp.*TL", replacement = "", x = spNames))
    Stroph = table(trophL)
    Stroph_cum <- cumsum(Stroph)
    L = max(trophL)
    
}else{
    PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
    Stroph_cum <- cumsum(Stroph)
    cdt_connected = cdt_preys = cdt_preds = FALSE
    
    # Random DAG
    while(!cdt_connected|!cdt_preys|!cdt_preds){
        PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
        
        # Each higher trophic level predates all inferior trophic levels
        for (l in 2:L){
          PredMat[1:Stroph_cum[(l-1)],(Stroph_cum[(l-1)]+1):Stroph_cum[l]] <- sample(c(0,1), replace=TRUE, prob=c(1-p/log(exp(1)-2+l), p/log(exp(1)-2+l)), size=Stroph_cum[l-1]*Stroph[l])
        }

        # Condition 1 : connected graph
        cdt_connected <- is.connected(graph_from_adjacency_matrix(PredMat))
        
        # Condition 2 : at least a prey in the previous trophic level for each predator
        cdt_preys <- all(sapply((Stroph[1]+1):S, function(s){
            l = trophL[s]  # trophic level
            preys = PredMat[ifelse(l>2,Stroph_cum[(l-2)]+1,1):Stroph_cum[(l-1)],s]  # preys of the trophic level just below
            return(sum(preys)>0)
        }))
        #cdt_preys <- all(colSums(PredMat)[-(1:Stroph[1])]>0)  # Analogous condition for all preys
        #cdt_preys <- TRUE
                
        # Condition 3 : no predator is very weak, no prey is overpredated (risks of being always absent)
         # Choose interaction strenghts
        #PredMat_w <- PredMat*runif(S**2,0,maxBI)  # M1 : uniform weights distribution
        PredMat_w <- sapply(1:S, function(j){sp <- PredMat[,j]
                                             sp[sp != 0] <- rdirichlet(n=1, alpha=rep(1,sum(sp!=0)))*ifelse(trophL[j]==L,2,2.5)*maxBI
                                             return(sp)})  # M2 : Dirichlet weights distribution
        colnames(PredMat_w) <- rownames(PredMat_w)
        cdt_preds <- all((colSums(PredMat_w)-colSums(t(PredMat_w)))[-(1:Stroph[1])] > maxBI) &
                     all((colSums(PredMat_w)-colSums(t(PredMat_w)))[1:Stroph[1]] > -3*maxBI)
        #cdt_preds <- TRUE
    }
}

# Build graph from PredMat
G = graph_from_adjacency_matrix(t(PredMat))


if (!loadSim){
    IntMat <- PredMat_w
    IntMat[lower.tri(IntMat)] <- -t(IntMat)[lower.tri(IntMat)]    # negative interactions for the preys
    diag(IntMat)[1:Stroph[1]] <- -d    # negative density-dependance for basal species
    
    write.csv2(IntMat, paste0(simPath,"InteractionMatrix.csv"))
}

##################################################################################################  
####### Build niche optima

if (loadSim){
    niche_optima <- read.csv2(paste0(simPath,"niche_optima.csv"))[,1]
    if (S!=length(niche_optima)) print("WARNING : Wrong number of species")
}else{
    envMin = 0; envMax = 1   # min-/maximal environmental abiotic values
    envs = seq(envMin, envMax, length.out = nEnv)
    
    niche_optima = round(runif(S, envMin+(envMax-envMin)/20, envMax-(envMax-envMin)/20),2)
    write.csv2(niche_optima, paste0(simPath,"niche_optima.csv"), row.names = FALSE)
}

##################################################################################################  
####### Run GLV abiotic GR
#### Run abiotic GLV GR
####
if (loadSim){
    glv.finalStates.abioticGR <- loadData(paste0(simPath,"glv_finalStates_abioticGR.csv"))
    nRep = length(glv.finalStates.abioticGR)
    envs = as.numeric(names(glv.finalStates.abioticGR[[1]]))
    nEnv = length(envs)
    
    glv.meanAbundances.abioticGR <- read.csv2(paste0(simPath,"glv.meanAbundances.abioticGR.csv"), row.names=1)[[1]]
    names(glv.meanAbundances.abioticGR) <- spNames
}else{
    
    glv.finalStates.abioticGR <- mclapply(1:nRep, function(x){
        final_state=sapply(envs, function(e){
                    glv.out.abioticGR <- simGLV(Stroph, IntMat, spNames, reduceGR = TRUE, env=e, niche_breadth=niche_breadthGR, niche_optima=niche_optima)
                    final_state <- t(data.frame(PA=as.numeric(glv.out.abioticGR[,ncol(glv.out.abioticGR)]>0)))
                    final_abundance <- t(data.frame(AB=glv.out.abioticGR[,ncol(glv.out.abioticGR)]))
                    return(list(PA=final_state,AB=final_abundance))})
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1)
    
    # Presence-absence
    glv.finalStates.abioticGR_PA <- lapply(glv.finalStates.abioticGR, function(final_state){
                    # Set NA values to 0 (computation errors due to too small values)
                    final_state$final_state_PA[is.na(final_state$final_state_PA)] <- 0
                    final_state$final_state_PA <- as.data.frame(final_state$final_state_PA, row.names = spNames)
                    colnames(final_state$final_state_PA) <- envs
                    return(final_state$final_state_PA)})
    
    glv.finalStates.abioticGR.merged <- lapply(1:length(glv.finalStates.abioticGR_PA), function(i)
        cbind(datasets = as.character(i), rowN = rownames(glv.finalStates.abioticGR_PA[[i]]), glv.finalStates.abioticGR_PA[[i]]))
    
    glv.finalStates.abioticGR.merged <- do.call(rbind, glv.finalStates.abioticGR.merged)
    
    write.csv2(glv.finalStates.abioticGR.merged, paste0(simPath,"glv_finalStates_abioticGR.csv"))
    
    #Abundance-> only needed to compute the fundamental niche
    glv.finalAbundances.abioticGR <- lapply(glv.finalStates.abioticGR, function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_AB[is.na(final_state$final_state_AB)] <- 0
        final_state$final_state_AB <- as.data.frame(t(final_state$final_state_AB), row.names = envs)
        colnames(final_state$final_state_AB) <- spNames
        return(final_state$final_state_AB)})
    
    glv.meanAbundances.abioticGR <- sapply(Reduce(rbind, glv.finalAbundances.abioticGR), function(x)weighted.mean(x, w=x>0))
    
    write.csv2(glv.meanAbundances.abioticGR, paste0(simPath,"glv.meanAbundances.abioticGR.csv"))
    
}




##################################################################################################  
#### Run abiotic GLV abiotic K basal
####

if (loadSim){
    glv.finalStates.abioticKbasal <- loadData(paste0(simPath,"glv_finalStates_abioticKbasal.csv"))
    nRep = length(glv.finalStates.abioticKbasal)
    
    glv.meanAbundances.abioticKbasal <- read.csv2(paste0(simPath,"glv.meanAbundances.abioticKbasal.csv"), row.names=1)[[1]]
    names(glv.meanAbundances.abioticKbasal) <- spNames
    
}else{

    glv.finalStates.abioticKbasal <- mclapply(1:nRep, function(x){
        final_state=sapply(envs, function(e){
            glv.out.abioticKbasal <- simGLV(Stroph, IntMat, spNames, reduceK = TRUE, env=e, niche_breadth = niche_breadth_Kbasal, niche_optima=niche_optima)
            final_state <- t(data.frame(PA=as.numeric(glv.out.abioticKbasal[,ncol(glv.out.abioticKbasal)]>0)))
            final_abundance <- t(data.frame(AB=glv.out.abioticKbasal[,ncol(glv.out.abioticKbasal)]))
            return(list(PA=final_state,AB=final_abundance))})
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1)
    
    # Presence-absence
    glv.finalStates.abioticKbasal_PA <- lapply(glv.finalStates.abioticKbasal, function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_PA[is.na(final_state$final_state_PA)] <- 0
        final_state$final_state_PA <- as.data.frame(final_state$final_state_PA, row.names = spNames)
        colnames(final_state$final_state_PA) <- envs
        return(final_state$final_state_PA)})
    
    glv.finalStates.abioticKbasal.merged <- lapply(1:length(glv.finalStates.abioticKbasal_PA), function(i)
                cbind(datasets = as.character(i), rowN = rownames(glv.finalStates.abioticKbasal_PA[[i]]), glv.finalStates.abioticKbasal_PA[[i]]))
    glv.finalStates.abioticKbasal.merged <- do.call(rbind, glv.finalStates.abioticKbasal.merged)
    write.csv2(glv.finalStates.abioticKbasal.merged, paste0(simPath,"glv_finalStates_abioticKbasal.csv"))
    
    
    #Abundance-> only needed to compute the fundamental niche
    glv.finalAbundances.abioticKbasal <- lapply(glv.finalStates.abioticKbasal, function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_AB[is.na(final_state$final_state_AB)] <- 0
        final_state$final_state_AB <- as.data.frame(t(final_state$final_state_AB), row.names = envs)
        colnames(final_state$final_state_AB) <- spNames
        return(final_state$final_state_AB)})
    
    glv.meanAbundances.abioticKbasal <- sapply(Reduce(rbind, glv.finalAbundances.abioticKbasal), function(x)weighted.mean(x, w=x>0))
    
    write.csv2(glv.meanAbundances.abioticKbasal, paste0(simPath,"glv.meanAbundances.abioticKbasal.csv"))
    
    
}


##################################################################################################  
####### Run Ricker  abiotic K basal
#######

sigma_ricker = 0.2

if(loadSim){
    ricker.finalStates.abioticKbasal <- loadData(paste0(simPath,"ricker_finalStates_abioticKbasal.csv"))
    nRep = length(ricker.finalStates.abioticKbasal)
    
    ricker.meanAbundances.abioticKbasal <- read.csv2(paste0(simPath,"ricker.meanAbundances.abioticKbasal.csv"), row.names=1)[[1]]
    names(ricker.meanAbundances.abioticKbasal) <- spNames
    
}else{

    ricker.finalStates.abioticKbasal <- mclapply(1:nRep, function(x){
        final_state=sapply(envs, function(e){
            ricker.out.abioticKbasal <- simRicker(Stroph, t(IntMat), spNames, reduceK = TRUE, env=e, niche_optima=niche_optima, sigma=sigma_ricker, death.t=10^-15, tend=1000)
            final_state <- t(data.frame(PA=as.numeric( ricker.out.abioticKbasal[,ncol(ricker.out.abioticKbasal)]>0)))
            final_abundance <- t(data.frame(AB=ricker.out.abioticKbasal[,ncol(ricker.out.abioticKbasal)]))
            return(list(PA=final_state,AB=final_abundance))})
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1)
    
    
    ricker.finalStates.abioticKbasal_PA <- lapply(ricker.finalStates.abioticKbasal, function(final_state){
        final_state$final_state_PA <- as.data.frame(final_state$final_state_PA, row.names = spNames)
        colnames(final_state$final_state_PA) <- envs
        return(final_state$final_state_PA)})
    
    
    ricker.finalStates.abioticKbasal.merged <- lapply(1:length(ricker.finalStates.abioticKbasal_PA), function(i)
               cbind(datasets = as.character(i), rowN = rownames(ricker.finalStates.abioticKbasal_PA[[i]]), ricker.finalStates.abioticKbasal_PA[[i]]))
    ricker.finalStates.abioticKbasal.merged <- do.call(rbind, ricker.finalStates.abioticKbasal.merged)
    write.csv2(ricker.finalStates.abioticKbasal.merged, paste0(simPath,"ricker_finalStates_abioticKbasal.csv"))
    
    #Abundance-> only needed to compute the fundamental niche
    ricker.finalAbundances.abioticKbasal <- lapply(ricker.finalStates.abioticKbasal, function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_AB[is.na(final_state$final_state_AB)] <- 0
        final_state$final_state_AB <- as.data.frame(t(final_state$final_state_AB), row.names = envs)
        colnames(final_state$final_state_AB) <- spNames
        return(final_state$final_state_AB)})
    
    ricker.meanAbundances.abioticKbasal <- sapply(Reduce(rbind, ricker.finalAbundances.abioticKbasal), function(x)weighted.mean(x, w=x>0))
    
    write.csv2(ricker.meanAbundances.abioticKbasal, paste0(simPath,"ricker.meanAbundances.abioticKbasal.csv"))
    
}

##################################################################################################  
####### Run SOI abiotic ER
#######

S = nrow(IntMat)  # number of species
I = 200           # total number of individuals
tend = 50         # ending time
m.vector = c(rep(1,Stroph[1]), rep(0.5,S-Stroph[1]))  # migration rates, lower for predator species
e.vector = c(rep(1,Stroph[1]), rep(1.5,S-Stroph[1]))  # extinction rates, higher for predator species

if (loadSim){
    soi.finalStates.abioticER <- loadData(paste0(simPath,"soi_finalStates_abioticER.csv"))
    nRep = length(soi.finalStates.abioticER)
}else{
    soi.finalStates.abioticER <- mclapply(1:nRep, function(x)sapply(envs, function(e){
                soi.out.abiotic <- simSOI(Stroph, I, IntMat, spNames, m.vector=m.vector, e.vector=e.vector, increaseE=TRUE, env=e, niche_optima=niche_optima, tend=tend)
                final_state <- t(data.frame(PA=as.numeric(rowSums(soi.out.abiotic[,(tend-10):tend])>I/10), row.names = spNames))
                return(final_state)}),mc.cores = detectCores()-1)
                                    
   soi.finalStates.abioticER <- lapply(soi.finalStates.abioticER, function(final_state){
                final_state <- as.data.frame(final_state, row.names = spNames)
                colnames(final_state) <- envs
                return(final_state)})
                                        
    soi.finalStates.abioticER.merged <- lapply(1:length(soi.finalStates.abioticER), function(i)
                cbind(datasets = as.character(i), rowN = rownames(soi.finalStates.abioticER[[i]]), soi.finalStates.abioticER[[i]]))
    soi.finalStates.abioticER.merged <- do.call(rbind, soi.finalStates.abioticER.merged)
    write.csv2(soi.finalStates.abioticER.merged, paste0(simPath,"soi_finalStates_abioticER.csv"))
}


##################################################################################################  
####### Run SOI abiotic ER basal
#######

S = nrow(IntMat)  # number of species
I = 200           # total number of individuals
tend = 50         # ending time
m.vector = c(rep(1,Stroph[1]), rep(0.5,S-Stroph[1]))  # migration rates, lower for predator species
e.vector = c(rep(1,Stroph[1]), rep(1.5,S-Stroph[1]))  # extinction rates, higher for predator species

if (loadSim){
    soi.finalStates.abioticERbasal <- loadData(paste0(simPath,"soi_finalStates_abioticERbasal.csv"))
    nRep = length(soi.finalStates.abioticERbasal)
}else{
    soi.finalStates.abioticERbasal <- mclapply(1:nRep, function(x)sapply(envs, function(e){
                soi.out.abioticERbasal <- simSOI(Stroph, I, IntMat, spNames, m.vector=m.vector, e.vector=e.vector, increaseE=TRUE, basal=TRUE, env=e, niche_optima=niche_optima, tend=tend)
                final_state <- t(data.frame(PA=as.numeric(rowSums(soi.out.abioticERbasal[,(tend-10):tend])>I/10), row.names = spNames))
                return(final_state)}),mc.cores = detectCores()-1)
                                         
    soi.finalStates.abioticERbasal <- lapply(soi.finalStates.abioticERbasal, function(final_state){
                final_state <- as.data.frame(final_state, row.names = spNames)
                colnames(final_state) <- envs
                return(final_state)})
    
    soi.finalStates.abioticERbasal.merged <- lapply(1:length(soi.finalStates.abioticERbasal), function(i)
                cbind(datasets = as.character(i), rowN = rownames(soi.finalStates.abioticERbasal[[i]]), soi.finalStates.abioticERbasal[[i]]))
    soi.finalStates.abioticERbasal.merged <- do.call(rbind, soi.finalStates.abioticERbasal.merged)
    write.csv2(soi.finalStates.abioticERbasal.merged, paste0(simPath,"soi_finalStates_abioticERbasal.csv"))
}

##################################################################################################  
####### Run VirtualCom
#######

fac_inter <- PredMat_w

comp_inter <- t(PredMat_w)
diag(comp_inter)[1:Stroph[1]] <- 1

env = envs                        # environment value at each site
niche_optima  = niche_optima      # niche means
niche_breadthVC = niche_breadthVC               # niche variance
comp_inter = comp_inter           # competition
fac_inter = fac_inter             # facilitation
beta_env = 1 ; beta_comp = 1 ; beta_fac = 1 ; beta_abun = 1  # filters weights
years = 20                        # years number of simulated time-steps
K = S                             # number of individuals in the community

if (loadSim){
    vc.finalStates.abiotic <- loadData(paste0(simPath,"vc_finalStates_abiotic.csv"))
    nRep = length(vc.finalStates.abiotic)
}else{
    vc.finalStates.abiotic <- lapply(1:nRep, function(x){
                                    final_state <- t(simulate_community(env, niche_optima, niche_breadthVC, comp_inter, fac_inter, 
                                       beta_env, beta_comp, beta_fac, beta_abun, years, K)[-(S+1)])
                                     colnames(final_state) <- envs
                                     return (final_state)})
    
    vc.finalStates.abiotic.merged <- lapply(1:length(vc.finalStates.abiotic), function(i)
                cbind(datasets = as.character(i), rowN = rownames(vc.finalStates.abiotic[[i]]), vc.finalStates.abiotic[[i]]))
    vc.finalStates.abiotic.merged <- do.call(rbind, vc.finalStates.abiotic.merged)
    write.csv2(vc.finalStates.abiotic.merged, paste0(simPath,"vc_finalStates_abiotic.csv"))
}


###################################################################################################################
####################### Fundamental niche
############################################################################################

###################################################################################################################
# GLV GR

# Intrinsic growth rates (null for preys). the density dep growth rate is instead negative. This is consistent with above simulations
b.fundNiche = lapply(1:L, function(l)rep(c(0.5,-0.2), times=c(Stroph_cum[max(1,l-1)], S-Stroph_cum[max(1,l-1)]))*(trophL==l))
names(b.fundNiche) <- paste0("TL",1:L)

# Strength of the noise added to the growth rates (must match the one defined above)
sd_noise.fundNiche = lapply(1:L, function(l) 0.1*(trophL == l))
names(sd_noise.fundNiche) <- paste0("TL",1:L)

niche_breadthGR = niche_breadthGR

# Compute the theoretical niche as a function of the above defined parameters and the mean abundances computed above.
# see notebook and rapport de stage for derivation of formulas

glv.fundNicheTh.abioticGR = t(sapply(envs, function(e)1-ptruncnorm(1-exp(-(e-niche_optima)**2/(2*niche_breadthGR**2)) - glv.meanAbundances.abioticGR %*% (IntMat*PredMat), 
                                                                   mean=Reduce('+', b.fundNiche), sd=Reduce('+', sd_noise.fundNiche), 
                                                                   a=ifelse(trophL==1,0,-Inf), b=ifelse(trophL>1,0,Inf))))
colnames(glv.fundNicheTh.abioticGR) = spNames
write.csv2(glv.fundNicheTh.abioticGR, paste0(simPath,"glv.fundNicheTh.abioticGR.csv"))

###################################################################################################################
#Everything again for GLV K  

niche_breadth_Kbasal=niche_breadth_Kbasal

glv.fundNicheTh.abioticKbasal = t(sapply(envs, function(e)1-ptruncnorm(1e-10*exp((e-niche_optima)**2/(2*niche_breadth_Kbasal**2))*(trophL==1) - glv.meanAbundances.abioticKbasal %*% (IntMat*PredMat), 
                                                                       mean = Reduce('+', b.fundNiche[-1])+(trophL==1), sd=Reduce('+', sd_noise.fundNiche[-1]), 
                                                                       a = ifelse(trophL==1,-Inf,-Inf), b = ifelse(trophL>1,0,Inf))))
#glv.fundNicheTh.abioticKbasal[,trophL==1] <- t(sapply(envs, function(e) -diag(IntMat.fundNiche$TL1)*exp(-(e-niche_optima)**2/(2*niche_breadth**2))) > 1e-10)[,trophL==1]
colnames(glv.fundNicheTh.abioticKbasal) = spNames
write.csv2(glv.fundNicheTh.abioticKbasal, paste0(simPath,"glv.fundNicheTh.abioticKbasal.csv"))


#################################################################################################
### Ricker

niche_breadth = 0.05 #defined inside the simRicker
sigma = sigma_ricker
mu = 0.02   #defined inside the simRicker


hitting_proba <- function(Y0, drift, sigma, a, b){
    return((1-exp(2*drift/sigma**2*(b-Y0)))/(1-exp(2*drift/sigma**2*(b-a))))
}

ricker.fundNicheTh.abioticKbasal = t(sapply(envs, function(e){
    K = rep(0.1,S)*pmax(exp(-(niche_optima-e)**2/(2*niche_breadth**2)), 1e-15)
    fundNicheTh =  (trophL==1)*(1-pnorm(1e-15, mean=t(-K*diag(IntMat)), sd=sigma/1000)/3) + # inadequate function, manual fitting
        (trophL==2)*t((1-hitting_proba(Y0=0, drift=t(ricker.meanAbundances.abioticKbasal %*% (IntMat * ((trophL<=2) %*% t((trophL<=2))))) - K*diag(IntMat) - mu*(1+diag(IntMat)),
                                       sigma=sigma, a=-15, b=8))) +
        (trophL==3)*t((1-hitting_proba(Y0=0, drift=t(ricker.meanAbundances.abioticKbasal %*% IntMat) - K*diag(IntMat) - mu*(1+diag(IntMat)),
                                       sigma=sigma, a=-15, b=8)))
    return(fundNicheTh)
}))
colnames(ricker.fundNicheTh.abioticKbasal) = spNames
write.csv2(ricker.fundNicheTh.abioticKbasal, paste0(simPath,"ricker.fundNicheTh.abioticKbasal.csv"))


##########################################################################################
# Additional plots

# 
# options(repr.plot.width=12, repr.plot.height=7, repr.plot.res=150)
# par(mfrow=c(2,3), mar=c(4,4,3,1))
# 
# glv.out.abioticGR <- simGLV(Stroph, IntMat, spNames, reduceGR = TRUE, env=0.6, niche_optima=niche_optima, niche_breadth=0.3, tend=200)
# groups <- glv.out.abioticGR[1,] ; suppressWarnings(groups[1:length(groups)] <- trophL)
# glv.out.abioticGR[glv.out.abioticGR>1] <- NA
# tsplot(glv.out.abioticGR, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="GLV - Abiotic effect on Growth Rates", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# glv.out.abioticKbasal <- simGLV(Stroph, IntMat, spNames, reduceK = TRUE, env=envs[1], niche_optima=niche_optima, tend = 200)
# groups <- glv.out.abioticKbasal[1,] ; suppressWarnings(groups[1:length(groups)] <- trophL)
# glv.out.abioticKbasal[glv.out.abioticKbasal>0.1] <- NA
# tsplot(glv.out.abioticKbasal, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="GLV - Abiotic effect on Carrying Capacities", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# ricker.out.abioticKbasal <- simRicker(Stroph, t(IntMat), spNames, niche_optima=niche_optima, reduceK=TRUE, env=0.6, sigma=0.1, death.t=10^-25, tend=2000)
# groups <- ricker.out.abioticKbasal[1,] ; groups[1:length(groups)] <- trophL
# ricker.out.abioticKbasal[ricker.out.abioticKbasal>0.6] <- NA
# tsplot(ricker.out.abioticKbasal, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="Ricker - Abiotic effect on Carrying Capacities", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# soi.out.abiotic <- simSOI(Stroph, I, IntMat, spNames, m.vector=m.vector, e.vector=e.vector, increaseE=TRUE, 
#                           env=envs[5], niche_optima=niche_optima, tend=100)
# groups <- soi.out.abiotic[1,] ; suppressWarnings(groups[1:length(groups)] <- trophL)
# tsplot(soi.out.abiotic, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="SOI - Abiotic effect on all Extinction Rates", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# soi.out.abioticERbasal <- simSOI(Stroph, I, IntMat, spNames, m.vector=m.vector, e.vector=e.vector, increaseE=TRUE, 
#                                  basal=TRUE, env=envs[1], niche_optima=niche_optima, tend=tend)
# groups <- soi.out.abioticERbasal[1,] ; suppressWarnings(groups[1:length(groups)] <- trophL)
# tsplot(soi.out.abioticERbasal, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="SOI - Abiotic effect on Basal Species' Extinction Rates", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# vc.out.abioticERbasal <- sapply(1:years, function(year)colSums(simulate_community(rep(0.6,K/2), niche_optima, niche_breadth=0.3, comp_inter, fac_inter, 
#                                                                                   beta_env, beta_comp, beta_fac, beta_abun, year, K/2)))
# groups <- vc.out.abioticERbasal[1,] ; suppressWarnings(groups[1:length(groups)] <- trophL)
# tsplot(vc.out.abioticERbasal, groups=groups, my.color.map = hsv(h = seq(0.3,0,length.out=L+1)), main="VirtualCom - Abiotic effect on Recruitment Probability", legend=F)
# legend("topright", lty=1, col=hsv(h = seq(0.3,0,length.out=L+1)), paste("TL", 1:L, sep=""), cex=1.2, y.intersp=2)
# 
# invisible(dev.copy(svg, file=paste0("Images/Species_dynamics_all_simulation_models.svg"), width=12, height=7))
# invisible(dev.off())
# 
# 
