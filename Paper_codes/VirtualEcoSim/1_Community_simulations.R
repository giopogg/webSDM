########################################################################################################################
########################################################################################################################
# Community simulation script from Giovanni Poggiato and Jérémy Andréoletti
# By running this script once, you simulate one 'replication' (i.e. 2550 communities) with
# the GLV-abioticGR, GLV-abioticKbasal and Ricker-abioticKbasal as described in
# Poggiato et al. 2023
# To obtain *exactly* the results of the paper, this script should be run 100 times. For each run, the variable 'job'
# should be specified with a number from 1 to 100.

# !!! Only the working directory and the parameter job should be specified to run the script

# Running the script will create a directory with all the simulated communities and model parameters inside
############################################################

# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with metanetwork.  If not, see <http://www.gnu.org/licenses/>


rm(list=ls())

#library("plot.matrix")
library("igraph")
library("gtools")
library("cheddar")
library("devtools") 
install_github("hallucigenia-sparsa/seqtime")  
library("seqtime") ###???
library("reshape2")
#library("gplots")
library("magrittr")
library("purrr")
library("readr")
library("matrixcalc")
library("vegan")
library("parallel")
source("1_Tools.R")
library("truncnorm")


# !!!! Here specify a number. Set to 1 for the first 'replication', to 2 for the second and so on...
job = 1

set.seed(as.numeric(gsub("job_", "", job))*100)  
# Parameters

# Choose the main simulation parameter values
S = 20       # number of species
L = 3         # number of trophic levels
nEnv = 51     # number of environmental samples
nRep = 50     # number of replicates
maxBI = 5     # maximum weight of biotic interaction
niche_breadthGR = 0.3
niche_breadth_Kbasal=0.05
niche_breadthVC=0.3 
# Name the directory in which to store the simulation results.
simPath = paste0("Simulations_S",S,"L",L,"_nEnv",nEnv,"_nRep",nRep,"_maxBI",maxBI,"_",job,"/")
dir.create(simPath, showWarnings = FALSE)




Stroph = round(2**(L:1)/sum(2**(L:1))*S)   # species trophic levels (predator:prey = 2)
names(Stroph) <- 1:L
trophL = rep(1:L, Stroph)                  # number of species in each trophic level

d = 1     # death rate due to intraspecific competition
p = 0.2   # interaction probability between predators and available preys (will be decreasing                    logarithmically with the trophic level)
spNames <- paste("Sp", 1:S, ".TL", trophL, sep='')  # species names


#################################################################################################
####### Build interaction matrix


PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
Stroph_cum <- cumsum(Stroph)
cdt_connected = cdt_preys = cdt_preds = FALSE
    
# Sample a random DAG (with constraints)
while(!cdt_connected|!cdt_preys|!cdt_preds){
    PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
        
        # Each higher trophic level predates all inferior trophic levels
    for (l in 2:L){
          PredMat[1:Stroph_cum[(l-1)],(Stroph_cum[(l-1)]+1):Stroph_cum[l]] <-
              sample(c(0,1), replace=TRUE, 
                     prob=c(1-p/log(exp(1)-2+l), p/log(exp(1)-2+l)),
                     size=Stroph_cum[l-1]*Stroph[l])
    }

        # Condition 1 : connected graph
    cdt_connected <- is.connected(graph_from_adjacency_matrix(PredMat))
        
    # Condition 2 : at least a prey in the previous trophic level for each predator
    cdt_preys <- all(sapply((Stroph[1]+1):S, function(s){
            l = trophL[s]  # trophic level
            preys = PredMat[ifelse(l>2,Stroph_cum[(l-2)]+1,1):Stroph_cum[(l-1)],s]  # preys of the trophic level just below
            return(sum(preys)>0)
        }))
    # Condition 3 : no predator is very weak, no prey is overpredated (risks of being always absent)
    # Choose interaction strenghts
    PredMat_w <- sapply(1:S,
                        function(j){sp <- PredMat[,j]
                                    sp[sp != 0] <- 
                                        rdirichlet(n=1,
                                                   alpha=rep(1,sum(sp!=0)))*ifelse(trophL[j]==L,2,2.5)*maxBI
                                    return(sp)})  # M2 : Dirichlet weights distribution
    colnames(PredMat_w) <- rownames(PredMat_w)
    cdt_preds <- all((colSums(PredMat_w)-colSums(t(PredMat_w)))[-(1:Stroph[1])] > maxBI) &
                     all((colSums(PredMat_w)-colSums(t(PredMat_w)))[1:Stroph[1]] > -3*maxBI)
        #cdt_preds <- TRUE
    
}

# Build graph from PredMat
G = graph_from_adjacency_matrix(t(PredMat))



IntMat <- PredMat_w
IntMat[lower.tri(IntMat)] <- -t(IntMat)[lower.tri(IntMat)]    # negative interactions for the preys
diag(IntMat)[1:Stroph[1]] <- -d    # negative density-dependance for basal species
    
write.csv2(IntMat, paste0(simPath,"InteractionMatrix.csv"))


##################################################################################################  
####### Build niche optima


envMin = 0; envMax = 1   # min-/maximal environmental abiotic values
envs = seq(envMin, envMax, length.out = nEnv)
    
# randomly sampled (but limits of the environmental gradient are avoided)
niche_optima = round(runif(S, envMin+(envMax-envMin)/20, envMax-(envMax-envMin)/20),2)
write.csv2(niche_optima, paste0(simPath,"niche_optima.csv"), row.names = FALSE)


##################################################################################################  
####### Run GLV abiotic GR
#### Run abiotic GLV GR
####

    
glv.finalStates.abioticGR <- try(
        mclapply(1:nRep, function(x){
        final_state=sapply(envs, function(e){
                    glv.out.abioticGR <- 
                        simGLV(Stroph, IntMat, spNames, 
                               reduceGR = TRUE, env=e, niche_breadth=niche_breadthGR, 
                               niche_optima=niche_optima)
                    final_state <- 
                        t(data.frame(PA=as.numeric(glv.out.abioticGR[,ncol(glv.out.abioticGR)]>0)))
                    
                    final_abundance <- t(data.frame(AB=glv.out.abioticGR[,ncol(glv.out.abioticGR)]))
                    return(list(PA=final_state,AB=final_abundance))
                    })
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1)
    )
    
# Clump to presence-absence
glv.finalStates.abioticGR_PA <- lapply(glv.finalStates.abioticGR, function(final_state){
                    # Set NA values to 0 (computation errors due to too small values)
                    final_state$final_state_PA[is.na(final_state$final_state_PA)] <- 0
                    final_state$final_state_PA <-
                        as.data.frame(final_state$final_state_PA, row.names = spNames)
                    colnames(final_state$final_state_PA) <- envs
                    return(final_state$final_state_PA)})

# Merge
glv.finalStates.abioticGR.merged <- lapply(1:length(glv.finalStates.abioticGR_PA),
                                           function(i) cbind(datasets = as.character(i),
                                                             rowN = rownames(glv.finalStates.abioticGR_PA[[i]]),
                                                             glv.finalStates.abioticGR_PA[[i]]))
    
glv.finalStates.abioticGR.merged <- do.call(rbind, glv.finalStates.abioticGR.merged)
    
# save
write.csv2(glv.finalStates.abioticGR.merged, paste0(simPath,"glv_finalStates_abioticGR.csv"))
    
#Abundance-> only needed to compute the fundamental niche
glv.finalAbundances.abioticGR <- lapply(glv.finalStates.abioticGR, function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_AB[is.na(final_state$final_state_AB)] <- 0
        final_state$final_state_AB <- as.data.frame(t(final_state$final_state_AB), row.names = envs)
        colnames(final_state$final_state_AB) <- spNames
        return(final_state$final_state_AB)})
    
glv.meanAbundances.abioticGR <- sapply(Reduce(rbind, glv.finalAbundances.abioticGR),
                                       function(x)weighted.mean(x, w=x>0))
    
write.csv2(glv.meanAbundances.abioticGR, paste0(simPath,"glv.meanAbundances.abioticGR.csv"))
    



##################################################################################################  
#### Run abiotic GLV abiotic K basal
####

glv.finalStates.abioticKbasal <- mclapply(1:nRep, function(x){
        final_state=sapply(envs, function(e){
            glv.out.abioticKbasal <- simGLV(Stroph, IntMat, spNames, reduceK = TRUE,
                                            env=e, niche_breadth = niche_breadth_Kbasal,
                                            niche_optima=niche_optima)
            final_state <- t(data.frame(PA=as.numeric(glv.out.abioticKbasal[,ncol(glv.out.abioticKbasal)]>0)))
            final_abundance <- t(data.frame(AB=glv.out.abioticKbasal[,ncol(glv.out.abioticKbasal)]))
            return(list(PA=final_state,AB=final_abundance))
            })
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1)
    
# Clump to presence-absence
glv.finalStates.abioticKbasal_PA <- lapply(glv.finalStates.abioticKbasal,
                                           function(final_state){
        # Set NA values to 0 (computation errors due to too small values)
        final_state$final_state_PA[is.na(final_state$final_state_PA)] <- 0
        final_state$final_state_PA <- as.data.frame(final_state$final_state_PA, row.names = spNames)
        colnames(final_state$final_state_PA) <- envs
        return(final_state$final_state_PA)})

# Merge
glv.finalStates.abioticKbasal.merged <- lapply(1:length(glv.finalStates.abioticKbasal_PA),
                                               function(i) cbind(datasets = as.character(i),
                                                                 rowN = rownames(glv.finalStates.abioticKbasal_PA[[i]]), 
                                                                 glv.finalStates.abioticKbasal_PA[[i]]))

glv.finalStates.abioticKbasal.merged <- do.call(rbind, glv.finalStates.abioticKbasal.merged)

# Save
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


##################################################################################################  
####### Run Ricker  abiotic K basal
#######

sigma_ricker = 0.2


ricker.finalStates.abioticKbasal <- try( mclapply(1:nRep, function(x){
        final_state=sapply(envs,
                           function(e){ 
                               ricker.out.abioticKbasal <- simRicker(Stroph, t(IntMat),
                                                                     spNames, reduceK = TRUE,
                                                                     env=e, niche_optima=niche_optima,
                                                                     sigma=sigma_ricker,
                                                                     death.t=10^-15, tend=1000)
                                final_state <- 
                                    t(data.frame(PA=as.numeric( ricker.out.abioticKbasal[,ncol(ricker.out.abioticKbasal)]>0)))
                                
                                final_abundance <- t(data.frame(AB=ricker.out.abioticKbasal[,ncol(ricker.out.abioticKbasal)]))
                                return(list(PA=final_state,AB=final_abundance))
                                })
        
        final_state_PA=apply(final_state, MARGIN = 2, FUN = function(x){x[[1]]} )
        final_state_AB=apply(final_state, MARGIN = 2, FUN = function(x){x[[2]]} )
        
        return(list(final_state_PA=final_state_PA,final_state_AB=final_state_AB))
    },mc.cores = detectCores()-1), TRUE)
    
# Only if simulations did not explode !
if(!any(unlist(lapply(ricker.finalStates.abioticKbasal, function(x) class(x)=="try-error")))){
    ricker.finalStates.abioticKbasal_PA <- lapply(ricker.finalStates.abioticKbasal, function(final_state){
        final_state$final_state_PA <- as.data.frame(final_state$final_state_PA, row.names = spNames)
        colnames(final_state$final_state_PA) <- envs
        return(final_state$final_state_PA)})
    
    
ricker.finalStates.abioticKbasal.merged <- lapply(1:length(ricker.finalStates.abioticKbasal_PA),
                                                  function(i)  cbind(datasets = as.character(i), 
                                                                     rowN = rownames(ricker.finalStates.abioticKbasal_PA[[i]]),
                                                                     ricker.finalStates.abioticKbasal_PA[[i]]))

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
if(!any(unlist(lapply(glv.finalStates.abioticGR, function(x) class(x)=="try-error")))){
    
glv.fundNicheTh.abioticGR = t(sapply(envs, function(e)1-ptruncnorm(1-exp(-(e-niche_optima)**2/(2*niche_breadthGR**2)) - glv.meanAbundances.abioticGR %*% (IntMat*PredMat), 
                                                                   mean=Reduce('+', b.fundNiche), sd=Reduce('+', sd_noise.fundNiche), 
                                                                   a=ifelse(trophL==1,0,-Inf), b=ifelse(trophL>1,0,Inf))))
colnames(glv.fundNicheTh.abioticGR) = spNames
write.csv2(glv.fundNicheTh.abioticGR, paste0(simPath,"glv.fundNicheTh.abioticGR.csv"))
}
###################################################################################################################
#Everything again for GLV K  

niche_breadth_Kbasal=niche_breadth_Kbasal

if(!any(unlist(lapply(glv.finalStates.abioticKbasal, function(x) class(x)=="try-error")))){

glv.fundNicheTh.abioticKbasal = t(sapply(envs, function(e)1-ptruncnorm(1e-10*exp((e-niche_optima)**2/(2*niche_breadth_Kbasal**2))*(trophL==1) - glv.meanAbundances.abioticKbasal %*% (IntMat*PredMat), 
                                                                       mean = Reduce('+', b.fundNiche[-1])+(trophL==1), sd=Reduce('+', sd_noise.fundNiche[-1]), 
                                                                       a = ifelse(trophL==1,-Inf,-Inf), b = ifelse(trophL>1,0,Inf))))
#glv.fundNicheTh.abioticKbasal[,trophL==1] <- t(sapply(envs, function(e) -diag(IntMat.fundNiche$TL1)*exp(-(e-niche_optima)**2/(2*niche_breadth**2))) > 1e-10)[,trophL==1]
colnames(glv.fundNicheTh.abioticKbasal) = spNames
write.csv2(glv.fundNicheTh.abioticKbasal, paste0(simPath,"glv.fundNicheTh.abioticKbasal.csv"))
}
#################################################################################################
### Ricker

niche_breadth = 0.05 #defined inside the simRicker
sigma = sigma_ricker
mu = 0.02   #defined inside the simRicker


hitting_proba <- function(Y0, drift, sigma, a, b){
    return((1-exp(2*drift/sigma**2*(b-Y0)))/(1-exp(2*drift/sigma**2*(b-a))))
}

if(!any(unlist(lapply(ricker.finalStates.abioticKbasal, function(x) class(x)=="try-error")))){
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

}
