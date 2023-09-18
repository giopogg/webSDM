### Functions to simulate communities
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


#############################################################################################
##### Global workflow for sensitivity analysis on Lotka-Volterra simulations
#############################################################################################

runGRSimulation = function(S = S, L = L, p = p, niche_breadthGR = niche_breadthGR, nEnv = nEnv,
                           nRep, simPath = simPath){
  
  
  maxBI = 5  # maximum weight of biotic interaction
  
  Stroph = round(2**(L:1)/sum(2**(L:1))*S)   # species trophic levels (predator:prey = 2)
  names(Stroph) <- 1:L
  trophL = rep(1:L, Stroph)                  # number of species in each trophic level
  d = 1     # death rate due to intraspecific competition
  spNames <- paste("Sp", 1:S, ".TL", trophL, sep='')  # species names
  
  
  ##################################################################################################  
  ####### Build interaction matrix
  
  
  PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
  Stroph_cum <- cumsum(Stroph)
  cdt_connected = cdt_preys = cdt_preds = FALSE
  
  # Random DAG
  while(!cdt_connected|!cdt_preys|!cdt_preds){
    PredMat <- matrix(data = 0, nrow=S, ncol=S, dimnames=list(spNames, spNames))
    
    # Each higher trophic level predates all inferior trophic levels
    for (l in 2:L){
      PredMat[1:Stroph_cum[(l-1)],
              (Stroph_cum[(l-1)]+1):Stroph_cum[l]] <- sample(c(0,1),
                                                             replace=TRUE,
                                                             prob=c(1-p/log(exp(1)-2+l),
                                                                    p/log(exp(1)-2+l)),
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
  
  niche_optima = round(runif(S, envMin+(envMax-envMin)/20, envMax-(envMax-envMin)/20),2)
  write.csv2(niche_optima, paste0(simPath,"niche_optima.csv"), row.names = FALSE)
  
  
  ##################################################################################################  
  ####### Run GLV abiotic GR
  #### Run abiotic GLV GR
  ####
  
  
  glv.finalStates.abioticGR <- mclapply(1:nRep,
                                        function(x){
                                          final_state = sapply(envs, function(e){
                                            glv.out.abioticGR <- simGLV(Stroph, IntMat, spNames, 
                                                                        reduceGR = TRUE, env = e,
                                                                        niche_breadth = niche_breadthGR,
                                                                        niche_optima = niche_optima)
                                            final_state <- t(data.frame(PA = as.numeric(glv.out.abioticGR[,ncol(glv.out.abioticGR)]>0)))
                                            final_abundance <- t(data.frame(AB = glv.out.abioticGR[, ncol(glv.out.abioticGR)]))
                                            return(list(PA = final_state,
                                                        AB = final_abundance))})
                                          final_state_PA = apply(final_state, MARGIN = 2,
                                                                 FUN = function(x){x[[1]]} )
                                          final_state_AB = apply(final_state, MARGIN = 2,
                                                                 FUN = function(x){x[[2]]} )
                                          
                                          return(list(final_state_PA = final_state_PA,
                                                      final_state_AB = final_state_AB))
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
  
  if(all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)])>0) &
     all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)]) < nrow(glv.finalStates.abioticGR.merged))){
     write.csv2(glv.finalStates.abioticGR.merged, paste0(simPath,"glv_finalStates_abioticGR.csv"))
  }else{ 
    cat(paste0("At least one species was never or was always present in",simPath,"\n"))
  }
  
  
  #Abundance-> only needed to compute the potential niche
  glv.finalAbundances.abioticGR <- lapply(glv.finalStates.abioticGR, function(final_state){
    # Set NA values to 0 (computation errors due to too small values)
    final_state$final_state_AB[is.na(final_state$final_state_AB)] <- 0
    final_state$final_state_AB <- as.data.frame(t(final_state$final_state_AB), row.names = envs)
    colnames(final_state$final_state_AB) <- spNames
    return(final_state$final_state_AB)})
  
  glv.meanAbundances.abioticGR <- sapply(Reduce(rbind, glv.finalAbundances.abioticGR),
                                         function(x)weighted.mean(x, w=x>0))
  
  if(all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)])>0) &
     all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)])<nrow(glv.finalStates.abioticGR.merged))){
  write.csv2(glv.meanAbundances.abioticGR, paste0(simPath,"glv.meanAbundances.abioticGR.csv"))
  }
  
  
  
  
  ###################################################################################################################
  ####################### Potential niche
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
  
  glv.fundNicheTh.abioticGR = t(sapply(envs,
                                       function(e)1-ptruncnorm(1-exp(-(e-niche_optima)**2/(2*niche_breadthGR**2)) - glv.meanAbundances.abioticGR %*% (IntMat*PredMat),
                                                               mean=Reduce('+', b.fundNiche),
                                                               sd = Reduce('+', sd_noise.fundNiche), 
                                                               a = ifelse(trophL==1,0,-Inf),
                                                               b=ifelse(trophL>1,0,Inf))))
  colnames(glv.fundNicheTh.abioticGR) = spNames
  if(all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)])>0) &
     all(colSums(glv.finalStates.abioticGR.merged[,-c(1:3)])<nrow(glv.finalStates.abioticGR.merged))){
  write.csv2(glv.fundNicheTh.abioticGR, paste0(simPath,"glv.fundNicheTh.abioticGR.csv"))
  }
}

#############################################################################################
##### Lotka-Volterra
#############################################################################################

#' @title Simulate time series with the generalized Lotka-Volterra model
#'
#' @description Simulate a community time series using the generalized Lotka-Volterra model, defined as
#' \eqn{\frac{dx}{dt}=x(b+Ax)}, where x is the vector of species abundances, A is the interaction matrix
#' and b the vector of growth rates.
#'
#' @param N species number
#' @param A interaction matrix
#' @param b growth rates
#' @param y initial abundances
#' @param tstart initial time point
#' @param tend final time point
#' @param tstep time step
#' @param death.t death threshold for species abundance
#' @return a matrix with species abundances as rows and time points as columns, column names give time points
#' @seealso \code{\link{ricker}} for the Ricker model
#' @examples
#' @export

glv <- function (N = 4, A, b = runif(N), y = runif(N), tstart = 0, tend = 100, tstep = 0.1, death.t=1e-10) 
{
  parms = cbind(b, A)
  parms = cbind(rep(N, N), parms)
  times <- seq(tstart, tend, by = tstep)
  
  rootfun <- function (t, y, pars) {
    return(y - death.t)
  }
  ## sets very low abundances to 0
  eventfun <- function(t, y, pars) {
    y[y - death.t<0] <- 0
    return(y)
  }
  
  commtime <- deSolve::lsoda(y, times, glvsolve, parms, 
                             rootfun=rootfun, events=list(func = eventfun, root = TRUE))
  time = commtime[, 1]
  commtime = commtime[, 2:ncol(commtime)]
  commtime = t(commtime)
  colnames(commtime) = time
  return(commtime)
}

# ==============================================================
# Equations (generalized Lotka-Volterra)
# ==============================================================

# matrix formulation of the ODE set
# t: current simulation time
# y: vector with current values of state variables (initial conditions)
# parms: parameter values
#
glvsolve<-function(t, y, parms){
  N=parms[1,1]  # species number
  b=parms[,2]   # vector of growth rates
  A=parms[,3:(N+2)] # interaction matrix
  if (length(parms[1,]) == (N+3)){
    c=parms[,(N+3)]
    b = b + c
  }
  else if (length(parms[1,]) > (N+3)){
    c=parms[,(N+3)]
    count = parms[1,(N+6)]
    times = parms[,(N+4)]
    durations = parms[,(N+5)]
    intervals = list()
    for (i in 1:count){
      tstart = times[i]
      tend = times[i] + durations[i]
      if ((t>tstart) & (t<tend)){
        b = b + c
        #print(b)
      }
    }
  }
  dydt <- y*(b+A %*% y)
  list(dydt)
}

simGLV <- function(Stroph, IntMat, spNames, reduceGR=F, reduceK=F, env=NULL,
                   niche_optima=NULL, niche_breadth=0.05, tend=200, S = nrow(IntMat)){
  #S = nrow(IntMat)
  b = rep(c(0.5,-0.2), times=c(Stroph[1], S-Stroph[1]))     # intrinsic growth rates (+ for basal species, - for predators)
  bNoise = rnorm(n=S, mean=0, sd=0.1)                       # add some noise
  while (any(b*(b+bNoise)<0)){bNoise = rnorm(n=S, mean=0, sd=0.1)}  # avoid wrong signs
  b = b + bNoise
  
  # Abiotic constraint on carrying capacities
  if (reduceK){
    diag(IntMat) <- diag(IntMat)*pmin(exp((niche_optima-env)**2/(2*niche_breadth**2)), 1e10)
  }
  
  # Abiotic constraint on growth rates
  if (reduceGR){
    if (reduceK){  # if both constraints are present, basal species are impacted only through K
      b[-1:Stroph[1]] <- b[-1:Stroph[1]] + pmax(exp(-(niche_optima-env)**2/(2*niche_breadth**2)), 1e-15)-1
    }else{
      b <- b + pmax(exp(-(niche_optima-env)**2/(2*niche_breadth**2)), 1e-15)-1
    }
  }
  yini = runif(S)    # initial abundances
  
  glv.out = glv(N=S, A=t(IntMat), b=b, y=yini, tend=tend)  # run simulation
  row.names(glv.out) <- spNames
  
  return(glv.out)
}
###############################################################################################
####### Ricker functions
###############################################################################################

ricker <- function (N, A, K = rep(0.1, N), y = runif(N), sigma = 0.05, mu=0.05, 
                    K.trend = NA, tend = 100, death.t = NA, tskip = 0, explosion.bound = 10^8, 
                    perturb = NULL) 
{
  if (length(y) != N) {
    stop("y needs to have N entries.")
  }
  if (nrow(A) != N || ncol(A) != N) {
    stop("A needs to have N rows and N columns.")
  }
  if (length(K) != N) {
    stop("K needs to have N entries.")
  }
  if (length(K.trend) > 1 && length(K.trend) != N) {
    stop("K.trend needs to have N entries.")
  }
  if (tskip >= tend) {
    stop("There are as many or more time points to skip than time points specified.")
  }
  out = matrix(nrow = N, ncol = tend - tskip)
  out[, 1] = y
  perturbCounter = 1
  durationCounter = 1
  K.copy = K
  perturbationOn = FALSE
  for (t in 2:tend) {
    if (sigma > 0) {
      b = rlnorm(N, meanlog = 0, sdlog = sigma)
    }
    else {
      b = rep(1, N)
    }
    if (!is.null(perturb)) {
      if (perturb$times[1] == 1) {
        stop("Please do not specify a perturbation at the first time point.")
      }
      applied = applyPerturbation(perturb, t = t, perturbCounter = perturbCounter, 
                                  durationCounter = durationCounter, perturbationOn = perturbationOn, 
                                  ori.growthrates = K.copy, abundances = y)
      y = applied$abundances
      K = applied$growthrates
      durationCounter = applied$durationCounter
      perturbCounter = applied$perturbCounter
      perturbationOn = applied$perturbationOn
    }
    if (length(K.trend) > 1) {
      K.onepercent = K.copy/100
      K.percent = K.trend * 100 * K.onepercent
      K = K + K.percent
      negative.K.indices = which(K < 0)
      K[negative.K.indices] = 0
    }
    
    # Change the dynamics to get more meaningful results in a predation context
    #y = b * y * exp(A %*% (y - K))
    y = ifelse(y==0, 0, b * y * exp(A %*% y - K*diag(A) - mu*(1+diag(A))))
    
    if (max(y) > explosion.bound) {
      print("Explosion!")
      res = c(-1, which(y == max(y)))
      return(res)
    }
    if (!is.na(death.t) && death.t > 0) {
      y[y < death.t] = 0
    }
    else if (length(y[y < 0]) > 0) {
      stop("Species below 0!")
    }
    if (t > tskip) {
      out[, t - tskip] = y
    }
  }
  return(out)
}

simRicker <- function(Stroph, IntMat, spNames, K=rep(0.1, S), sigma=0.05, death.t=NA,
                      reduceK=F, env=NULL, niche_optima=NULL, niche_breadth=0.05, tskip=0, tend=200){
  S=nrow(IntMat)
  
  # Abiotic constraint on carrying capacities
  if (reduceK){
    K <- K*pmax(exp(-(niche_optima-env)**2/(2*niche_breadth**2)), 1e-15)
  }
  
  yini = abs(rnorm(S, mean=ifelse(trophL==1, K, mean(K)), sd=ifelse(trophL==1, K, mean(K))/3))    # initial abundances
  
  R_IntMat = ricker(N=S, A=IntMat, K=K, y=yini, sigma=sigma, K.trend=NA, tend=tend, 
                    death.t=death.t, tskip=tskip, explosion.bound=10^8, perturb=NULL)
  row.names(R_IntMat) <- spNames
  
  return(R_IntMat)
}


###############################################################################################
####### SOI functions
###############################################################################################

simSOI <- function(Stroph, I, IntMat, spNames, m.vector = runif(nrow(IntMat)), e.vector = runif(nrow(IntMat)), 
                   increaseE=FALSE, basal=FALSE, env=NULL, niche_optima=NULL, niche_breadth=0.3, tend=200){
  S=nrow(IntMat)
  
  # Abiotic constraint on extinction rates
  if (increaseE){
    if (basal){
      e.vector[1:Stroph[1]] <- e.vector[1:Stroph[1]]*exp((niche_optima[1:Stroph[1]]-env)**2/(2*niche_breadth**2))
    }else{
      e.vector <- e.vector*exp((niche_optima-env)**2/(2*niche_breadth**2))
    }
  }
  
  soi.out = seqtime::soi(N=S, I=I, A=t(IntMat), m.vector=m.vector, e.vector=e.vector, K=2, tend=tend)
  
  row.names(soi.out) <- spNames
  
  return(soi.out)
}


###############################################################################################
####### Virtual Com functions
###############################################################################################


simulate_community <- function(
  env = runif(500, 0, 100), #environment value at each site
  niche_optima  = seq(2, 98, 5), #niche means
  niche_breadth = 20, #niche variance
  comp_inter = NA, fac_inter = NA, #competition or facilitation
  beta_env = 1, beta_comp = 5, beta_fac = 0, beta_abun = 0, #filters weights
  years = 20, #number of simulated time-steps
  K = 40  #number of individuals in the community
) {
  sim_com <- function(env, niche_breadth, niche_optima, comp_inter, fac_inter, beta_env,
                      beta_comp, beta_fac, beta_abun, years, K)
  {
    n_sp = length(niche_optima)
    
    if (length(comp_inter) == 1) comp_inter = matrix(comp_inter, n_sp, n_sp) #comp_inter and fac_inter are input                                                                                    matrices of the simulations
    if (length(fac_inter)  == 1) fac_inter  = matrix(fac_inter, n_sp, n_sp)
    species_comp <- comp_inter #species_comp is the competition matrix
    species_fac <- fac_inter #species_fac is the facilitation matrix
    
    log_p_env <- sapply(niche_optima, dnorm, mean = env, sd = niche_breadth, log = TRUE)
    log_p_env <- log_p_env - log(dnorm(0) / 10)
    
    community <- factor(
      x      = sample(seq_along(niche_optima), K, replace = TRUE),
      levels = seq_len(n_sp)
    )
    abund <- table(community)
    
    for (j in seq_len(years)) {
      for (k in seq_len(K)) {
        #all different filters
        f_comp <- 1 - colSums(species_fac[community, ]) / K
        p_comp <- 1 - colSums(species_comp[community, ]) / K
        
        p_all <- exp(beta_env * log_p_env + 
                       beta_comp * log(p_comp) - 
                       beta_fac * log(f_comp) + 
                       log(1 + beta_abun * abund))
        
        p_all <- ifelse(is.na(p_all), min(p_all, na.rm = TRUE), p_all)
        
        if (all(is.na(p_all)) || identical(min(p_all), max(p_all))) p_all = NULL
        if (any(is.infinite(p_all))) {
          community[sample(K, 1)] <- sample(seq_len(n_sp)[p_all == Inf], 1) #sampling the new community
        } else {
          community[sample(K, 1)] <- sample(n_sp, 1, prob = p_all)
        }
        abund <- table(community)
      }
    }
    return (as.integer(abund) > 0)
  }
  ans <- mclapply(
    env, sim_com, niche_breadth, niche_optima, comp_inter, fac_inter,
    beta_env, beta_comp, beta_fac, beta_abun, years, K, mc.cores = detectCores()
  )
  ans <- do.call(rbind, ans)
  ans <- cbind(ans, env)
  sp_labs <- paste0(
    "sp_", gsub(" ", 0, format(seq_along(niche_optima), width = 2))
  )
  colnames(ans) <- c(sp_labs, "env")
  return (as.data.frame(ans))
}


loadData <- function (filePath){
  df.merged <- read.csv2(filePath, check.names = FALSE)
  df.list <- lapply(unique(df.merged$datasets), function(i) {
    df <- subset(df.merged, datasets == i, select = -c(datasets))
    rownames(df) <- df$rowN
    return(df[-(1:2)])
  })
  return(df.list)
}