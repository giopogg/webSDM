### Functions to run statistical analyses on the simulated communities, Giovanni Poggiato and Jérémy Andréoletti

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

plotDistributions <- function (SIM, CV=T, prob.cov=T, plotprey=FALSE, plotpred=FALSE, RN=T,
                               main="", filename=NULL){
  
  nEnv = nrow(SIM$X)/(nRep*nbMerge)
  envs = SIM$X[1:nEnv, ifelse(intercept, -1, -(K+1))]
  if(is.null(SIM$fundNiche) & !RN) stop("No theoretical fundamental niche available")
  
  proba = list()

  if(CV){
    # if realised niche
    if(RN){
    #tSDM
    if(prob.cov){

        proba$pred = SIM$pCV.mean.stan_prob
        proba$est.02 = SIM$pCV.qinf.stan_prob
        proba$est.97 = SIM$pCV.qsup.stan_prob
      
    }else{
    
      proba$pred = SIM$pCV.mean.stan_bin
      proba$est.02 = SIM$pCV.qinf.stan_bin
      proba$est.97 = SIM$pCV.qsup.stan_bin
      
    }
    }else{#if FN
      
      proba$pred = SIM$pFundCV.mean.stan
      proba$est.02 = SIM$pFundCV.qinf.stan
      proba$est.97 = SIM$pFundCV.qsup.stan
      
    }
    #SDMs (realised or fundamental niche is the same!)
      proba$SDM.pred = SIM$SDM.pCV.mean.stan
      proba$SDM.est.02 = SIM$SDM.pCV.qinf.stan
      proba$SDM.est.97 = SIM$SDM.pCV.qsup.stan
      
  }else{#if not CV
    
    # if realised niche
    if(RN){
      #tSDM
      if(prob.cov){
        
        proba$pred = SIM$p.mean.stan_prob
        proba$est.02 = SIM$p.qinf.stan_prob
        proba$est.97 = SIM$p.qsup.stan_prob
        
      }else{
        
        proba$pred = SIM$p.mean.stan_bin
        proba$est.02 = SIM$p.qinf.stan_bin
        proba$est.97 = SIM$p.qsup.stan_bin
        
      }
    }else{#if FN
      
      proba$pred = SIM$pFund.mean.stan
      proba$est.02 = SIM$pFund.qinf.stan
      proba$est.97 = SIM$pFund.qsup.stan
      
    }
    #SDMs (realised or fundamental niche is the same!)
    proba$SDM.pred = SIM$SDM.p.mean.stan
    proba$SDM.est.02 = SIM$SDM.p.qinf.stan
    proba$SDM.est.97 = SIM$SDM.p.qsup.stan
    
  }
    
    
  realizedNiche = SIM$prob
  
  if(RN){
  Niche = SIM$prob
  }else{
  Niche = SIM$fundNiche
  }
  

  # build niche table
  table=data.frame(env=envs,rbind(data.frame(stack(Niche),type=rep("True",nrow(stack(Niche)))),
                                      data.frame(stack(as.data.frame(proba$pred)),type=rep("tSDM",nrow(stack(Niche)))),
                                      data.frame(stack(as.data.frame(proba$SDM.pred)),type=rep("SDM",nrow(stack(Niche))))))
    
  if(plotprey){
    Preys=data.frame()
    for(s in 1:S){
      for (prey in community$trophic.links[community$trophic.links$consumer == spNames[nameOrder[s]],]$resource){
          intStrength = abs(IntMat[prey,spNames[nameOrder[s]]])/(strengthBI*2)  # strength of the biotic interaction
          Preys = rbind(Preys, data.frame(focal=paste0("Y", gsub("Sp|\\..*", "", spNames[nameOrder[s]])),env=envs,
                                          value=realizedNiche[,paste0("Y", gsub("Sp|\\..*", "", prey))]*intStrength,
                                          predprey=paste0("Y", gsub("Sp|\\..*", "", prey))))
          
        }
      }
    }
  if(plotpred){
    Pred=data.frame()
    for(s in 1:S){
      for (pred in community$trophic.links[community$trophic.links$ressource == spNames[nameOrder[s]],]$consumer) {
          intStrength = abs(IntMat[prey,spNames[nameOrder[s]]])/(strengthBI*2)  # strength of the biotic interaction
          Pred = rbind(Pred, data.frame(focal=paste0("Y", gsub("Sp|\\..*", "", spNames[nameOrder[s]])),env=envs,
                                        value=realizedNiche[,paste0("Y", gsub("Sp|\\..*", "", prey))]*intStrength,
                                        predprey=paste0("Y", gsub("Sp|\\..*", "", prey))))
          
        }
      }
      
    }
    
  for(s in 1:S){
      assign(paste0("tableFoc",s), table[which(table$ind==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),])
      if(plotprey){  
        assign(paste0("tablePrey",s), Preys[which(Preys$focal==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),])
      }
      if(plotpred){
        assign(paste0("tablePrey",s), Pred[which(Pred$focal==paste0("Y", gsub("Sp|\\..*", "", spNames[s]))),])
      }
      if(RN){
        cols = c("grey75",alpha("mediumorchid",1),"darkblue")
      }else{
        cols = c(alpha("royalblue1",0.5),"plum4","darkblue")
      }
      assign(paste0("p",s),ggplot() +
               geom_line(data=get(paste0("tableFoc",s))[get(paste0("tableFoc",s))$type=="True",],aes(x=env,y=values),lwd=3, col=cols[1])+
               geom_line(data=get(paste0("tableFoc",s))[get(paste0("tableFoc",s))$type=="tSDM",],aes(x=env,y=values),lwd=2, lty=1, col=cols[2])+
               geom_line(data=get(paste0("tableFoc",s))[get(paste0("tableFoc",s))$type=="SDM",],aes(x=env,y=values),lwd=1, lty=2, col=cols[3]) +
               theme_classic()+ggtitle(spNames[s]) + geom_vline(xintercept = SIM$B[paste0("Y", gsub("Sp|\\..*", "", spNames[s])),2],lty=2,col="red",lwd=1,alpha=0.2) +
               xlab("Env") + ylab("Probability of presence")
      )
      
      if(plotprey){
        assign(paste0("p",s), get(paste0("p",s)) + geom_line(data=get(paste0("tablePrey",s)), aes(x=env,y=value, group=predprey),color="green",linetype=2,alpha=0.8))
      }
      if(plotpred){
        assign(paste0("p",s), get(paste0("p",s)) + geom_line(data=get(paste0("tablePred",s)), aes(x=env,y=value, group=predprey),color="red",linetype=2,alpha=0.8))
        
      }
      
    }
    
    #p=eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p",1:S),collapse=","),"),nrow=4,top=main)")))
  eval(parse(text=paste0("grid.arrange(grobs=list(",paste(paste0("p",1:S),collapse=","),"),nrow=4,top=main)")))
  eval(parse(text=paste0("p=arrangeGrob(grobs=list(",paste(paste0("p",1:S),collapse=","),"),nrow=4,top=main)")))
  if(!is.null(filename)) ggsave(filename = filename,p,width=20,height=10, dpi = 150, units = "in")

  #return(p)
}




compute_TL_laplacian <- function(G){
  # recursive function on connected components of G
  if (igraph::vcount(G) == 1) return(setNames(0, igraph::V(G)$name))
  A = as.matrix(igraph::get.adjacency(G))
  names_loc = rownames(A)
  u  = igraph::degree(G)
  v =  igraph::degree(G,mode ='in') - igraph::degree(G,mode = 'out')
  
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u, ncol=length(u)) - A - t(A) # matrix is made symmetric!
  
  TL_vec = Matrix::solve(L,v)
  TL_vec = c(0,TL_vec)
  TL_vec = TL_vec - min(TL_vec)
  names(TL_vec) = igraph::V(G)$name
  return(TL_vec)
}






# Trophic SDM
## Parameters:
# Y : Community data must be a sites x species matrix
# X : Environmental data must be a sites x predictor matrix
# G : Trophic network, igraph
# formulas : list (names of each element are the name of vertices) of one-hand formulas whose names are the same as in X. 
# sp.formula : formula for composite variable. It HAS TO BE a string but written as the rhs of a formula. It can include quadratic terms as formula (also poly works). Be careful ONLY "richness". For example, you could set sp.formula="richness" or sp.formula="richness+richness^2" or sp.formula="sin(richness^2)" (no, you won't) 
# sp.partition : list where the edges (i.e. their name as string) are grouped. Each element of the list is a group. E.g. sp.partition = list(c("Y_1","Y_2"),"Y_3",c("Y_4","Y_5","Y_6"))
# penalisation: how to penalise the likelihood to reduce the dimension. Each SDM method has its own penalisation type, some SDM method have no available penalisation approach. For now, we have "horshoe" for stan_glm, "elasticnet" for glm and  "coeff.signs" (positive prey interactions and negative predator interactions) for glm and stan_glm
# method: which SDM method to use (available for now: glm (frequentist), stan_glm (full bayesian MCMC), bayesglm (variational approx of the posterior))
# family: the family parameter of the glm function (see glm). family=gaussian for gaussian data or family=binomial(logit) or binomial(probit) for presence-absence data
# fitPreds : logical. Should predators be included as covariates?
# iter : the number of MCMC samples (if bayesian approach)

trophicSDM = function(Y,X,G,formulas=NULL,sp.formula=NULL,sp.partition=NULL,penal=NULL,method="stan_glm",family,fitPreds=FALSE,iter=1000,run.parallel=T,chains=2,verbose=F){
  
  # checks & errors
  if(!is_igraph(G)) stop("G is not an igraph object")
  if(!method %in% c("glm","stan_glm","bayesglm")) stop("the selected method has not yet been implemented or it has been misspelled")
  if(!(is.null(penal) || penal %in% c("horshoe","elasticnet","coeff.signs"))) stop("the selected penalisation has not yet been implemented or it has been misspelled")
  
  # laplacian sorting of the graph (component by component)
  sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=T)]
  
  # initialize empty lists of models
  m = form.all = as.list(vector(length=vcount(G)))
  names(m) = names(form.all) = names(sortedV)
  
  # core part: loop on the species to fit their distribution
  if(run.parallel){
    m=mclapply(1:vcount(G),FUN = function(j){
    if(verbose) print(paste("--- Species", names(sortedV[j]), "---"))
    
    # call a function that does a SDM with j as focal species, automatically finds covariates from G and formula.foc
    temp.mod = SDMfit(focal=names(sortedV[j]), Y, X, G, sp.formula, sp.partition,
                      formula.foc=paste(as.character(formulas[[names(sortedV[j])]]),collapse=" "),
                      method=method, penal=penal, family=family, fitPreds=fitPreds, iter=iter,verbose=verbose,chains=chains)
    
    # assign models and formulas (that now includes biotic variables too)
    temp.mod
  },mc.cores=detectCores())
  
  form.all=lapply(m, `[[`, 2)
  m=lapply(m, `[[`, 1)
  names(m) = names(form.all) = names(sortedV)
  
  
  }else{
  
  
  for(j in 1:vcount(G)){
    if(verbose) print(paste("--- Species", names(sortedV[j]), "---"))
    
    # call a function that does a SDM with j as focal species, automatically finds covariates from G and formula.foc
    temp.mod = SDMfit(focal=names(sortedV[j]), Y, X, G, sp.formula, sp.partition,
                      formula.foc=paste(as.character(formulas[[names(sortedV[j])]]),collapse=" "),
                      method=method, penal=penal, family=family, fitPreds=fitPreds, chains=chains,iter=iter,verbose=verbose)

    # assign models and formulas (that now includes biotic variables too)
    m[[names(sortedV[j])]]=temp.mod$m
    form.all[[names(sortedV[j])]]=temp.mod$form.all
  }
    
 }
  # return values
  list(model=m,form.all=form.all,form.env=formulas,sp.formula=NULL,sp.partition=NULL,data=list(X=X,Y=Y),method=method,penal=penal,iter=iter,family=family,G=G)
}


# SDMfit : the core function that fit a (classic) SDM on the focal species, with its neighbours (preys if bottom-up approach, predators if top-down approach or both) and the environmental covariates as predictors.
## Parameters:
# focal : the vertex of the graph (i.e. the species) to be modelled
# other parameters: see trophicSDM


# nice stuff about formulas https://www.datacamp.com/community/tutorials/r-formula-tutorial

SDMfit=function(focal,Y,X,G,formula.foc,sp.formula,sp.partition,method="bayesglm",family=NULL,penal=NULL,fitPreds=FALSE,iter=1000,chains=2,verbose=T){
  
  # for the particular case of the stan_glm model with "coeff.signs" constraint
  # another stan-for-R package is used (brsm) and formulas have to be adapted
  useBRMS = !is.null(penal) && penal=="coeff.signs" & method=="stan_glm"
  
  
  ## Build environmental part of the formula
    if(useBRMS){
      # "tempX" is an intermediary variable allowing to set lower/upper-bounds on priors
      form.env = "y ~ tempX"
      form.env.brms = paste("tempX",as.character(formula.foc))
    }else form.env = paste("y",as.character(formula.foc))

  
  
  ## Build final formula
  preys = names(neighbors(G,focal,mode=c("out")))
  preds = names(neighbors(G,focal,mode=c("in")))
  
  if (fitPreds){
    # define a separate category for pairs of species that predate each other (2-species loop),
    # as they must only be included once in the regression
    predsANDpreys = intersect(preds, preys)
    preds = setdiff(preds, predsANDpreys)
    preys = setdiff(preys, predsANDpreys)
  }
  
  ### Include preys as covariates
  if (length(preys)>0){
    formulas = buildFormula(form.init=form.env, species=preys, type="Preys", sp.formula=sp.formula, sp.partition=sp.partition, useBRMS)
    form.all = formulas$form.all
    if (useBRMS){
      # intermediary species variables are themselves defined from observed species variables
      form.preys.brms = formulas$form.brms
    }
  }else{
    # the focal species is basal
    form.all = as.character(form.env)
  }
  
  ### Include predators as covariates
  if (fitPreds & length(preds)>0){
    formulas = buildFormula(form.init=form.all, species=preds, type="Preds", sp.formula, sp.partition, useBRMS)
    form.all = formulas$form.all
    if (useBRMS) form.preds.brms = formulas$form.brms
  }
  
  ### Include species that are at the same time preys and predators
  if (fitPreds && length(predsANDpreys)>0){
    formulas = buildFormula(form.init=form.all, species=predsANDpreys, type="PreysAndPreds", sp.formula, sp.partition, useBRMS)
    form.all = formulas$form.all
    if (useBRMS) form.predsANDpreys.brms = formulas$form.brms
  }
  
  
  ## Build the data matrix to be given to the fit functions
  data = data.frame(X,Y)
  data = dplyr::rename(data,y=all_of(focal))  # rename the focal column to be called "y"
  
  
  ## Fit models depending on the chosen method and penalisation
  
  ### GLM
  if(method=="glm"){
    if(is.null(penal)){
      m = glm(form.all, data=data, family=family)
    }else{
      if(penal=="horshoe") stop("Horshoe penalisation is not possible with glm")
      if(penal=="elasticnet"){
        y = as.matrix(model.frame(form.all,data=data)$y)
        x = as.matrix(model.frame(form.all,data=data)[,-1])
        
        m = glmnet(y=y,x=x,family=family$family,
                   lambda=cv.glmnet(y=y,x=x,family=family$family)$lambda.1se)
      }
      if(penal=="coeff.signs"){
        # Constrain the signs of trophic interaction coefficients to remain positive (preys->preds) or negative (preds->preys)
        # species that are at the same time preys and predators are not constrained
        y = as.matrix(model.frame(form.all,data=data)$y)
        x = as.matrix(model.frame(form.all,data=data)[,-1])
        
        # Extract the regression variables and set lower and upper limits
        variables = sapply(str_extract_all(form.all, "X[0-9]*|Y[0-9]*|I\\(([^()]*|\\(([^()]*|\\([^()]*\\))*\\))*\\)")[[1]],
                           function(x)str_extract(x,"X[0-9]*|Y[0-9]*")[[1]])
        lower.limits <- sapply(variables, function(VAR){if (VAR %in% preys) return(0)
          if (grepl("X", VAR) | VAR %in% c(preds,predsANDpreys)) return(-Inf)
          warning(paste("Variable", VAR, "not present in any category"), call.=F)})
        upper.limits <- sapply(variables, function(VAR){if (VAR %in% preds) return(0)
          if (grepl("X", VAR) | VAR %in% c(preys,predsANDpreys)) return(Inf)
          warning(paste("Variable", VAR, "not present in any category"), call.=F)})
        
        m = glmnet(y=y,x=x,family=family$family, lambda=0, thres = 1E-10, 
                   lower.limits=lower.limits, upper.limits=upper.limits)
      }
    }
  }
  
  
  ### stan_glm
  if(method=="stan_glm"){
    refresh=ifelse(verbose,1000,0)
    
    if(is.null(penal)){
      m = stan_glm(form.all, data=data, family=family, iter=iter, prior=normal(0,10),chains=chains,refresh=refresh)
    }else{
      if(penal=="elasticnet")  stop("Elastic net is not stan_glm")
      if(penal=="horshoe"){
        HS_prior = hs(df = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)
        m = stan_glm(form.all, data = data, family = family, prior = HS_prior,
                     iter=iter,chains=chains, refresh = refresh)
      }
      if(penal=="coeff.signs"){
        # constrain the signs of trophic interaction coefficients to remain positive (preys->preds) or negative (preds->preys)
        variables = strsplit(as.character(form.all), "y ~ | \\+ ")[[1]][-1]
        # choose priors depending on the nature of the variables
        priors = Reduce(rbind, lapply(variables, function(VAR){
          VAR <- gsub("\\^", "E", gsub('\\(|\\)', "", VAR))  # correct name for polynomial terms
          if (grepl("X", VAR)) return(set_prior("normal(0,10)", class="b", nlpar=VAR))
          if (sub("temp","",VAR) %in% preys) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))
          if (sub("temp","",VAR) %in% preds) return(set_prior("normal(0,5)", class="b", ub=0, nlpar=VAR))
          if (sub("temp","",VAR) %in% predsANDpreys) return(set_prior("normal(0,5)", class="b", nlpar=VAR))
          if (grepl("richnessPreys",VAR)) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))
          if (grepl("richnessPreds",VAR)) return(set_prior("normal(0,5)", class="b", ub=0, nlpar=VAR))
          if (grepl("richnessPreysAndPreds",VAR)) return(set_prior("normal(0,5)", class="b", nlpar=VAR))
          warning(paste("Variable", sub("temp","",VAR), "not present in any category"), call.=F)
        }))
        
        # each temporary variable depends on our initial variables ("false" non-linear model)
        form.all = bf(form.all, nl=TRUE) + lf(form.env.brms)
        if(length(preys)>0) form.all = form.all + lf(form.preys.brms)[[1]]
        if(fitPreds & length(preds)>0) form.all = form.all + lf(form.preds.brms)[[1]]
        if(fitPreds && length(predsANDpreys)>0) form.all = form.all + lf(form.predsANDpreys.brms)[[1]]
        
        if(verbose) print(c(Formula=form.all))
        
        m = brm(form.all, data=data, family=family, iter=iter, prior=priors, control=list(adapt_delta = 0.9))
      }
    }
  }
  
  # create matrix of coefficient with postmean, post quantiles and whole samples
  
  # return
  list(m=m,form.all=form.all)
}


## Build formula function


buildFormula <- function(form.init, species, type, sp.formula=NULL, sp.partition=NULL, useBRMS){
  #if no composite variables are used
  if(is.null(sp.formula)){
    if(useBRMS){
      # "tempYi" are intermediary variables allowing to set lower/upper-bounds on priors
      form.all = paste(as.character(form.init), paste(paste0("temp",species), collapse="+"), sep="+")
      # intermediary species variables are defined hierarchically from observed species variables
      form.brms = lapply(species, function(sp)as.formula(paste0("temp",sp," ~ 0 + ",sp)))
    } else form.all = paste(as.character(form.init), paste(species,collapse="+"), sep="+")
  }else{
    # if no group definition
    if(is.null(sp.partition)){
      # we replace "richness" with the sum of all species (and all eventual other terms), both in linear and eventually other terms
      form.all=as.character(form.init)
      
      if(grepl("richness",sp.formula)){
        # replace richness by the observed species variables
        if(useBRMS){
          # hierarchical definition
          form.rich = "richness"
          form.brms = paste0("richness ~ 0 + ", gsub("richness",paste0("I(",paste(paste0("temp",species),collapse = "+"),")"),sp.formula))
        } else form.rich = gsub("richness", paste0("I(",paste(species,collapse = "+"),")"), sp.formula)
        form.all = paste(form.all,form.rich,sep="+")
      }
      
      if(grepl("diversity",sp.formula)){
        
        # not yet available
        
      }
    }else{
      # all as above but repeated for every cluster
      form.all=as.character(form.init)
      
      if(grepl("richness",sp.formula)){
        form.rich = NULL
        for(j in 1:length(sp.partition)){
          species.j = species[species %in% sp.partition[[j]]]
          
          if(length(species.j)>1){        # at least 2 species in the cluster
            if(is.null(form.rich)){
              # create form.rich
              if(useBRMS){
                form.rich = paste0("richness",type,j)
                form.brms = list(paste0("richness",type, j, " ~ 0 + ", gsub("richness", paste0("(",paste(species.j,collapse="+"),")"), paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))))
              } else form.rich = gsub("richness", paste0("(",paste(species.j,collapse="+"),")"), paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))
            }else{
              # add the new formula to the already existing form.rich
              if(useBRMS){
                form.rich = paste(form.rich, paste0("richness",type, j), sep="+")
                form.brms = c(form.brms, paste0("richness",type, j, " ~ 0 + ", gsub("richness", paste0("(",paste(species.j,collapse="+"),")"), paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))))
              } else form.rich = paste(form.rich, gsub("richness", paste0("(",paste(species.j,collapse="+"),")"), paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")")),sep="+")
            }
          }else if(length(species.j)==1){  # only 1 species in the cluster
            if(is.null(form.rich)){
              if(useBRMS){
                form.rich = paste0("richness",type,j)
                form.brms = list(paste0("richness",type, j, " ~ 0 + ", species.j))
              } else form.rich = species.j
            }else{
              if(useBRMS){
                form.rich = paste(form.rich, paste0("richness",type, j), sep="+")
                form.brms = c(form.brms, paste0("richness",type, j, " ~ 0 + ", species.j))
              } else form.rich = paste(form.rich, species.j, sep="+")
            }
          }
        }
        form.all = as.formula(paste(form.all,form.rich,sep="+"))
      }
    }
  }
  if (useBRMS){
    return (list(form.all=form.all, form.brms=lapply(form.brms,as.formula)))
  }
  return (list(form.all=form.all))
}


###################################################################################################
### trophicSDM_predict
#trophicSDM_predict : the function that predicts the output of the trophicSDM model, eventually
#depending on a new set of environmental covariates. Several parameters are available depending 
#on the users' needs.

## Parameters:
# m : the fitted model
# Xnew : Environmental data must be a sites x predictor matrix. If NULL, the fitted environmental predictors are used (i.e. in-sample prediction)
# binary.resp : only useful with presence-absence variable. TRUE if we also want binary-output, FALSE if only probabilities are given.
# prob.cov : TRUE if the predicted probabilities of presence of the preys are used to predict species. FALSE if the randomly (depending of predicted probabilities of course) generated presence-absence are used
# mode : type of graph neighbors used for prediction, "out" for preys only, "in" for predators only, "all" for both
# pred.sample : the total number of sample from the predictive distribution of each species at each site
# error_prop_sample : fine tuning parameter for error propagation. If 1, for each of the predictive samples of the preys, we draw one sample of the predicted distribution for the focal species. 
#                     If greater than one, we sample pred.sample/error_prop_sample of the predictive distribution of each species (N.B. carefully to guarantee that the samples are consistent across them),
#                     and then for each of the selected samples, we draw error_prop_sample from the predictive distribution of the focal species, for each selected sample of the predictors.

trophicSDM_predict=function(m,Xnew,binary.resp=NULL,prob.cov=NULL,mode="all",pred_samples=1000,error_prop_sample=10,verbose=F){
  
  # checks & errors
  if(m$method=="glm" & (pred_samples != 1 | error_prop_sample != 1)){stop("glm requires pred_sample and error_prop_sample both =1")}
  
  if(is.null(Xnew)) Xnew = m$data$X
  
  family = m$family
  n = nrow(Xnew)
  S = ncol(m$data$Y)
  G = m$G
  sp.prediction = as.list(vector(length=vcount(G)))
  sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=T)]
  names(sp.prediction) = sortedV$name
  
  if(family$family %in% c("bernoulli", "binomial")){
    
    if(is.null(binary.resp))  message("Please specify the parameter binary.resp. TRUE if you also want binary predictions, FALSE otherwise")
    if(is.null(prob.cov))  message("Please specify the parameter prob.cov. TRUE if you want to use probabilities instead of presence-absence as a proxy for other species predictions. FALSE otherwise")
    
    Neighb = lapply(sortedV, function(sV)neighbors(G,sV,mode))
    
    # core loop on species (as in trophicSDM)
    for(j in 1:vcount(G)){
      if(verbose){print(paste("--- Species", names(sortedV[j]), "---")); print(m$model[[names(sortedV[j])]]$formula)}
      #warning("", call.=F)
      
      # neighbor species
      neighb.sp = Neighb[[sortedV$name[j]]]
      # neighbor species that have already been predicted
      known.sp = sortedV[which(sapply(sp.prediction, is.list))]
      known.neighb = neighb.sp[which(neighb.sp %in% known.sp)]
      # neighbor species that haven't been predicted yet
      unknown.neighb = neighb.sp[which(!(neighb.sp %in% known.sp))]
      
      # select the predictive samples from the preys (that have already been predicted!!)
      
      # create newdata: will include the (eventually new) environmental predictors and the predicted distribution of the preys
      # it's an array where we have the predictive samples on the third dimension.
      # we will therefore apply SDMpredict to each of the matrices (i.e. first two dimensions) in order to correctly propagate the error
      newdata = array(data=NA,dim=c(n,(ncol(Xnew)+length(known.neighb)+length(unknown.neighb)),pred_samples/error_prop_sample))
      colnames(newdata) <- c(colnames(Xnew),names(known.neighb),names(unknown.neighb))
      
      # fill the abiotic variables
      newdata[,1:ncol(Xnew),] <- as.matrix(Xnew)
      
      # fill the biotic part of new data
      
      ## species that have already been predicted
      if (length(known.neighb)>0) for(k in 1:length(known.neighb)){
        if(prob.cov){
          # if prob.cov=TRUE, use the species predicted probabilities of presence
          newdata[,ncol(Xnew)+k,] <- sp.prediction[[names(known.neighb[k])]]$predictions.prob[,seq(from=1,by=error_prop_sample,to=pred_samples)]# here select the random samples!  
        }else{
          newdata[,ncol(Xnew)+k,] <- tryCatch(sp.prediction[[names(known.neighb[k])]]$predictions.bin[,seq(from=1,by=error_prop_sample,to=pred_samples)], error=function(cond){message(paste("dim(sp.prediction[[names(known.neighb[k])]]$predictions.bin):", sp.prediction[[names(known.neighb[k])]]$predictions.bin));return(NA)})# here select the random samples!
        }
      }
      
      ## species that have not already been predicted so another method has to be chosen
      if (length(unknown.neighb)>0) for(k in 1:length(unknown.neighb)){
        if(prob.cov){
          # if prob.cov=TRUE, choose a method to define the unknown species probability of presence
          # M1 : uniform draw
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- runif(n=length(newdata[,ncol(Xnew)+length(known.neighb)+k,]))# here select the random samples!  
          # M2 : mean observed prevalence
          newdata[,ncol(Xnew)+length(known.neighb)+k,] <- mean(m$model[[names(sortedV[j])]]$data[,names(unknown.neighb)[k]])
          # M3 : true Presence/Absence (for tests only)
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- SIMlist.test$GLV_abioticGR$Y[,names(unknown.neighb)[k]]
          # M4 : true realized niche (for tests only)
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- SIMlist$GLV_abioticGR$prob[,names(unknown.neighb)[k]]  
        }else{
          # if prob.cov=FALSE, choose a method to define the unknown species presence/absence
          # M1 : uniform sampling
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- sample(c(0,1),size=dim(newdata[,ncol(Xnew)+length(known.neighb)+k,]), replace=T)# here select the random samples!
          # M2 : sampling based on mean observed prevalence
          newdata[,ncol(Xnew)+length(known.neighb)+k,] <- rbinom(n=length(newdata[,ncol(Xnew)+length(known.neighb)+k,]), size=1, prob=mean(m$model[[names(sortedV[j])]]$data[,names(unknown.neighb)[k]]))
          # M3 : true Presence/Absence (for tests only)
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- SIMlist.test$GLV_abioticGR$Y[,names(unknown.neighb)[k]]
          # M4 : true realized niche (for tests only)
          # newdata[,ncol(Xnew)+length(known.neighb)+k,] <- SIMlist$GLV_abioticGR$prob[,names(unknown.neighb)[k]]
        }
      }
      
      # apply the function SDMpredict to each layer of the array (MARGIN=3)
      pred.array = apply(newdata, MARGIN=3, FUN = function(x){SDMpredict(focal=names(sortedV[j]),m=m,newdata=x,pred_samples=error_prop_sample,binary.resp=binary.resp,prob.cov=prob.cov)})
      # unlist and format
      sp.prediction[[names(sortedV[j])]] = list(predictions.prob=NULL,predictions.bin=NULL)
      sp.prediction[[names(sortedV[j])]]$predictions.prob = do.call(cbind,lapply(pred.array,FUN=function(x) x$predictions.prob))
      if(!prob.cov | binary.resp) {
        sp.prediction[[names(sortedV[j])]]$predictions.bin = do.call(cbind,lapply(pred.array,FUN=function(x) x$predictions.bin))
      }
    }
    if(!binary.resp){
      # if binary.resp=FALSE, remove the binary predictions (that were eventually needed for computations if prob.cov=FALSE)
      for(j in 1:length(sp.prediction)){
        sp.prediction[[j]]$predictions.bin=NULL
      }
    }
    return(list(sp.prediction=sp.prediction,binary.resp=binary.resp,prob.cov=prob.cov,pred_samples=pred_samples,error_prop_sample=error_prop_sample))
    
  }
  
  # Just the same of above but predictions are directly given (no need to choose between prob or not)
  if(family$family=="gaussian"){
    
    for(j in 1:vcount(G)){ # eventually modify with apply
      # check if basal or not
      basal=ifelse(length(neighbors(G,v=sortedV[j],mode="out"))==0,TRUE, FALSE)
      
      if(basal){
        sp.prediction[[names(sortedV[j])]]=SDMpredict(focal=names(sortedV[j]),m=m,newdata=Xnew,pred_samples=pred_samples,binary.resp=binary.resp,prob.cov=prob.cov)
      }else{
        neigh.sp = neighbors(G,v=sortedV[j],mode="out")
        newdata=array(data=NA,dim=c(n,(ncol(Xnew)+length(neigh.sp)),pred_samples/error_prop_sample))
        colnames(newdata)=c(colnames(Xnew),names(neigh.sp))
        
        # fill the abiotic variables
        newdata[,1:ncol(Xnew),]=as.matrix(Xnew)
        
        for(k in 1:length(neigh.sp)){
          newdata[,ncol(Xnew)+k,]=sp.prediction[[names(neigh.sp[k])]][,seq(from=1,by=error_prop_sample,to=pred_samples)]# here select the random samples!
        }
        
        pred.array=apply(newdata,MARGIN=3,FUN = function(x){list(SDMpredict(focal=names(sortedV[j]),m=m,newdata=x,pred_samples=error_prop_sample,binary.resp=binary.resp,prob.cov=prob.cov))})
        sp.prediction[[names(sortedV[j])]] = do.call(cbind,lapply(pred.array,FUN=function(x) x[[1]]))
      }
    }
    
    return(list(sp.prediction=sp.prediction,pred_samples=pred_samples,error_prop_sample=error_prop_sample))
    
  }
  
}

## Parameters:
# m : the model
# focal : the focal species to be predicted
# newdata: the corresponding predictors (formatted in trophicSDM_predict)
# pred_samples : the number of samples to be predicted. Notice that it corresponds to the third dimension of the newdata array.
# for the other parameters see trophicSDM_predict


SDMpredict=function(m,focal=focal,newdata,pred_samples,binary.resp,prob.cov){
  
  family=m$family$family
  method=m$method
  sortedV= topo_sort(G,mode="in")
  n=nrow(m$data$X)
  S=ncol(m$data$Y)
  
  # If presence/absence data
  
  if(family %in% c("bernoulli", "binomial")){
    # initialize vectors
    predictions.prob=matrix(data=NA,nrow=n,ncol=pred_samples)
    if(!prob.cov | binary.resp) {predictions=matrix(data=NA,nrow=n,ncol=pred_samples)}else{predictions=NULL}
    
    if(method=="stan_glm"){
      
      # retrieve the expected predicted probability of presence
      ## ! for the particular case of the stan_glm model with a constraint on the coefficients
      ## ! signs another stan-for-R package is used (brsm) and the code has to be adapted
      if (class(m$model[[focal]])[1]=='brmsfit'){  # brms package
        predictions.prob=t(posterior_epred(m$model[[focal]], newdata=as.data.frame(newdata), nsamples=pred_samples))
      }else{  # rstan package
        predictions.prob=t(posterior_epred(m$model[[focal]], newdata=as.data.frame(newdata), draws=pred_samples))
      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov | binary.resp) {
        if (class(m$model[[focal]])[1]=='brmsfit'){  # brms package
          predictions=t(posterior_predict(m$model[[focal]], newdata=as.data.frame(newdata), nsamples=pred_samples))
        }else{  # rstan package
          predictions=t(posterior_predict(m$model[[focal]], newdata=as.data.frame(newdata), draws=pred_samples))
        }
      }
    }
    
    
    if(method=="glm")  {
      if(pred_samples!=1) stop("pred_sample must be 1 if method=glm!")
      
      # retrieve the expected predicted probability of presence
      predictions.prob=as.matrix(predict(m$model[[focal]],type = "response",newdata=as.data.frame(newdata)))
      
      # retrieve a unique presence/absence prediction
      if(binary.resp & !prob.cov) {
        predictions=as.matrix(sapply(as.vector(predictions.prob),FUN=function(x) rbinom(prob=x,size=1,n=1)))
      }                 
    }
    return(list(predictions.prob=predictions.prob,predictions.bin=predictions))
    
  }
  
  # If gaussian data
  
  if(family=="gaussian"){
    
    if(method=="stan_glm"){
      predictions=t(posterior_predict(m$model[[focal]], newdata=as.data.frame(newdata),draws=pred_samples))
    }
    
    if(method=="bayesglm"){
      m.sim=sim(m$model[[focal]],n.sims=pred_samples)
      m.sim=cbind(coef(m.sim),sigma=sigma.hat(m.sim))
      data.temp=as.matrix(cbind(1,model.frame(m$form.all[[focal]][-2], data=as.data.frame(newdata))))
      predictions=apply(m.sim,FUN=function(x) data.temp%*%x[1:(length(x)-1)]+rnorm(n=1,mean=0,sd=x[length(x)]),MARGIN=1)
      
    }
    
    if(method=="glm")  {
      
      if(pred_samples!=1) {
        
        stop("pred_sample must be 1 if method=glm!")}
      
      predictions=as.matrix(predict(m$model[[focal]],type = "response",newdata=as.data.frame(newdata)))
      
    }
    
    return(predictions) 
  }
}



# Cross validation function. Computes prediction from the model in K-fold CV and, if asked, also computes predictions of the fundamental niche
# INPUT:
# m: the fitted model
# K: the number of folds. can be NULL if index is specified
# index: a vector containing the fold assigned to each of the sites

tSDM_CV = function(m,K,partition=NULL, fundNiche=F,prob.cov=T,iter=m$iter,
                   pred_samples=m$iter,run.parallel,verbose=F,
                   error_prop_sample=10,fitPreds=F){
  
  n = nrow(m$data$X)
  S = length(m$model)
  # Create partition (if needed)
  if(!is.null(partition)){if(length(partition) != n) stop("partition must be a vector of length n (the number of sites")
  }else{partition <- sample(1:K,size=n,replace=TRUE,prob=rep(n/K,K))}
  
  
  preds = array(dim=c(n,pred_samples,S))
  
  if(fundNiche){ predsFund = array(dim=c(n,pred_samples,S))}
  
  for(i in 1:K){
    
    train = which(partition != i)
    test = which(partition == i)
    Y = m$data$Y[train,]
    X = m$data$X[train,]
    
    m_K = trophicSDM(Y=Y,X=X,G=m$G,formulas=m$form.env,penal=m$penal,method=m$method,
                     family=m$family,fitPreds=fitPreds,iter=iter,
                     run.parallel = run.parallel, verbose = verbose,
                     chains=2)
    
    pred_K = trophicSDM_predict(m=m_K,Xnew = m$data$X[test,],
                                binary.resp=F,prob.cov=prob.cov,
                                mode=ifelse(fitPreds,"all","out"),pred_samples=pred_samples,
                                error_prop_sample=error_prop_sample)
    
    preds[test,,] = abind(lapply(pred_K$sp.prediction,function(x) x$predictions.prob),along=3)
    
    if(fundNiche){
      
      for(s in 1:S){
        
        sp_mod = m$model[[s]]
        
        
        #fix all species to one (i.e. present)
        newdata =  cbind(m$data$X[test,],
                         data.frame(matrix(1,nrow=length(test),ncol=ncol(m$data$Y),
                                           dimnames= list(NULL,colnames(m$data$Y))))
        )
        pred_temp = SDMpredict(m=SIM$m_stan,focal=names(m$model)[[s]],newdata = newdata,
                               pred_samples = pred_samples, binary.resp = F,prob.cov = prob.cov)
        
        predsFund[test,,s] = pred_temp$predictions.prob
        
        
      }
      
    }
    
    print(paste0("Fold ", i, " out of ", K,"\n"))
    
  }
  
  meanPred = apply(preds,mean,MARGIN = c(1,3))
  Pred975 = apply(preds, quantile, MARGIN=c(1,3),0.975)
  Pred025 = apply(preds, quantile, MARGIN=c(1,3),0.025)
  
  colnames( meanPred ) = colnames(Pred975) = colnames(Pred025) = names(pred_K$sp.prediction)
  # Compute Joint TSS and AUC
  
  auc=tss=vector(length=S)
  
  eval = mclapply(1:S,function(x){
    eval=dismo::evaluate(p=meanPred[which(m$data$Y[,colnames(meanPred)[x]]==1),x], a=meanPred[which(m$data$Y[,colnames(meanPred)[x]]==0),x])
    data.frame(auc=eval@auc,tss=max(eval@TPR+eval@TNR-1))
  })
  
  CVmetrics = do.call(rbind,eval)
  CVmetrics$species = colnames(meanPred)
  
  if(fundNiche){
    
    FN.meanPred = apply(predsFund,mean,MARGIN = c(1,3))
    FN.Pred975 = apply(predsFund, quantile, MARGIN=c(1,3),0.975)
    FN.Pred025 = apply(predsFund, quantile, MARGIN=c(1,3),0.025)
    colnames( FN.meanPred ) = colnames(FN.Pred975) = colnames(FN.Pred025) = names(pred_K$sp.prediction)
    
  }
  
  list(meanPred = meanPred, Pred975 = Pred975, Pred025 = Pred025, CVmetrics = CVmetrics, partition = partition, fundNiche = list(FN.meanPred, FN.Pred975, FN.Pred025))

}



tSDM_CV_SIMUL = function(mod, K, fundNiche = F, prob.cov = T,iter,
                         pred_samples, error_prop_sample = 10,
                         fitPreds = F,run.parallel = T, nEnv, verbose=F, envCV=F, chains=2){
  
  X = mod$data$X[1:nEnv,"X1"]
  S = length(mod$model)
  
  if(envCV){
  partition = rep(1:K,each = round(nEnv/K))
  
  if(length(partition)<nEnv) partition = c(partition, rep(K,nEnv-length(partition)))
  }else{
    partition = sample(1:K,size=nEnv,replace=T)
  }
  index_all = vector(length=nrow(mod$data$X))
  
  for(i in 1:K){
    
    index_all[ which(mod$data$X[,"X1"] %in% X[partition==i])] = i

  }
  
  preds = array(dim=c(nEnv,pred_samples,S))
  
  if(fundNiche){ predsFund = array(dim=c(nEnv,pred_samples,S))}
  
  for(i in 1:K){

    train = which(index_all != i)
    test = which(partition == i)

    m_K = trophicSDM(Y=mod$data$Y[train,],X=mod$data$X[train,],G=mod$G,formulas=mod$form.env,penal=mod$penal,method=mod$method,
                     family=mod$family,fitPreds=fitPreds,iter=iter,run.parallel = run.parallel,verbose=verbose,chains=chains)
    

    pred_K = trophicSDM_predict(m=m_K,Xnew = mod$data$X[test,],
                                binary.resp=F,prob.cov=prob.cov,
                                mode=ifelse(fitPreds,"all","out"),pred_samples=pred_samples,
                                error_prop_sample=error_prop_sample, verbose = verbose)
    
    preds[test,,] = abind(lapply(pred_K$sp.prediction,function(x) x$predictions.prob),along=3)
    
    if(fundNiche){
      
      for(s in 1:S){
        sp_mod = mod$model[[s]]
        
        
        #fix all species to one (i.e. present)
        newdata =  cbind(mod$data$X[test,],
                         data.frame(matrix(1,nrow=length(test),ncol=ncol(mod$data$Y),
                                           dimnames= list(NULL,colnames(mod$data$Y))))
        )
        pred_temp = SDMpredict(m=m_K,focal=names(mod$model)[[s]],newdata = newdata,
                               pred_samples = pred_samples, binary.resp = F,prob.cov = prob.cov)
        
        predsFund[test,,s] = pred_temp$predictions.prob
        
      }
      
    }
    
    print(paste0("Fold ", i, " out of ", K,"\n"))
    
  }
  
  meanPred = apply(preds,mean,MARGIN = c(1,3))
  Pred975 = apply(preds, quantile, MARGIN=c(1,3),0.975)
  Pred025 = apply(preds, quantile, MARGIN=c(1,3),0.025)
  
  colnames( meanPred ) = colnames(Pred975) = colnames(Pred025) = names(pred_K$sp.prediction)

  
  if(fundNiche){
    
    FN.meanPred = apply(predsFund,mean,MARGIN = c(1,3))
    FN.Pred975 = apply(predsFund, quantile, MARGIN=c(1,3),0.975)
    FN.Pred025 = apply(predsFund, quantile, MARGIN=c(1,3),0.025)
    colnames( FN.meanPred ) = colnames(FN.Pred975) = colnames(FN.Pred025) = names(pred_K$sp.prediction)
    
  }else{
    FN.meanPred = FN.Pred975 =  FN.Pred025 = NULL
  }
  
  list(meanPred = meanPred, Pred975 = Pred975, Pred025 = Pred025, partition.env = partition, index_all=index_all,
       fundNiche = list(FN.meanPred = FN.meanPred, FN.Pred975 = FN.Pred975, FN.Pred025 = FN.Pred025))

}







GraphToCommunity <- function(G){
  nodes = data.frame(node=as_ids(V(G)))
  properties = list(title="Random network")
  trophic.links = data.frame(resource=gsub(".*\\|", "", as_ids(E(G))), 
                             consumer=gsub("\\|.*", "", as_ids(E(G))))
  return(Community(nodes, properties, trophic.links))
}




