##########################################################################################################
############ Functions trophic-SDMs
##########################################################################################################



#### Notes:

# stanarm, stan_glm


# https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
# https://cran.r-project.org/web/packages/rstanarm/vignettes/binomial.html
# http://mc-stan.org/rstanarm/articles/priors.html
# http://mc-stan.org/rstanarm/reference/priors.html
# https://avehtari.github.io/modelselection/diabetes.html
# https://mc-stan.org/rstanarm/reference/stan_glm.html

# Sparsification:

#https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-information-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.full
#http://mc-stan.org/rstanarm/reference/priors.html



#required libraries

library(igraph)
library(dplyr)
library(rstanarm)
library(arm)
library(glmnet)
library(formula.tools)
library(abind)



###########################################################################################################################
############ Fit functions
############################################################################################################################


### trophicSDM : the function to fit a trophicSDM. Given the interaction graph, the community matrix
# and the environmental layer, it fits the cascade of SDM. The environmental formulas, the SDM method,
# the penalisation and eventually the form of the composite functions can be chosen by the user

## Parameters:
# Y : Community data must be a sites x species matrix
# X : Environmental data must be a sites x predictor matrix
# G : Trophic network, igraph
# formulas : list (names of each element are the name of vertices) of one-hand formulas whose names are the same as in X. If NULL, all variables in X are included, linearly.
# method: which SDM method to use (available for now: glm (frequentist), stan_glm (full bayesian MCMC), bayesglm (variational approx of the posterior))
# penalisation: how to penalise the likelihood to reduce the dimension. Each SDM method has its own penalisation type, some SDM method have no available penalisation approach. For now, we have "horshoe" for stan_glm and "elasticnet" for glm
# family: the family parameter of the glm function (see glm). family=gaussian for gaussian data or family=binomial(logit) or binomial(probit) for presence-absence data
# iter : the number of MCMC samples (if bayesian approach)

trophicSDM = function(Y,X,G,formulas=NULL,penal=NULL,method="stan_glm",family,iter=1000){
  
  #checks & errors
  if(!is_igraph(G)) stop("G is not an igraph object")
  if(!is.dag(G)) stop("G is not acyclic, DAG legacy is not over yet...")
  if(!method %in% c("glm","stan_glm","bayesglm")) stop("the selected method has not yet been implemented or it has been misspelled")
  if(!(is.null(penal) || penal %in% c("horshoe","elasticnet"))) stop("the selected penalisationhas not yet been implemented or it has been misspelled")
  
  #topological sorting of the graph. mode="in" as we assume top-down control, "out otherwise
  sortedV= topo_sort(G,mode="in")
  
  #initialize empty lists of models
  m=form.all=as.list(vector(length=vcount(G)))
  names(m)=names(form.all)=names(sortedV)
  
  #core part: loop on the species to fit their distribution
  for(j in 1:vcount(G)){ 
    
    #call a function that does an SDM with j as focal species, automatically finds covariates from
    # G and formula.foc
    temp.mod = SDMfit(focal=names(sortedV[j]),Y,X,G,
                        formula.foc=ifelse(!is.null(formulas),as.character(formulas[[names(sortedV[j])]]),NULL),
                        method=method,penal=penal,family=family,iter=iter)
    
    #assign models and formulas (that now includes biotic variables too)
    m[[names(sortedV[j])]]=temp.mod$m
    form.all[[names(sortedV[j])]]=temp.mod$form.all
  }
  #return values
  list(model=m,form.all=form.all,data=list(X=X,Y=Y),method=method,penal=penal,iter=iter,family=family,G=G)
}




### SDMfit : the core function that fit a (classic) SDM on the focal species, with its preys (predators if top-down approach)
# and the environmental covariates as predictors.

## Parameters:
# focal : the vertex of the graph (i.e. the species) to be modelled
# other parameters: see trophicSDM


# nice stuff about formulas https://www.datacamp.com/community/tutorials/r-formula-tutorial

SDMfit=function(focal,Y,X,G,formula.foc,method="bayesglm",family=NULL,penal=NULL,iter=iter){
  
  #build environmental part of the formula
  if(is.null(formula.foc)){ 
    #all env variables included
    form.env= as.formula(paste("y~",paste(colnames(X),collapse="+")))
  }else{form.env = as.formula(paste("y",as.character(formula.foc)))}

  #build final formula
  preys=names(neighbors(G,focal,mode=c("out")))
  if(length(preys)>0){
  form.all=as.formula(paste(as.character(form.env),paste(preys,collapse = "+"),sep="+"))
  }else{ 
    #the focal species is basal
    form.all=form.env
  }
  
  #build the data matrix to be given to the fit functions
  data = data.frame(X,Y)
  data = dplyr::rename(data,y=focal)  #rename the focal column to be called "y"
  
  
  # Fit models depending on the chosen method and penalisation
  
  ### GLM
  if(method=="glm"){
    if(is.null(penal)){
      
      m= glm(form.all,data=data,family=family)
    }else{
    if(penal=="horshoe") stop("Horshoe penalisation is not possible with glm")
    
    
    if(penal=="elasticnet"){
      y=as.matrix(model.frame(form.all,data=data)$y)
      x=as.matrix(model.frame(form.all,data=data)[,-1])
    
      m = glmnet(y=y,x=x,family=family$family,
                 lambda=cv.glmnet(y=y,x=x,family=family$family)$lambda.1se)
    }
      
    }
  }
  ### bayesglm
  
  if(method=="bayesglm"){
    if(!is.null(penal)) {stop("Penalisation is not possible with bayesglm")
    }else{
    m = bayesglm(form.all,data=data,family=family)
    }
  }
  
  ### stan_glm
  if(method=="stan_glm"){
    
  
  if(is.null(penal)){

           m = stan_glm(form.all, data = data,
                  family = family,iter=iter,
                  prior=normal(0,100))
        }else{
  if(penal=="elasticnet")  stop("Elastic net is not stan_glm")
  if(penal=="horshoe"){

      HS_prior = hs(df = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)
    
        m = stan_glm(form.all, data = data,
                    family = family,
                    prior = HS_prior,iter=iter)
      }
  }
  }

  #create matrix of coefficient with postmean, post quantiles and whole samples
  
  #return
  list(m=m,form.all=form.all)
}


###########################################################################################################################
############ Predict functions
############################################################################################################################



### trophicSDM_predict : the function that predicts the output of the trophicSDM model, eventually
# depending on a new set of environmental covariates. Several parameters are available depending ont the users' needs.

## Parameters:
# m : the fitted model
# Xnew : Environmental data must be a sites x predictor matrix. If NULL, the fitted environemntal predictors are used (i.e. in-sample prediction)
# binary.resp : only useful with presence-absence variable. TRUE if we also want binary-output, FALSE if only probabilities are given.
# prob.cov : TRUE if the predicted probabilities of presence of the preys are used to predict species. FALSE if the randomly (depending of predicted probabilities of course) generated presence-absence are used
# pred.sample : the total number of sample from the predictive distribution of each species at each site
# error_prop_sample : fine tuning parameter for error propagation. If 1, for each of the predictive samples of the preys, we draw one sample of the predicted distribution for the focal species. 
#                     If greater than one, we sample pred.sample/error_prop_sample of the predictive distribution of each species (N.B. carefully to guarantee that the samples are consistent across them), and then for each of the selected samples,
#                     we draw error_prop_sample from the predictive distribution of the focal species, for each selected sample of the predictors.

trophicSDM_predict=function(m,Xnew,binary.resp=NULL,prob.cov=NULL,pred_samples=1000,error_prop_sample=10){
  
  #checks & errors
  if(m$method=="glm" & (pred_samples!= 1 | error_prop_sample != 1)){stop("glm requires pred_sample and error_prop_sample both =1")}
  sortedV= topo_sort(G,mode="in")
  
  if(is.null(Xnew)) Xnew=m$data$X
  
  family=m$family
  n=nrow(Xnew)
  S=ncol(m$data$Y)
  G=m$G
  sp.prediction=as.list(vector(length=vcount(G)))
  names(sp.prediction)=names(sortedV)
  
  if(family$family=="binomial"){
    
    if(is.null(binary.resp))  message("Please specify the parameter binary.resp. TRUE if you also want binary predictions, FALSE otherwise")
    if(is.null(prob.cov))  message("Please specify the parameter prob.cov. TRUE if you want to use probabilities instead of presence-absence as a proxy for other species predictions. FALSE otheriwse")
    
  # core loop on species (as in trophicSDM)
  for(j in 1:vcount(G)){ 
    
    #check if basal or not
    basal=ifelse(length(neighbors(G,v=sortedV[j],mode="out"))==0,TRUE, FALSE)
    
    if(basal){
    
    #if basal, then just call a classic predict function depending on the (eventually new) environmental predictors
    sp.prediction[[names(sortedV[j])]]=SDMpredict(focal=names(sortedV[j]),m=m,newdata=Xnew,pred_samples=pred_samples,binary.resp=binary.resp,prob.cov=prob.cov)
    
    }else{
      
      # otherwise, select the predictive samples from the preys (that have already been predicted!!)
      neigh.sp = neighbors(G,v=sortedV[j],mode="out")
      #create newdata: will include the (eventually new) environmental predictors and the predicted distribution of the preys
      # it's an array where we have the predictive samples of the preys on the third dimension.
      # we will therefore apply SDMpredict to each of the matrices (i.e. first two dimensions) in order to correctly propagate the error
      newdata=array(data=NA,dim=c(n,(ncol(Xnew)+length(neigh.sp)),pred_samples/error_prop_sample))
      colnames(newdata)=c(colnames(Xnew),names(neigh.sp))
      
      
      #fill the abiotic variables
      newdata[,1:ncol(Xnew),]=as.matrix(Xnew)
      
      # fill the biotic part of new data
      for(k in 1:length(neigh.sp)){
      if(prob.cov){
        
        # if prob.cov= TRUE, use the species predicted probabilities of presence
        newdata[,ncol(Xnew)+k,]=sp.prediction[[names(neigh.sp[k])]]$predictions.prob[,seq(from=1,by=error_prop_sample,to=pred_samples)]#here select the random samples!
        }else{
          
        # if prob.cov= TRUE, use the species predicted probabilities of presence
        newdata[,ncol(Xnew)+k,]=sp.prediction[[names(neigh.sp[k])]]$predictions.bin[,seq(from=1,by=error_prop_sample,to=pred_samples)]#here select the random samples!
      }
        
      }
      # apply the function SDMpredict to each layer of the array (MARGIN=3)
      pred.array=apply(newdata,MARGIN=3,FUN = function(x){SDMpredict(focal=names(sortedV[j]),m=m,newdata=x,pred_samples=error_prop_sample,binary.resp=binary.resp,prob.cov=prob.cov)})
      #unlist and format
      sp.prediction[[names(sortedV[j])]]=list(predictions.prob=NULL,predictions.bin=NULL)
      sp.prediction[[names(sortedV[j])]]$predictions.prob = do.call(cbind,lapply(pred.array,FUN=function(x) x$predictions.prob))
      if(!prob.cov | binary.resp) {

      sp.prediction[[names(sortedV[j])]]$predictions.bin = do.call(cbind,lapply(pred.array,FUN=function(x) x$predictions.bin))
              }
    }
  }
  if(!binary.resp){
    #if response=FALSE, remove the binary predictions (that were eventually needed for computations if prob.cov=FALSE)
    for(j in 1:length(sp.prediction)){
      sp.prediction[[j]]$predictions.bin=NULL
    }
  }
  return(list(sp.prediction=sp.prediction,binary.resp=F,pred_samples=1000,error_prop_sample=10,prob.cov=FALSE))
  
  }
  
  #Just the same of above but predictions are directly given (no need to choose between prob or not)
  if(family$family=="gaussian"){
    
    for(j in 1:vcount(G)){ #eventually modify with apply
      #check if basal or not
      basal=ifelse(length(neighbors(G,v=sortedV[j],mode="out"))==0,TRUE, FALSE)
      
      if(basal){
        sp.prediction[[names(sortedV[j])]]=SDMpredict(focal=names(sortedV[j]),m=m,newdata=Xnew,pred_samples=pred_samples,binary.resp=binary.resp,prob.cov=prob.cov)
      }else{
        neigh.sp = neighbors(G,v=sortedV[j],mode="out")
        newdata=array(data=NA,dim=c(n,(ncol(Xnew)+length(neigh.sp)),pred_samples/error_prop_sample))
        colnames(newdata)=c(colnames(Xnew),names(neigh.sp))
        
        #fill the abiotic variables
        newdata[,1:ncol(Xnew),]=as.matrix(Xnew)
        
        for(k in 1:length(neigh.sp)){
            newdata[,ncol(Xnew)+k,]=sp.prediction[[names(neigh.sp[k])]][,seq(from=1,by=error_prop_sample,to=pred_samples)]#here select the random samples!
          }
          
        pred.array=apply(newdata,MARGIN=3,FUN = function(x){list(SDMpredict(focal=names(sortedV[j]),m=m,newdata=x,pred_samples=error_prop_sample,binary.resp=binary.resp,prob.cov=prob.cov))})
        sp.prediction[[names(sortedV[j])]] = do.call(cbind,lapply(pred.array,FUN=function(x) x[[1]]))
      }
    }
    
    return(list(sp.prediction=sp.prediction,pred_samples=1000,error_prop_sample=10))
    
  }
  
}



### SDMpredict : the function that predicts the distribution of a single species. 
#                Simply recalls the "predict" function of the right model

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
  if(family=="binomial"){
    #initialize vectors
     predictions.prob=matrix(data=NA,nrow=n,ncol=pred_samples)
     if(!prob.cov | binary.resp) {predictions=matrix(data=NA,nrow=n,ncol=pred_samples)}

  
  if(method=="stan_glm"){
       predictions.prob=t(posterior_epred(m$model[[focal]], newdata=as.data.frame(newdata),draws=pred_samples))
       if(!prob.cov | binary.resp) {
       predictions=t(posterior_predict(m$model[[focal]], newdata=as.data.frame(newdata),draws=pred_samples))
       }
  }
    
  if(method=="bayesglm"){
    
       coef.sim=coef(sim(m$model[[focal]],n.sims=pred_samples))
       
       data.temp=as.matrix(cbind(1,model.frame(m$form.all[[focal]][-2], data=as.data.frame(newdata))))
       predictions.prob=apply(coef.sim,FUN=function(x) 1/(1+exp(-data.temp%*%x)),MARGIN=1)
       if(binary.resp & !prob.cov) {
         predictions=matrix(sapply(as.vector(predictions.prob),FUN=function(x) rbinom(prob=x,size=1,n=1)),nrow=n)
       }
  }
  
  if(method=="glm")  {
      if(pred_samples!=1) {
        stop("pred_sample must be 1 if method=glm!")}
        predictions.prob=as.matrix(predict(m$model[[focal]],type = "response",newdata=as.data.frame(newdata)))
                         
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

