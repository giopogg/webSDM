
# Cross validation function. Computes prediction from the model in K-fold CV and, if asked, also computes predictions of the fundamental niche
# INPUT:
# m: the fitted model
# K: the number of folds. can be NULL if index is specified
# index: a vector containing the fold assigned to each of the sites

trophicSDM_CV = function(tSDM,K,partition=NULL, prob.cov=T, iter = NULL,
                         pred_samples = NULL, run.parallel = T, verbose=F){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM needs to be a trophicSDMfit object")

  # set MCMC parameters
  if(tSDM$model.call$method == "glm"){ iter = 1; pred_samples = 1
  }else{
    if(is.null(iter)) iter = tSDM$model.call$iter
    if(is.null(pred_samples)) pred_samples = tSDM$model.call$iter
  }

  n = tSDM$data$n
  S = tSDM$data$S
  G = tSDM$data$G

  # Create partition (if needed)
  if(!is.null(partition)){
    if(length(partition) != n) stop("partition must be a vector of length n (the number of sites")
  }else{
    partition <- sample(1:K,size=n,replace=TRUE,prob=rep(n/K,K))
    }


  preds = array(dim=c(n,pred_samples,S))

  for(i in 1:K){

    train = which(partition != i)
    test = which(partition == i)
    Y = tSDM$data$Y[train,]
    X = tSDM$data$X[train,]

    m_K = trophicSDM(Y = Y, X = X, G = tSDM$data$G, env.formula = tSDM$model.call$form.env,
                     penal = tSDM$model.call$penal, method = tSDM$model.call$method, mode = tSDM$model.call$mode,
                     family = tSDM$model.call$family, iter = iter,
                     run.parallel = run.parallel, verbose = verbose,
                     chains=2)

    pred_K = predict(m_K, Xnew = tSDM$data$X[test,],
                     prob.cov = prob.cov, pred_samples = pred_samples)

    preds[test,,] = abind(pred_K,along=3)

    print(paste0("Fold ", i, " out of ", K,"\n"))

  }

  if(tSDM$model.call$method == "glm"){

    meanPred = preds[,1,]

    colnames( meanPred ) = names(pred_K)

  }
  if(tSDM$model.call$method == "stan_glm"){

  meanPred = apply(preds,mean,MARGIN = c(1,3))
  Pred975 = apply(preds, quantile, MARGIN=c(1,3),0.975)
  Pred025 = apply(preds, quantile, MARGIN=c(1,3),0.025)

  colnames( meanPred ) = colnames(Pred975) = colnames(Pred025) = names(pred_K$sp.prediction)
  }


  # Compute Joint TSS and AUC
  # if(tSDM$model.call$family$family == "binomial"){
  # auc = tss = vector(length = S)
  #
  # eval = mclapply(1:S,function(x){
  #   eval = dismo::evaluate(p = meanPred[which(tSDM$data$Y[,colnames(meanPred)[x]]==1),x],
  #                          a = meanPred[which(tSDM$data$Y[,colnames(meanPred)[x]]==0),x] )
  #
  #   data.frame(auc = eval@auc, tss = max(eval@TPR+eval@TNR-1))
  # })
  #
  # CVmetrics = do.call(rbind,eval)
  # CVmetrics$species = colnames(meanPred)
  # }

  list(meanPred = meanPred, Pred975 = Pred975, Pred025 = Pred025, partition = partition)

}
