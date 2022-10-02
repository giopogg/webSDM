
# Cross validation function. Computes prediction from the model in K-fold CV and, if asked, also computes predictions of the fundamental niche
# INPUT:
# m: the fitted model
# K: the number of folds. can be NULL if index is specified
# index: a vector containing the fold assigned to each of the sites

trophicSDM_CV = function(m,K,partition=NULL, fundNiche=F,prob.cov=T,iter=m$iter,
                   pred_samples=m$iter,run.parallel,verbose=F){

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

    m_K = trophicSDM(Y=Y,X=X,G=m$G,env.formula=m$form.env,penal=m$penal,method=m$method,mode=m$mode,
                     family=m$family,iter=iter,
                     run.parallel = run.parallel, verbose = verbose,
                     chains=2)

    pred_K = trophicSDM_predict(m=m_K,Xnew = m$data$X[test,],
                                binary.resp=F,prob.cov=prob.cov
                                ,pred_samples=pred_samples
    )

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
