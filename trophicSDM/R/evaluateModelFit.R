evaluateModelFit = function(tSDM, Ynew = NULL, Ypredicted = NULL){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  if(is.null(Ynew)){

    message("You did not provide Ynew, the observed species distribution Y is used as default.")
    Ynew = tSDM$data$Y
  }

  if(is.null(Ypredicted)) {

    message("You did not provide Ypredicted, species predictions are obtained using predict()")

    if(tSDM$model.call$method == "glm") pred_samples =1
    if(tSDM$model.call$method == "stan_glm") pred_samples = tSDM$model.call$iter/10


    Ypredicted = predict(tSDM, pred_samples = pred_samples, fullPost = F)

    if(tSDM$model.call$family$family == "binomial"){
      if(tSDM$model.call$method == "glm") {
        Ypredicted = do.call(cbind, Ypredicted)
      }
      if(tSDM$model.call$method == "stan_glm"){
        Ypredicted = do.call(cbind, lapply(Ypredicted, function(x) x$predictions.mean))
      }

    }
    if(tSDM$model.call$family$family == "gaussian"){
      if(tSDM$model.call$method == "glm") {
        Ypredicted = do.call(cbind, lapply(Ypredicted))
      }
      if(tSDM$model.call$method == "stan_glm"){
        Ypredicted = do.call(cbind, lapply(Ypredicted, function(x) x$predictions.mean))
      }
    }
  }

  if(!all(colnames(Ypredicted) %in% tSDM$data$sp.name) | !all(colnames(Ynew) %in% tSDM$data$sp.name)){
    stop("colnames of Ynew and Ypredicted must be the same of tSDM$data$sp.name)")
  }

  S = tSDM$data$S

  # Compute Joint TSS and AUC
  if(tSDM$model.call$family$family == "binomial"){
    auc = tss = vector(length = S)

    eval = mclapply(1:S,function(x){
      eval = dismo::evaluate(p = Ypredicted[which(Ynew[,colnames(Ypredicted)[x]]==1),x],
                             a = Ypredicted[which(Ynew[,colnames(Ypredicted)[x]]==0),x] )

      data.frame(auc = eval@auc, tss = max(eval@TPR+eval@TNR-1))
    })

    metrics = do.call(rbind,eval)
    metrics$species = tSDM$data$sp.name
  }

  if(tSDM$model.call$family$family == "gaussian"){

    R2 = vector(length = S)

    eval = sapply(1:S, function(s) cor(Ynew[,s], Ypredicted[,s])^2)
    metrics = data.frame(R2 = eval, species = tSDM$data$sp.name)
  }
  return(metrics)

}
