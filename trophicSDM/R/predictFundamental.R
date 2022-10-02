predictFundamental = function(tSDM, Xnew = NULL, pred_samples = NULL, run.parallel = T,
                              verbose = F, fullPost = T){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  if(is.null(Xnew)) Xnew = tSDM$data$X

  if(is.null(pred_samples)){

    if(tSDM$model.call$method=="glm") pred_samples = 1
    if(tSDM$model.call$method=="stan_glm") pred_samples = tSDM$model.call$iter

  }

  family = tSDM$model.call$family
  n = nrow(Xnew)
  S = tSDM$data$S
  G = tSDM$data$G
  mode = tSDM$model.call$mode

  sp.prediction = as.list(vector(length=vcount(G)))
  names(sp.prediction) = tSDM$data$sp.name

  for(s in 1:S){

    sp_mod = tSDM$model[[s]]


    #fix all species to one (i.e. present)
    newdata =  cbind(Xnew,
                     data.frame(matrix(1, nrow=n, ncol=ncol(tSDM$data$Y),
                                       dimnames= list(NULL,colnames(tSDM$data$Y))))
    )

    pred_temp = predict(sp_mod,newdata = newdata,
                        pred_samples = pred_samples, prob.cov = T)

    sp.prediction[[s]] = pred_temp$predictions.prob

  }

  if(tSDM$model.call$method == "stan_glm" & !fullPost ){

    # predictions.prob
    for(j in 1:length(sp.prediction)){
      predictions.mean = rowMeans(sp.prediction[[j]])
      predictions.q975 = apply(sp.prediction[[j]], 1, quantile, 0.975)
      predictions.q025 = apply(sp.prediction[[j]], 1, quantile, 0.025)

      sp.prediction[[j]] = list(predictions.mean = predictions.mean,
                                                 predictions.q025 = predictions.q025,
                                                 predictions.q975 = predictions.q975
      )

    }
  }

  return(sp.prediction)

}
