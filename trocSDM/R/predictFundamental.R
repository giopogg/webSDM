#' Computes predicted values for a local model, i.e., a fitted SDMfit object
#' This is sequentially called, for each species, by the function trophicSDM.predict
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param Xnew a matrix specifying the environmental covariates for the predictions to be made. If NULL (default), predictions are done on the training dataset (e.g. by setting Xnew = tSDM$data$X).
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param run.parallel Whether to use parallelise code when possible. Can speed up computation time.
#' @param verbose Whether to print advances of the algorithm.
#' @param fullPost Optional parameter for stan_glm only. Whether to give back the full posterior predictive distribution (default, fillPost = TRUE) or just the posterior mean, and 2.5% and 97.5% quantiles.
#' @return A list containing for each species the predicted value at each sites. If method = "stan_glm", then each element of the list is a sites x pred_samples matrix containing the posterior predictive distribution of the species at each sites.
#' @export
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @examples
#'



predictFundamental = function(tSDM, Xnew = NULL, pred_samples = NULL, run.parallel = T,
                              verbose = F, fullPost = T){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  if(is.null(Xnew)) Xnew = tSDM$data$X

  if(is.null(pred_samples)){

    if(tSDM$model.call$method=="glm") pred_samples = 1
    if(tSDM$model.call$method=="stan_glm") pred_samples = tSDM$model.call$iter/10

  }

  family = tSDM$model.call$family
  n = nrow(Xnew)
  S = tSDM$data$S
  G = tSDM$data$G
  mode = tSDM$model.call$mode

  sp.prediction = as.list(vector(length=vcount(G)))
  names(sp.prediction) = names(tSDM$model)

  for(s in 1:S){

    sp_mod = tSDM$model[[s]]


    #fix all species to oneif species are modeled as a function of their preys (i.e., all preys are present) , elsewhere they are fixed to 0 (i.e. all predators are absent).

    if(mode == "prey"){
    newdata =  cbind(Xnew,
                     data.frame(matrix(1, nrow=n, ncol=ncol(tSDM$data$Y),
                                       dimnames= list(NULL,colnames(tSDM$data$Y)))))
    }else{
      newdata =  cbind(Xnew,
                       data.frame(matrix(0, nrow=n, ncol=ncol(tSDM$data$Y),
                                         dimnames= list(NULL,colnames(tSDM$data$Y)))))
    }

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
