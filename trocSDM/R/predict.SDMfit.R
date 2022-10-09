#' Computes predicted values for a local model, i.e., a fitted SDMfit object
#' This is sequentially called, for each species, by the function trophicSDM.predict
#' @param object A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param newdata A matrix containing both environmental covariates and the biotic variables that the local model uses to predict the species distribution.
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param prob.cov If set to FALSE, it gives back also predicted presence-absences (which is then used by trophicSDM.predict to predict the predators).
#' @param ... 	additional arguments
#' @return A list containing for each species the predicted value at each sites. If method = "stan_glm", then each element of the list is a sites x pred_samples matrix containing the posterior predictive distribution of the species at each sites.
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @examples
#'
#' @method predict SDMfit
#' @export

predict.SDMfit=function(object, newdata, pred_samples, prob.cov,...){

  SDM = object 
  if(class(SDM) != "SDMfit") stop("you need to provide a SDMfit object")
  family=SDM$family$family
  method=SDM$method

  # If presence/absence data

  if(family %in% c("bernoulli", "binomial")){

    if(method=="stan_glm"){

      # retrieve the expected predicted probability of presence
      ## ! for the particular case of the stan_glm model with a constraint on the coefficients
      ## ! signs another stan-for-R package is used (brsm) and the code has to be adapted
      if (class(SDM$model)[1]=='brmsfit'){  # brms package
        predictions.prob=t(posterior_epred(SDM$model, newdata=as.data.frame(newdata), nsamples=pred_samples))
      }else{  # rstan package
        predictions.prob=t(posterior_epred(SDM$model, newdata=as.data.frame(newdata), draws=pred_samples))
      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov) {
        if (class(SDM$model)[1]=='brmsfit'){  # brms package
          predictions=t(posterior_predict(SDM$model, newdata=as.data.frame(newdata), nsamples=pred_samples))
        }else{  # rstan package
          predictions=t(posterior_predict(SDM$model, newdata=as.data.frame(newdata), draws=pred_samples))
        }
      }else{
        predictions = NULL
        }
    }


    if(method=="glm")  {
      if(pred_samples!=1) stop("pred_sample must be 1 if method=glm!")

      # retrieve the expected predicted probability of presence
      if(is.null(SDM$penal)){
      predictions.prob=as.matrix(predict(SDM$model,type = "response",newdata=as.data.frame(newdata)))
      }else{

        newx = model.frame(gsub(".*y","", SDM$form.all), as.data.frame(newdata))

        predictions.prob=as.matrix(predict(SDM$model,type = "response",
                                           newx= as.matrix(newx)))

      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov){
        predictions=as.matrix(sapply(as.vector(predictions.prob),FUN=function(x) rbinom(prob=x,size=1,n=1)))
      }else{
        predictions = NULL
        }
    }

    return(list(predictions.prob = predictions.prob, predictions.bin = predictions))

  }

  # If gaussian data

  if(family=="gaussian"){

    if(method=="stan_glm"){
      predictions=t(posterior_predict(SDM$model, newdata=as.data.frame(newdata),draws=pred_samples))
    }

    if(method=="glm")  {

      if(pred_samples!=1) {

        stop("pred_sample must be 1 if method=glm!")}

      predictions=as.matrix(predict(SDM$model,type = "response",newdata=as.data.frame(newdata)))

    }

    return(predictions)
  }
}
