#' Predicts with a local model
#' 
#' Computes predicted values for a local model, i.e., a fitted SDMfit object This is sequentially called, for each species, by the function trophicSDM.predict
#' @param object A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param newdata A matrix containing both environmental covariates and the biotic variables that the local model uses to predict the species distribution.
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param prob.cov Only for presence-absence data. If set to FALSE, it gives back also predicted presence-absences (which is then used by trophicSDM.predict to predict the predators).
#' @param ... 	additional arguments
#' @return A list containing for each species the predicted value at each sites. If method = "stan_glm", then each element of the list is a sites x pred_samples matrix containing the posterior predictive distribution of the species at each sites. If prob.cov = TRUE, it returns a list containing:
#' \itemize{
#' \item{predictions.prob}{Predicted probabilities of presence}
#' \item{predictions.prob}{Predicted presence-absences}
#' }
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # In order to predict non-basal species, we need to also provide
#' # the predicted occurrences of its preys. Here we compute the probability of
#' # presence of species Y4 at environemntal conditions c(0.5,0.5)
#' # when its prey Y3 is present.
#' predict(m$model$Y4, newdata = data.frame(X_1 = 0.5, X_2 = 0.5, Y3 = 1), pred_samples = 10)
#' @importFrom  rstanarm posterior_epred
#' @importFrom  brms posterior_epred
#' @importFrom  rstanarm posterior_predict
#' @importFrom  brms posterior_predict
#' @importFrom stats predict model.frame rbinom
#' @method predict SDMfit
#' @export

predict.SDMfit = function(object, newdata, pred_samples = NULL, prob.cov = TRUE,...){

  SDM = object 
  if(!inherits(SDM, "SDMfit")) stop("you need to provide a SDMfit object")
  
  if(is.null(pred_samples)){
    
    if(SDM$method=="glm") pred_samples = 1
    if(SDM$method=="stan_glm") pred_samples = round(SDM$iter/10)
    
  }
  
  family = SDM$family$family
  method = SDM$method

  # If presence/absence data

  if(family %in% c("bernoulli", "binomial")){

    if(method=="stan_glm"){

      # retrieve the expected predicted probability of presence
      ## ! for the particular case of the stan_glm model with a constraint on the coefficients
      ## ! signs another stan-for-R package is used (brsm) and the code has to be adapted
      if (inherits(SDM$model,'brmsfit')){  # brms package
        predictions.prob=t(posterior_epred(SDM$model, newdata=as.data.frame(newdata), nsamples=pred_samples))
      }else{  # rstan package
        predictions.prob=t(posterior_epred(SDM$model, newdata=as.data.frame(newdata), draws=pred_samples))
      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov) {
        if(inherits(SDM$model,'brmsfit')){  # brms package
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

  if(SDM$family$family =="gaussian"){

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
