#' Predicts species potential niche
#' 
#' Computes predicted values of the potential niches of species from the fitted trophicSDMfit model at environmental conditions specified by \code{Xnew}. Predictions are obtained by setting preys to present when mode = "prey" or setting predators to absent when mode = "predator".
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param Xnew a matrix specifying the environmental covariates for the predictions to be made. If NULL (default), predictions are done on the training dataset (e.g. by setting Xnew = tSDM$data$X).
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param verbose Whether to print advances of the algorithm.
#' @param fullPost Optional parameter for stan_glm only. Whether to give back the full posterior predictive distribution (default, fullPost = TRUE) or just the posterior mean, and 2.5% and 97.5% quantiles.
#' @return A list containing for each species the predicted value at each sites. If method = "stan_glm", then each element of the list is a sites x pred_samples matrix containing the posterior predictive distribution of the species at each sites.
#' @export
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @importFrom  igraph V decompose neighbors vcount
#' @importFrom stats quantile
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 100,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # Obtain 100 draws from the posterior predictive distribution of species potential niche
#' # (pred_samples = 50)
#' # Since we don't specify Xnew, the function sets Xnew = X by default
#' Ypred = predictPotential(m, fullPost = TRUE, pred_samples = 50)
#' # We can ask the function to only give back posterior mean and 95% credible intervals with
#' # fullPost = FALSE
#' \donttest{
#' Ypred = predictPotential(m, fullPost = FALSE, pred_samples = 50)
#' }
#' #' We can now evaluate species probabilities of presence for the enviromental
#' # conditions c(0.5, 0.5)
#' predictPotential(m, Xnew = data.frame(X_1 = 0.5, X_2 = 0.5), pred_samples = 50)
#' 
#' # If we fit the model using in a frequentist  way (e.g. glm)
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' # We are obliged to set pred_samples = 1 
#' # (this is done by default if pred_samples is not provided)
#' # In the frequentist case, fullPost is useless.
#' Ypred = predictPotential(m, pred_samples = 1)

predictPotential = function(tSDM, Xnew = NULL, pred_samples = NULL,
                              verbose = FALSE, fullPost = TRUE){

  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )

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
                        pred_samples = pred_samples, prob.cov = TRUE)

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
