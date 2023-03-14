#' Computes an approximation of loo for the whole model
#'
#' Only works if method = 'stan_glm'. The global loo is computed by summing the loo of all the local models (since the likelihood factorises, the log-likelihood can be summed)This is an implementation of the methods described in Vehtari, Gelman, and Gabry (2017) and Vehtari, Simpson, Gelman, Yao, and Gabry (2019).
#' @param x A trophicSDMfit object obtained with trophicSDM()
#' @param ... 	additional arguments
#' @return The value of the loo for the whole model
#' @author Giovanni Poggiato
#' @importFrom brms loo log_lik
#' @importFrom rstanarm loo
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' m = trophicSDM(Y,X,G, env.formula,
#'                family = binomial(link = "logit"), penal = NULL, iter = 50,
#'                mode = "prey", method = "stan_glm")
#' \donttest{brms::loo(m)}
#' @method loo trophicSDMfit
#' @export
loo.trophicSDMfit = function(x, ...){

  tSDM = x
  
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM needs to be a trophicSDMfit object")

  if(tSDM$model.call$method != "stan_glm") stop("loo is available only for stan_glm method")

  return(do.call(sum,lapply(tSDM$model, function(x) loo(x$model)$estimates[1,1])))

}
