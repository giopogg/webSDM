#' Prints a SDMfit object
#' @param x A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param ... 	additional arguments
#' @return Prints a summary of the local SDM
#' @author Giovanni Poggiato
#' @method print SDMfit
#' @export
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 100,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' m$model$Y1
print.SDMfit = function(x, ...){

  SDM = x
  
  if(!inherits(SDM, "SDMfit")) stop("SDM is not an object of class SDMfit" )

  summary(SDM)
  
  #Just to fix pkgdown
  invisible(x)
}

