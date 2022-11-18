#' Prints a fitted trophicSDM model
#' 
#' @param x A trophicSDMfit object obtained with trophicSDM()
#' @param ... 	additional arguments
#' @return Prints a summary of the fitted trophic SDM
#' @author Giovanni Poggiato
#' @method print trophicSDMfit
#' @export
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")

print.trophicSDMfit = function(x, ...){

  tSDM = x
  
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )

  summary(tSDM)
  
  #Just to fix pkgdown
  invisible(x)
}
