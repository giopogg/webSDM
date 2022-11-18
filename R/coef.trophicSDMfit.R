#' Gets regression coefficients from a fitted trophicSDM model.
#'
#' Gets regression coefficients (eventually standardised) of a fitted trophicSDM. p-values or credible intervals are returned when available. 
#' @param object A trophicSDMfit object obtained with trophicSDM()
#' @param standardise Whether to standardise regression coefficients. Default to FALSE. If TRUE, coefficients are standardised using the latent variable standardisation (see Grace et al. 2018) for more details.
#' @param level The confidence level of credible intervals, only available for stan_glm method. Default to 0.95.
#' @param ... 	additional arguments
#' @return A list containing, for each species, the inferred coefficients (with credible intervals or p-values when available).
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @author Giovanni Poggiato
#' @importFrom stats coef var
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # unstandardised regression coefficients
#' coef(m)
#' #standardised regression coefficients with 90% credible intervals
#' coef(m, standardised = TRUE, level = 0.9)
#' # Run the same model using glm as fitting method
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' # Now we have p-values instead of credible intervals
#' coef(m)
#' 
#' # Notice that unstandardised coefficients are always accessible
#' # in the fitted model:
#' m$coef
#' @method coef trophicSDMfit
#' @export

coef.trophicSDMfit = function(object, standardise = FALSE, level = 0.95, ...){
  
  tSDM = object

  if(!inherits(tSDM, "trophicSDMfit")) stop("object is not of class trophicSDMfit" )
  if(!is.null(tSDM$model.call$penal)){if(tSDM$model.call$penal == "coeff.signs"){stop("This function is not available for coeff.signs penalisation")}}
  
  lapply(tSDM$model, function(x) coef(x, standardise = standardise, level = level))

}
