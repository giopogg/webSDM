#' Gets regression coefficients from a fitted trophicSDM model.
#' @param object A trophicSDMfit object obtained with trophicSDM()
#' @param standardise Whether to standardise regression coefficients. Default to FALSE. If TRUE, coefficients are standardised using the latent variable standardisation (see Grace et al. 2018) for more details.
#' @param level The confidence level of credible intervals, only available for stan_glm method. Default to 0.95.
#' @param ... 	additional arguments
#' @return A list containing, for each species, the inferred coefficients (with credible intervals or p-values when available).
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @author Giovanni Poggiato
#' @examples
#' 
#' @method coef trophicSDMfit
#' @export
coef.trophicSDMfit = function(object, standardise = F, level = 0.95, ...){
  
  tSDM = object

  if(class(tSDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  lapply(tSDM$model, function(x) coef(x, standardise = T, level = level))

}
