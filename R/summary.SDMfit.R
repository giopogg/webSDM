#' Summary of a fitted SDMfit model
#' @param object A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param ... 	additional arguments
#' @return Prints a summary of the local SDM
#' @author Giovanni Poggiato
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' summary(m$model$Y1)
#' @method summary SDMfit
#' @export
summary.SDMfit = function(object, ...){
  
  SDM = object

  if(!inherits(SDM, "SDMfit")) stop("SDM is not an object of class SDMfit" )
  cat("================================================================== \n")

  model = paste0("Local SDMfit for species ", SDM$sp.name, " with ",
                 ifelse(is.null(SDM$penal), "no", SDM$penal), " penalty ", SDM$penal,
                 ", fitted using ", SDM$method,
                 " \n")
  cat(model)
  cat("================================================================== \n")
  cat("* Useful S3 methods\n")
  cat("    print(), coef(), plot(), predict()\n")
  cat(paste0("    $model gives the ",class(SDM$model)[1], " class object \n"))
  cat("================================================================== \n")
  summary(SDM$model)
  
  #Just to fix pkgdown
  invisible(object)

}
