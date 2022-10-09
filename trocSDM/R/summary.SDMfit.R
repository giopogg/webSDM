#' Prints the summary of a fitted SDMfit model
#' @param object A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param ... 	additional arguments
#' @author Giovanni Poggiato
#' @examples
#'
#' @method summary SDMfit
#' @export
summary.SDMfit = function(object, ...){
  
  SDM = object

  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )
  cat("==================================================================\n")

  model = paste0("Local SDMfit for species ", SDM$sp.name, " with ",
                 ifelse(is.null(tSDM$model.call$penal), "no", tSDM$model.call$penal), " penalty ", SDM$penal,
                 ", fitted using ", SDM$method,
                 " \n")
  cat(model)
  cat("==================================================================\n")
  cat("* Useful S3 methods\n")
  cat("    print(), coef(), plot(), predict()\n")
  cat(paste0("    $model gives the ",class(SDM$model)[1], " class object \n"))
  cat("==================================================================\n")
  summary(SDM$model)

}
