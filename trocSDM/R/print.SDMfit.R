#' Prints a SDMfit object
#' @param x A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param ... 	additional arguments
#' @author Giovanni Poggiato
#' @method print SDMfit
#' @export
#' @examples
#'
print.SDMfit = function(x, ...){

  SDM = x
  
  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  summary(SDM)
}

