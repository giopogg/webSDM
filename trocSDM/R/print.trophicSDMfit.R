#' Prints a fitted trophicSDM model
#' @param x A trophicSDMfit object obtained with trophicSDM()
#' @param ... 	additional arguments
#' @author Giovanni Poggiato
#' @method print trophicSDMfit
#' @export
#' @examples
#'

print.trophicSDMfit = function(x, ...){

  tSDM = x
  
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  summary(tSDM)
}
