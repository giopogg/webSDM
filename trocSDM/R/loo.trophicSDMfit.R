#' Computes an approximation of loo for the whole model
#' 
#' The global loo is computed by summing the loo of all the local models (since the likelihood factorises, the log-likelihood can be summed)
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @author Giovanni Poggiato
#' @method loo trophicSDMfit
#' @export
#' @examples
#' 

loo.trophicSDMfit = function(tSDM){
  
  if(class(tSDM) != "trophicSDMfit") stop("tSDM needs to be a trophicSDMfit object")
  
  if(tSDM$model.call$method != "stan_glm") stop("loo is available only for stan_glm method")
  
  return(do.call(sum,lapply(tSDM$model, function(x) loo(x$model)$estimates[1,1])))
  
}