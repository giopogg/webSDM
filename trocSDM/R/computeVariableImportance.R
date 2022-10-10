#' Computes variable importance of (groups of) variables of fitted a trophicSDM model.
#' 
#' Computes variable importance of (groups of) variables of fitted a trophicSDM model, for each species. Variable importance are computed as the standardised regression coefficients (summed across species of the same group). Standardisation is done using latent variable standardisation described in Grace et al. 2018.
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param groups A list where each element is group. Each group is specified as a vector containing species or environmental covariates names of a given group. Each element of the list (i.e. each group) has to be named.
#' @return A groups x species matrix containing variable importance for each groups of variables and each species.
#' @author Giovanni Poggiato
#' #' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @examples
#'
#' @export



computeVariableImportance = function(tSDM, groups = NULL){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class SDMfit" )

  if(is.null(groups)) {
    groups = as.list(c(colnames(tSDM$data$X)[-1], colnames(tSDM$data$Y)))
    names(groups) = c(colnames(tSDM$data$X)[-1], colnames(tSDM$data$Y))
  }

  if(is.null(names(groups))) stop("groups should be a lists with names")

  n = tSDM$data$n
  p = tSDM$data$p
  S = tSDM$data$S


  VI = matrix(NA,nrow = length(groups), ncol = S, dimnames = list(names(groups),names(tSDM$model)))

  for(j in 1:S){

    coef_temp = coef(tSDM$model[[j]], standardise = T)[-1,"estimate"]

    for(k in 1:length(groups)){

      sel = names(coef_temp) %in% groups[[k]]
      if(any(sel)) VI[k,j] = sum(abs(coef_temp[sel])) else VI[k,j] = 0

    }
  }

  return(VI)

}
