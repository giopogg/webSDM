#' Computes variable importance of (groups of) variables of fitted a trophicSDM model.
#' 
#' Computes variable importance of (groups of) variables of fitted a trophicSDM model, for each species. Variable importance are computed as the standardised regression coefficients (summed across species of the same group). Standardisation is done using latent variable standardisation described in Grace et al. 2018.
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param groups A list where each element is group. Each group is specified as a vector containing species or environmental covariates names of a given group. Each element of the list (i.e. each group) has to be named.
#' @return A groups x species matrix containing variable importance for each groups of variables and each species.
#' @author Giovanni Poggiato
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @importFrom stats coef
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 100,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' #Compute the importance of each variable
#' computeVariableImportance(m)
#' #Compute the importance of three different set of variables
#' computeVariableImportance(m, groups =list("X" = c("X_1","X_2"), 
#'                                              "Ybasal" = c("Y1","Y2","Y3"),
#'                                              "Ypredator"= c("Y4", "Y5", "Y6")))
#' @export

computeVariableImportance = function(tSDM, groups = NULL){

  if(!is.null(tSDM$model.call$penal)) {if(tSDM$model.call$penal == "coeff.signs"){stop("This function is not available for coeff.signs penalisation")}}
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class SDMfit" )
  if(!is.null(tSDM$model.call$sp.formula)) warning("If you use composite variables, you should group together species that belong to the same composite variable. For example, if sp.formula = 'richness' and sp.partition = NULL, you should put all species in the same group in the argument 'groups'. If you define a partition of species in sp.partition, then species in the same group in sp.partition should put all species in the same group in the argument 'groups'")
  
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

      sel = unique(unlist(sapply(groups[[k]], grep, names(coef_temp))))
      if(length(sel) > 0) VI[k,j] = sum(abs(coef_temp[sel])) else VI[k,j] = 0

    }
  }

  return(VI)

}
