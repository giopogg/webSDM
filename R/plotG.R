#' Plots the metaweb G
#' 
#' Plots the metaweb G used to fit the trophicSDM model
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @author Giovanni Poggiato
#' @return A ggnet object
#' @export
#' @import ggplot2
#' @importFrom igraph layout_with_sugiyama
#' @importFrom GGally ggnet2
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 100,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' \donttest{
#' plotG(m)
#' }

plotG = function(tSDM){

  #########Checks
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )

  G = tSDM$data$G

  layout = layout_with_sugiyama(G)$layout
  rownames(layout) = tSDM$data$sp.name
  
  ggnet2(G, mode = layout, arrow.size = 10, node.alpha = 0.5, label=T, arrow.gap = 0.04)

}
