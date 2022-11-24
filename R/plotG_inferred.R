#' Plot the metaweb G according to the inferred coefficients
#' 
#' Plot the metaweb G with links colored accordingly to the inferred prey-predator regression coefficients of a fitted trophicSDM model. Plots the metaweb G, where each predator-prey link is colored according to whether the related regression coefficient if inferred as positive (in red), negative (in blue) or non-significant (dashed grey line) according to the confidence level specified in "level". Estimates of the significant standardised regression coefficients are pasted on the links. Only works if species are modeled as a function of their preys or predators without composite variables (i.e., the function only works if tSDM is fitted with sp.formula = NULL and sp.partition = NULL)
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param level The confidence level used to decide whether regression coefficients are non-significant or not. Default to 0.9.
#' @return A ggnet object
#' @author Giovanni Poggiato
#' @import ggplot2
#' @importFrom igraph E get.edge.ids edge_attr edge_attr<- layout_with_sugiyama
#' @importFrom GGally ggnet2
#' @importFrom stats coef
#' @export
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' \donttest{
#' plotG_inferred(m)
#' }

plotG_inferred = function(tSDM, level = 0.90){

  #########Checks
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )
  if(!is.null(tSDM$model.call$penal)){ if(tSDM$model.call$penal == "coeff.signs"){stop("This function is not available for coeff.signs oenalisation")}}
  
  if(!is.null(tSDM$model.call$sp.formula)) stop("plotG_inferred only works without composite variables")
  G = tSDM$data$G
  S = tSDM$data$S
  # Take biotic coefficients only
  all_bio_list = lapply(tSDM$model, function(x) {
    all_coef = coef(x, standardise = T, level = level)

    # In the bayesian case set to zeros links that are not significant through credible intervals
    if(x$method == "stan_glm") {all_coef = apply(all_coef, 1,
                                                function(x) ifelse(x[2] <0 & x[3]>0, 0, x[1]))
    }else{
    # In the frequentist case set to zeros link that are not significant through p-values
    if(x$method == "glm" & is.null(x$penal)) {all_coef =  apply(all_coef, 1,
                                                               function(x) ifelse(x[2] < 1-level, x[1], 0))
    }else{ all_coef = all_coef[,"estimate"]}
    }
    #select only biotic coeff
    all_coef[unlist(sapply(tSDM$data$sp.name, function(x) grep(x, names(all_coef))))]
    
  }
  )

  all_bio = unlist(all_bio_list)

  #Assign to each edge the inferred coefficient as attribute
  idx = vector()
  for(i in 1:length(all_bio)){
    tmp = igraph::get.edge.ids(G, c(gsub("\\..*","",names(all_bio[i])),
                          gsub(".*\\.","",names(all_bio[i]))),
                     directed = F)

    idx = c(idx,tmp )

  }
  edge_attr(G)$weight = all_bio[order(idx)]

  layout = layout_with_sugiyama(G)$layout
  rownames(layout) = tSDM$data$sp.name

  edge.color_loc = sapply(1:length(igraph::E(G)), function(x) ifelse(edge_attr(G)$weight[x]>0, "#CC0000", "#0000CC"))
  edge.color_loc[which(edge_attr(G)$weight == 0)] = "grey"
  
  ggnet2(G, mode = layout, arrow.size = 8, node.alpha = 0.5, label=T, arrow.gap = 0.04,
         edge.label = as.character(signif(edge_attr(G)$weight,1)), edge.color =  edge.color_loc,
         edge.alpha = ifelse(edge.color_loc == "grey", 0.5, 1),
         edge.lty = ifelse(edge.color_loc == "grey", 2, 1))

}
