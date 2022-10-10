#' Plots the metaweb G
#' 
#' Plots the metaweb G used to fit the trophicSDM model
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @author Giovanni Poggiato
#' @return A ggnet object
#' @export
#' @examples
#'

plotG = function(tSDM){

  #########Checks
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  G = tSDM$data$G

  layout = layout_with_sugiyama(G)$layout
  rownames(layout) = tSDM$data$sp.name
  ggnet2(G, mode = layout, arrow.size = 10, node.alpha = 0.5, label=T, arrow.gap = 0.04)

}
