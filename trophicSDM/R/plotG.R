plotG = function(tSDM, beta = 0.01){

  #########Checks
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  G = tSDM$data$G

  layout = layout_with_sugiyama(G)$layout
  rownames(layout) = tSDM$data$sp.name
  ggnet2(G, mode = layout, arrow.size = 10, node.alpha = 0.5, label=T, arrow.gap = 0.04)

}
