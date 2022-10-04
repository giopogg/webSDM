# Plot the inferred effects

plotG_inferred = function(tSDM, level = 0.90){

  #########Checks
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  if(!is.null(tSDM$model.call$sp.formula)) stop("plotAlpha only works without composite variables")
  G = tSDM$data$G

  # Take biotic coefficients only
  all_bio_list = lapply(tSDM$model, function(x) {
    all_coef = coef(x, standardise = T, level = level)

    # In the bayesian case set to zeros links that are not significant through credible intervals
    if(x$method == "stan_glm") all_coef = apply(all_coef, 1,
                                                function(x) ifelse(x[2] <0 & x[3]>0, 0, x[1]))

    # In the frequentist case set to zeros link that are not significant through p-values
    if(x$method == "glm" & is.null(x$penal)) all_coef =  apply(all_coef, 1,
                                                               function(x) ifelse(x[2] < 1-level, 0, x[1]))

    #select only biotic coeff
    all_coef[unlist(sapply(tSDM$data$sp.name, function(x) grep(x, names(all_coef))))]
  }
  )

  all_bio = unlist(all_bio_list)

  #Assign to each edge the inferred coefficient as attribute
  idx = vector()
  for(i in 1:length(all_bio)){
    tmp = get.edge.ids(G, c(gsub("\\..*","",names(all_bio[i])),
                          gsub(".*\\.","",names(all_bio[i]))),
                     directed = F)

    idx = c(idx,tmp )

  }
  edge.attributes(G)$weight = all_bio[order(idx)]

  layout = layout_with_sugiyama(G)$layout
  rownames(layout) = tSDM$data$sp.name

  edge.color_loc = sapply(1:S, function(x) ifelse(edge.attributes(G)$weight[x]>0, "positive", "negative"))
  edge.color_loc[which(edge.attributes(G)$weight == 0)] = "non-significant"


  edge.color_loc = sapply(1:S, function(x) ifelse(edge.attributes(G)$weight[x]>0, "#CC0000", "#0000CC"))
  edge.color_loc[which(edge.attributes(G)$weight == 0)] = "grey"


  ggnet2(G, mode = layout, arrow.size = 8, node.alpha = 0.5, label=T, arrow.gap = 0.04,
         edge.label = as.character(signif(edge.attributes(G)$weight,1)), edge.color =  edge.color_loc,
         edge.alpha = ifelse(edge.color_loc == "grey", 0.5, 1),
         edge.lty = ifelse(edge.color_loc == "grey", 2, 1))



}
