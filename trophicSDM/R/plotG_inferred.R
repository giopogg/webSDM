# Plot the inferred effects

plotG_inferred = function(tSDM, level = 0.95){

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

  # Now plot
  G.mnet = build_metanet(tSDM$data$G)

  G.mnet = compute_TL(G.mnet)

  ggmetanet(metanetwork = G.mnet, beta = 0.1)


}
