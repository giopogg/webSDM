plotG = function(tSDM){

  #########Checks
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  G.mnet = build_metanet(tSDM$data$G)

  G.mnet = compute_TL(G.mnet)

  ggmetanet(metanetwork = G.mnet, beta = 0.1)

}
