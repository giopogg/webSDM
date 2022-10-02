# level only for stan_glm
coef.trophicSDMfit = function(tSDM, standardise = F, level = 0.95){

  if(class(tSDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  lapply(tSDM$model, function(x) coef(x, standardise = T, level = level))

}
