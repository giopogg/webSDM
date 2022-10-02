summary.trophicSDMfit = function(tSDM){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

    model = paste0("A trophicSDM fit with penalty ", tSDM$model.call$penal,
    ", fitted using ", tSDM$model.call$method,
    " \n with a ", ifelse(tSDM$model.call$mode == "prey", "bottom-up", "top-down"), " approach \n",
    " \n Number of species : ", tSDM$data$S,
    " \n Number of links : ", length(igraph::E(tSDM$data$G)),
    " \n")
    cat("==================================================================\n")
    cat(model)
    cat("==================================================================\n")
    cat("* Useful fields\n")
    cat("    $coef \n")
    cat("* Useful S3 methods\n")
    cat("    print(), coef(), plot(), predict(), evaluateModelFit() \n")
    cat("    predictFundamental(), plotG(), plotVP(), plotBiotic() \n")
    cat("* Local models (i.e. single species SDM) can be accessed through \n")
    cat("    $model\n")
}
