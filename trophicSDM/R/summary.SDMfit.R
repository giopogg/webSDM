summary.SDMfit = function(SDM){

  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )


  cat("==================================================================\n")

  model = paste0("The local SDM fit for species ", SDM$sp.name, " with penalty ", SDM$penal,
                 ", fitted using ", SDM$method,
                 " \n")
  cat(model)
  cat("==================================================================\n")
  cat("* Useful S3 methods\n")
  cat("    print(), coef(), plot(), predict()\n")
  cat(paste0("    $model gives the ",class(SDM$model)[1], " class object \n"))
  cat("==================================================================\n")
  summary(SDM$model)

}
