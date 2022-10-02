print.trophicSDMfit = function(tSDM){

  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )

  summary(tSDM)
}
