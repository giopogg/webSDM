#' Plots the regression coefficients of a fitted trophicSDM model
#' 
#' @param x A trophicSDMfit object obtained with trophicSDM()
#' @param species A vector of species names to be plot. If NULL (default), all species are plotted.
#' @param ... 	additional arguments
#' @author Giovanni Poggiato
#' @examples
#' 
#' @method plot trophicSDMfit
#' @export
plot.trophicSDMfit = function(x, species = NULL, ...){

  tSDM = x
  #########Checks
  if(class(tSDM) != "trophicSDMfit") stop("tSDM is not an object of class trophicSDMfit" )


  if(!is.null(species) &
     !all(species %in% tSDM$data$sp.name)) stop("species must be either NULL or a subset of tSDM$data$sp.name")

  ##########

  # if is.null assign to the whole species pool
  if(is.null(species)) species = tSDM$data$sp.name

  plist = lapply(tSDM$model[names(tSDM$model) %in% species], function(x) plot(x, plot = F))

  nCol = 2
  nRow = ifelse(length(species) < 6, ceiling(length(species)/2), 3)
  nPage =  ceiling(length(species)/6)

  devAskNewPage(TRUE)
  for (i in 1:nPage) do.call(grid.arrange, c(plist[((i-1)*6+1):((i-1)*6+6)], ncol = nCol, nrow = nRow))
  devAskNewPage(options("device.ask.default")[[1]])

}
