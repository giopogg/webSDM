#' Plots the regression coefficients of a fitted trophicSDM model
#' 
#' Plots the regression coefficients of a fitted trophicSDM model. A subset of species to be plotted can be specified in the parameter\code{species}.
#' @param x A trophicSDMfit object obtained with trophicSDM()
#' @param species A vector of species names to be plot. If NULL (default), all species are plotted.
#' @param ... additional arguments
#' @return A plot of the regression coefficients of the fitted tropic SDM
#' @author Giovanni Poggiato
#' @import ggplot2
#' @importFrom  grDevices devAskNewPage
#' @importFrom gridExtra grid.arrange
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # Plot just the first three species
#' \donttest{
#' plot(m, species = c("Y1","Y2","Y3"))
#' }
#' # If species = NULL (default), all species are plotted.
#' @method plot trophicSDMfit
#' @export

plot.trophicSDMfit = function(x, species = NULL, ...){
  
  tSDM = x
  #########Checks
  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )
  
  if(!is.null(species) &
     !all(species %in% tSDM$data$sp.name)) stop("species must be either NULL or a subset of tSDM$data$sp.name")
  
  ##########
  
  # if is.null assign to the whole species pool
  if(is.null(species)) species = tSDM$data$sp.name
  
  plist = lapply(tSDM$model[names(tSDM$model) %in% species], function(x) plot(x))
  
  nCol = 2
  nRow = ifelse(length(species) < 6, ceiling(length(species)/2), 3)
  nPage =  ceiling(length(species)/6)
  
  if(tSDM$data$S>6){
  devAskNewPage(TRUE)
  }
  for (i in 1:nPage) do.call(grid.arrange,
                             c(plist[((i-1)*6+1):ifelse(((i-1)*6+6)>length(plist),
                                                        length(plist),
                                                        ((i-1)*6+6))],
                               ncol = nCol, nrow = nRow))
  devAskNewPage(options("device.ask.default")[[1]])
  
}
