#' Plots the regression coefficients of a local model
#' 
#' Plots the regression coefficients of a local SDMfit model
#' @param x A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param level the confidence level of the confidence intervals
#' @param ... 	additional arguments
#' @return A plot of the regression coefficients of the fitted local SDM
#' @author Giovanni Poggiato
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # Plot species Y6
#' \donttest{
#' plot(m$model$Y6)
#' }
#' @import ggplot2
#' @import brms
#' @import broom.mixed
#' @importFrom jtools plot_summs
#' @importFrom gridExtra grid.arrange
#' @importFrom stats coef
#' @importFrom utils globalVariables
#' @method plot SDMfit
#' @export
plot.SDMfit = function(x, level = 0.95,...){

  SDM = x
  
  if(!inherits(SDM, "SDMfit")) stop("SDM is not an object of class SDMfit" )

  if(SDM$method == "glm"){

    if(is.null(SDM$penal)){
      p = plot_summs(SDM$model, robust = TRUE,plot.distributions = TRUE, inner_ci_level = level, omit.coefs = NULL) +
        ggtitle(paste0("Species : ", SDM$sp.name))

      return(p)

    }else{

      table = data.frame(x = as.vector(coef(SDM$model)), y = rownames(coef(SDM$model)))

      p = ggplot(data = table) + geom_point(aes(x,y)) + labs(x = "", y = "") +
        geom_vline(xintercept=0, lty =2, alpha = 0.5) +
        theme_classic() +  theme(axis.text.y = element_text(face="bold", size=13)) +
        ggtitle(paste0("Species : ", SDM$sp.name))

      return(p)

    }
  }

  if(SDM$method == "stan_glm"){
    p = plot(SDM$model, prob = level, plotfun = "areas") + ggtitle(paste0("Species : ", SDM$sp.name))

    return(p)
  }
}
