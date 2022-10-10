#' Plots the regression coeffiecients of a local model
#' 
#' Plots the regression coeffiecients of a local SDMfit model
#' @param x A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param level the confidence level of the confidence intervals
#' @param ... 	additional arguments
#' @author Giovanni Poggiato
#' @examples
#'
#' @method plot SDMfit
#' @export
plot.SDMfit = function(x, level = 0.95, ...){

  SDM = x
  
  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  if(SDM$method == "glm"){

    if(is.null(SDM$penal)){
      p = plot_summs(SDM$model, robust = TRUE,plot.distributions = TRUE, inner_ci_level = level, omit.coefs = NULL) +
        ggtitle(paste0("Species : ", SDM$sp.name))

      if(plot) print(p)

      return(p)

    }else{

      table = data.frame(x = as.vector(coef(SDM$model)), y= rownames(coef(SDM$model)))

      p = ggplot(data = table) + geom_point(aes(x,y)) + labs(x = "", y = "") +
        geom_vline(xintercept=0, lty =2, alpha = 0.5) +
        theme_classic() +  theme(axis.text.y = element_text(face="bold", size=13)) +
        ggtitle(paste0("Species : ", SDM$sp.name))

      print(p)
      
      return(p)

    }
  }

  if(SDM$method == "stan_glm"){
    p = plot(SDM$model, prob = level, plotfun = "areas") + ggtitle(paste0("Species : ", SDM$sp.name))

    print(p)

    return(p)
  }
}
