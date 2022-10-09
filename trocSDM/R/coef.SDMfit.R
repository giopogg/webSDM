#' Gets regression coefficients from a local model, i.e. a SDMfit object.
#' @param object A SDMfit object, typically obtained with trophicSDM() and available in the field $model of a trophicSDMfit object
#' @param standardise Whether to standardise regression coefficients. Default to FALSE. If TRUE, coefficients are standardised using the latent variable standardisation (see Grace et al. 2018) for more details.
#' @param level The confidence level of credible intervals, only available for stan_glm method. Default to 0.95.
#' @param ... 	additional arguments
#' @return A table containing the inferred coefficients (with credible intervals or p-values when available).
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @author Giovanni Poggiato
#' @examples
#'
#' @method coef SDMfit
#' @export


coef.SDMfit = function(object, standardise = F, level = 0.95, ...){

  SDM = object
  
  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  if(SDM$method == "glm"){

    if(!standardise){

      if(is.null(SDM$penal)){

        table = cbind(estimate = SDM$model$coefficients, p_val = summary(SDM$model)$coefficients[,4])

        return(table)

      }else{

        if(SDM$penal == "elasticnet"){

          temp_coef = as.matrix(coef(SDM$model))
          names(temp_coef) = rownames(coef(SDM$model))
          colnames(temp_coef) = "estimate"

          return(temp_coef)

        }
      }
    }else{

      if(is.null(SDM$penal)){

        # standardise using the same formula in piecewiseSEM::coefs
        # It does a classical standardisation for the x, then it divides by the sd of the latent y
        # It comes from var(y*) = var(\beta*X) + pi^2/3. See Grace et al. 2018 Ecosphere for more details

        temp_coef =
          SDM$model$coefficients[-1]*apply(dplyr::select(SDM$data, -y), 2, sd)/
          sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

        temp_coef = c("(Intercept)"= coef(SDM)[1], temp_coef)

        table = cbind(estimate = temp_coef, p_val = summary(SDM$model)$coefficients[,4])

        return(table)

      }else{

        if(SDM$penal == "elasticnet"){

          temp_coef = coef(SDM)
          # standardise using the same formula in piecewiseSEM::coefs
          # It does a classical standardisation for the x, then it divides by the sd of the latent y
          # It comes from var(y*) = var(\beta*X) + pi^2/3. See Grace et al. 2018 Ecosphere for more details
          temp_coef =
            temp_coef[-1]*apply(dplyr::select(SDM$data, -y), 2, sd)/
            as.numeric(sqrt(var(predict(SDM$model, link = T,
                                        newx = as.matrix(dplyr::select(SDM$data,-y)))) + pi^2/3))

          temp_coef = c(coef(SDM)[1], temp_coef)
          temp_coef = as.matrix(data.frame(estimate = temp_coef))
          return(temp_coef)

        }
      }

    }
  }

  if(SDM$method == "stan_glm"){

    if(!standardise){

      table = cbind(mean = SDM$model$coefficients, posterior_interval(SDM$model, prob = level))
      return(table)

    }else{

      # standardise using the same formula in piecewiseSEM::coefs
      # It does a classical standardisation for the x, then it divides by the sd of the latent y
      # It comes from var(y*) = var(\beta*X) + pi^2/3. See Grace et al. 2018 Ecosphere for more details

      # standardise the mean
      estimate =
        SDM$model$coefficients[-1]*apply(dplyr::select(SDM$data, -y), 2, sd)/
        sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

      # standardise the quantiles
      table_quantile = posterior_interval(SDM$model, prob = level)[-1,]*apply(dplyr::select(SDM$data, -y), 2, sd)/
        sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

      table = cbind(estimate, table_quantile)

      table = rbind("(Intercept)" = SDM$model$coefficients[1], table)
      return(table)

    }
  }
}
