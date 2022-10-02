
# level only for stan_glm

coef.SDMfit = function(SDM, standardise = F, level = 0.95){

  if(class(SDM) != "SDMfit") stop("SDM is not an object of class SDMfit" )

  if(SDM$method == "glm"){

    if(!standardise){

      if(is.null(SDM$penal)){

        table = cbind(estimate = SDM$model$coefficients, p_val = summary(SDM$model)$coefficients[,4])

        return(table)

      }else{

        if(SDM$penal == "elasticnet"){

          temp_coef = as.vector(coef(SDM$model))
          names(temp_coef) = rownames(coef(SDM$model))
          return(temp_coef)

        }
      }
    }else{

      if(is.null(SDM$penal)){

        # standardise using the same formula in piecewiseSEM::coefs
        # It does a classical standardisation for the x, then it divides by the sd of the latent y
        # It comes from var(y*) = var(\beta*X) + pi^2/3. See Grace et al. 2018 Ecosphere for more details

        temp_mean =
          SDM$model$coefficients[-1]*apply(dplyr::select(SDM$data, -y), 2, sd)/
          sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

        temp_coef = c("(Intercept)"= coef(SDM)[1], temp_mean)

        table = cbind(temp_coef, p_val = summary(SDM$model)$coefficients[,4])

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
      temp_mean =
        SDM$model$coefficients[-1]*apply(dplyr::select(SDM$data, -y), 2, sd)/
        sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

      # standardise the quantiles
      table_quantile = posterior_interval(SDM$model, prob = level)[-1,]*apply(dplyr::select(SDM$data, -y), 2, sd)/
        sqrt(var(predict(SDM$model, link = T)) + pi^2/3)

      table = cbind(temp_mean, table_quantile)

      table = rbind("(Intercept)" = SDM$model$coefficients[1], table)
      return(table)

    }
  }
}
