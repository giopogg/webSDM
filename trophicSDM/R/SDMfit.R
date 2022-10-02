# SDMfit : the core function that fit a (classic) SDM on the focal species, with its neighbours (preys if bottom-up approach, predators if top-down approach or both) and the environmental covariates as predictors.
## Parameters:
# focal : the vertex of the graph (i.e. the species) to be modelled
# other parameters: see trophicSDM


# nice stuff about formulas https://www.datacamp.com/community/tutorials/r-formula-tutorial

SDMfit=function(focal, Y, X, G, formula.foc, sp.formula, sp.partition,
                mode = "prey", method="stan_glm",family = NULL, penal = NULL,
                iter = 1000, chains = 2,
                verbose = T){

  # for the particular case of the stan_glm model with "coeff.signs" constraint
  # another stan-for-R package is used (brsm) and formulas have to be adapted
  useBRMS = !is.null(penal) && penal=="coeff.signs" & method=="stan_glm"


  ## Build environmental part of the formula
  if(useBRMS){
    # "tempX" is an intermediary variable allowing to set lower/upper-bounds on priors
    form.env = "y ~ tempX"
    form.env.brms = paste("tempX",as.character(formula.foc))
  }else form.env = paste("y",as.character(formula.foc))


  ## Build final formula
  neigh = names(neighbors(G,focal,mode=ifelse(mode=="prey","out","in")))

  ### Include neigh as covariates
  if (length(neigh)>0){
    formulas = buildFormula(form.init=form.env,
                            species = neigh, sp.formula=sp.formula,
                            sp.partition=sp.partition,
                            useBRMS)

    form.all = formulas$form.all

    if (useBRMS){
      # intermediary species variables are themselves defined from observed species variables
      form.neigh.brms = formulas$form.brms
    }
  }else{
    # the focal species is basal
    form.all = as.character(form.env)
  }


  ## Build the data matrix to be given to the fit functions
  data = data.frame(X,Y)
  data = dplyr::rename(data,y=all_of(focal))  # rename the focal column to be called "y"


  ## Fit models depending on the chosen method and penalisation

  ### GLM
  if(method=="glm"){
    if(is.null(penal)){
      m = glm(form.all, data=data, family=family)
    }else{
      if(penal=="horshoe") stop("Horshoe penalisation is not possible with glm")
      if(penal=="elasticnet"){
        y = as.matrix(model.frame(form.all,data=data)$y)
        x = as.matrix(model.frame(form.all,data=data)[,-1])

        m = glmnet(y=y,x=x,family=family$family,
                   lambda=cv.glmnet(y=y,x=x,family=family$family)$lambda.1se)
      }
      if(penal=="coeff.signs"){

        stop("coeff.signs not available for glm")

      }
    }
  }


  ### stan_glm
  if(method == "stan_glm"){

    refresh = ifelse(verbose,1000,0)

    if(is.null(penal)){
      m = stan_glm(form.all, data=data, family=family, iter=iter,
                   prior=normal(0,10),chains=chains,refresh=refresh)
    }else{
      if(penal=="elasticnet")  stop("Elastic net is not stan_glm")
      if(penal=="horshoe"){

        HS_prior = hs(df = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)
        m = stan_glm(form.all, data = data, family = family, prior = HS_prior,
                     iter = iter,chains = chains, refresh = refresh)

      }

      if(penal=="coeff.signs"){
        # constrain the signs of trophic interaction coefficients to remain positive (preys->preds) or negative (preds->preys)
        variables = strsplit(as.character(form.all), "y ~ |\\+")[[1]][-1]
        # choose priors depending on the nature of the variables
        priors = Reduce(rbind, lapply(variables, function(VAR){
          VAR <- gsub("\\^", "E", gsub('\\(|\\)', "", VAR))  # correct name for polynomial terms
          if (grepl("X", VAR)) return(set_prior("normal(0,10)", class="b", nlpar=VAR))
          if(mode=="prey"){
            if (sub("temp","",VAR) %in% neigh) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))
            if (grepl("richness",VAR)) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))
            if (grepl("customform",VAR)) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))

          }else{
            if (sub("temp","",VAR) %in% neigh) return(set_prior("normal(0,5)", class="b", ub=0, nlpar=VAR))
            if (grepl("richness",VAR)) return(set_prior("normal(0,5)", class="b", ub=0, nlpar=VAR))
            if (grepl("customform",VAR)) return(set_prior("normal(0,5)", class="b", lb=0, nlpar=VAR))

          }
          warning(paste("Variable", sub("temp","",VAR), "not present in any category"), call.=F)
        }))

        # each temporary variable depends on our initial variables ("false" non-linear model)
        form.all = bf(form.all, nl=TRUE) + lf(form.env.brms)
        if(length(neigh)>0) form.all = form.all + lf(form.neigh.brms)[[1]]

        if(family$family == "binomial"){
        m = brm(form.all, data=data, family=bernoulli(link = family$link),
                iter=iter, prior=priors, control=list(adapt_delta = 0.9))
        }else{
          m = brm(form.all, data=data, family=family, iter=iter, prior=priors, control=list(adapt_delta = 0.9))

        }
      }
    }
  }

  # create matrix of coefficient with postmean, post quantiles and whole samples

  # return
  SDM = list(model = m, form.all = form.all, method = method,family = family, penal = penal,
             iter = iter, chains = chains, sp.name = focal)
  class(SDM) = "SDMfit"

  if("glm" %in% class(SDM$model))  SDM$data = SDM$model$model
  if("glmnet" %in% class(SDM$model)) SDM$data = data.frame(y = y,x)


  SDM$coef = coef(SDM) # add coefficients by calling coef.SDMfit

  return(list(m = SDM, form.all = form.all))
}
