#' Fitting a single-species SDM
#'
#' SDMfit is used to fit a single species SDM, what we call a 'local model' of trophicSDM. It returns an object of class 'SDMfit'. Requires basically the same inputs of trophicSDM, with the requirement to specify with the parameter 'focal' the species that is modeled by the SDMfit.
#' @param focal the name of the species to be modeled
#' @param Y The sites x species matrix containing observed species distribution (e.g. presence-absence).
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param G The species interaction network (aka metaweb). Needs to be an igraph object. Links must go from predator to preys. It needs to be a directed acyclic graph.
#' @param formula.foc The formula for the abiotic part of the species distribution model.
#' @param method which SDM method to use. For now the available choices are: "glm" (frequentist) or "stan_glm" (full Bayesian MCMC, default). Notice that using "glm" does not allow error propagation when predicting.
#' @param mode  "prey" if bottom-up control (default), "predators" otherwise. Notice that G needs to be such that links point from predators to prey.
#' @param family the family parameter of the glm function (see glm). family=gaussian(link ="identity") for gaussian data or family=binomial(link = "logit") or binomial(link = "probit") for presence-absence data.
#' @param iter (for method="stan_glm" only) Number of iterations for each MCMC chain if stan_glm is used
#' @param chains (for method="stan_glm" only) Number of MCMC chains (default to 2)
#' @param penal (optional, default to NULL) Penalisation method to shrink regression coefficients.If NULL (default), the model does not penalise the regression coefficient. For now, available penalisation method are "horshoe" for stan_glm, "elasticnet" for glm and  "coeff.signs" (prey coefficients are set to positive and predator coefficients to negative) for glm and stan_glm.
#' @param sp.formula (optional) It allows to specify a particular definition of the biotic part of the model, e.g., using composite variables (e.g., richness), or an interaction of the biotic and abiotic component. More details in 'Details'.
#' @param sp.partition (optional) a list to specify groups of species that are used to compute composite variables, e.g., a species can be modeled as a function of the richness of each group of preys. It has to be a list, each element is a vector containing the names of species in the group.
#' @param verbose Whether to print algorithm progresses
#' @export
#' @return A list containing 'm', a "SDMfit" object and 'form.all', a string describing the formula of the SDMfit object. The "SDM" fit object contains:
#'    \item{model}{The output of the function used to fit the SDM. E.g., an object of class "glm" is method = "glm", an object of class "stanreg" if method = "stan_glm".}
#'    \item{Y}{A numeric vector of standard errors on parameters}
#'
#'    \item{form.all}{The formula used to fit the SDM (both abiotic and biotic terms)}
#'
#'    \item{method, family, penal, iter, chains}{The input parameters used to fit the SDM.}
#'
#'    \item{sp.name}{The name of the species modeled}
#'
#'    \item{data}{The model.frame data.frame used to fit the model}
#'
#'    \item{coef}{The inferred coefficients (with credible intervals or p-values when available)}
#'
#'    \item{AIC}{The AIC of the local model}
#'
#'    \item{log.lik}{The log.likelihood of the local model}
#'
#' @details "sp.formula" and "sp.partition" can be combined to define any kind of composite variables for the biotic part of the formula. "sp.formula" can be :
#' \itemize{
#' \item A string defining a formula as function of "richness". E.g., sp.formula="richness+I(richness)^2" (species are modeled as a function of a quadratic polynomial of their prey richness), "I(richness>0)" (species are modeled as a function of a dummy variable that is equal to 1 when at least one species is present). Importantly, when group of preys (or predators) are specified by "sp.partition", species are modeled as a function of the composite variable specified by "sp.formula" for each of their prey groups.\cr
#' \item A more flexible option is to specify sp.formula as a list (whose names are species' names) that contains for each species the definition biotic part of the model. Notice that, in this case, the function does not check that the model is a DAG. This allow to define any kind of composite variable, or to model interactions between environmental covariates and preys (or predators).
#' }
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @importFrom igraph neighbors
#' @importFrom dplyr rename all_of
#' @importFrom stats glm model.matrix logLik deviance model.frame coef
#' @importFrom glmnet glmnet
#' @importFrom rstanarm stan_glm
#' @importFrom rstanarm hs normal
#' @importFrom brms set_prior brm bf lf bernoulli
#' @importFrom glmnet cv.glmnet
#' @examples
#' data(Y,X,G)
#' # Run a local model (i.e. a SDM) for species Y6
#' mySDM = SDMfit("Y6", Y, X, G, "~X_1 + X_2", mode = "prey",
#'        method = "stan_glm", family = binomial(link = "logit"))
#' mySDM$m
#' @export

SDMfit = function(focal, Y, X, G, formula.foc, sp.formula = NULL, sp.partition = NULL,
                  mode = "prey", method = "stan_glm", family, penal = NULL,
                  iter = 1000, chains = 2,
                  verbose = TRUE){

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

  data = data.frame(X,Y)
  data = dplyr::rename(data,y=all_of(focal))  # rename the focal column to be called "y"

  ### Include neigh as covariates
  if (length(neigh)>0){
    formulas = buildFormula(form.init=form.env,
                            species = neigh, sp.formula=sp.formula,
                            sp.partition=sp.partition,
                            useBRMS)

    form.all = formulas$form.all

    #### Remove duplicate variables in formula
    if(!useBRMS){
      temp.mf = model.frame(form.all, data)
      temp.mf.all = temp.mf[!duplicated(as.list(temp.mf))]

      if(length(colnames(temp.mf.all)) != length(colnames(temp.mf))){
        message("Formula was modified since it led to identical columns of the design matrix (e.g. Y_1 or Y_1^2 for binary data)")}

      form.all = paste0("y ~ ", paste(colnames(temp.mf.all)[-1], collapse = "+"))
    }else{

      for(i in 1:length(formulas$form.brms)){
        form.temp = formulas$form.brms[[i]]
        form.temp = sub(".*~", "y ~", as.character(form.temp))

        temp.mf = model.frame(form.temp, data)
        temp.mf.all = temp.mf[!duplicated(as.list(temp.mf))]

        if(length(colnames(temp.mf.all)) != length(colnames(temp.mf))){
          message("Formula was modified since it led to identical columns of the design matrix (e.g. Y_1 or Y_1^2 for binary data)")}

        formulas$form.brms[[i]] = as.formula(paste0(sub("()","",formulas$form.brms[[i]][2]),
                                                    "~ 0 +", paste(colnames(temp.mf.all)[-1], collapse = "+")))
      }
    }

    if (useBRMS){
      # intermediary species variables are themselves defined from observed species variables
      form.neigh.brms = formulas$form.brms
    }
  }else{
    # the focal species is basal
    form.all = as.character(form.env)
  }


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

  if(!useBRMS){
    SDM$data = as.data.frame(cbind(y = SDM$data$y, model.matrix(as.formula(SDM$form.all),SDM$data)))

    SDM$coef = coef(SDM) # add coefficients by calling coef.SDMfit
  }

  # Add AIC and log.lik
  if("glm" %in% class(SDM$model)){
    SDM$AIC = SDM$model$aic
    SDM$log.lik = logLik(SDM$model)
  }
  if("glmnet" %in% class(SDM$model)){
    tLL = -deviance(SDM$model)
    SDM$log.lik = tLL/2
    k = ncol(SDM$data)
    SDM$AIC <- -tLL+2*k
  }
  if(SDM$method == "stan_glm"){
    SDM$log.lik = sum(colMeans(log_lik(SDM$model)))
    k = ncol(SDM$data) # includes intercepts
    SDM$AIC <- -2*SDM$log.lik + 2*k
  }

  return(list(m = SDM, form.all = form.all))
}
