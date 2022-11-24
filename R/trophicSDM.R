#' Fitting a trophic Species distribution model
#'
#' trophicSDM is used to fit a trophic species distribution model. Requires the species distribution data Y (the sites x species matrix), explanatory variables X and a directed acyclic graph G containing species interactions (i.e., the metaweb, with links going from predators to prey). The function fits the distribution of each species as a function of their preys (with mode = "prey", by default) or predators (if set mode = "predator").
#' @param Y The sites x species matrix containing observed species distribution (e.g. presence-absence).
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param G The species interaction network (aka metaweb). Needs to be an igraph object. Links must go from predator to preys. It needs to be a directed acyclic graph.
#' @param env.formula The definition of the abiotic part of the model. It can be :
#' \itemize{
#' \item a string specifying the formula (e.g. "~ X_1 + X_2"). In this case, the same environmental variables are used for every species.
#' 
#' \item A list that contains for each species the formula that describes the abiotic part of the model. In this case, different species can be modeled as a function of different environmental covariates. The names of the list must coincide with the names of the species.
#' }
#' @param method which SDM method to use. For now the available choices are: \code{"glm"} (frequentist) or \code{"stan_glm"} (full bayesian MCMC, default). Notice that using "glm" does not allow error propagation when predicting.
#' @param mode  "prey" if bottom-up control (default), "predators" otherwise. Notice that G needs to be such that links point from predators to prey.
#' @param family the family parameter of the glm function (see \code{glm}). \code{gaussian(link = "identity")} for gaussian data. \code{binomial(link = "logit")} or \code{binomial(link = "probit")} for presence-absence data.
#' @param iter (for \code{"stan_glm"} only) Number of iterations for each MCMC chain if stan_glm is used
#' @param chains (for \code{"stan_glm"} only) Number of MCMC chains (default to 2)
#' @param penal Penalisation method to shrink regression coefficients. If \code{NULL} (default), the model does not penalise the regression coefficient. For now, available penalization method are  \code{"horshoe"} for method stan_glm, \code{"elasticnet"} for method glm. It is also possible to constrain the sign of biotic coefficients (prey coefficients are set to positive and predator coefficients to negative) by setting \code{"coeff.signs"} for methods glm and stan_glm.
#' @param sp.formula (optional) It allows to specify a particular definition of the biotic part of the model, e.g., using composite variables (e.g., richness), or an interaction of the biotic and abitic component. More details in 'Details'.
#' @param sp.partition (optional) a list to specify groups of species that are used to compute composite variables, e.g., a species can be modeled as a function of the richness of each group of preys. It has to be a list, each element is a vector containing the names of species in the group. More details in 'Details'.
#' @param run.parallel Whether species models are fitted in parallel (can speed computational up time). Default to \code{FALSE}.
#' @param verbose Whether to print algorithm progresses
#' @return A "trophicSDMfit" object, containing:
#'    \item{model}{A list containing the local models (i.e. a SDM for each species). Each local model is an object of class "SDMfit". See \code{?SDMfit} for more informations.}
#'    \item{Y}{A numeric vector of standard errors on parameters}
#'
#'    \item{form.all}{A list describing each species formula (both biotic and abiotic terms)}
#'
#'    \item{data}{A list containing all the data used to fit the model}
#'
#'    \item{model.call}{A list containing the modeling choices of the fitted model (e.g. method, penalisation...)}
#'
#'    \item{coef}{A list containing, for each species, the inferred coefficients (with credible intervals or p-values when available)}
#'    \item{MCMC.diag}{MCMC convergence metrics, only available for MCMC methods}
#'
#'    \item{AIC}{Model's AIC}
#'
#'    \item{log.lik}{Model's log.likelihood}
#'
#' @details "sp.formula" and "sp.partition" can be combined to define any kind of composite variables for the biotic part of the formula. "sp.formula" can be :
#' \itemize{
#' \item A string defining a formula as function of "richness", e.g., \code{"richness+I(richness)^2"} (species are modeled as a function of a quadratic polynomial of their prey richness), \code{"I(richness>0)"} (species are modeled as a function of a dummy variable that is equal to 1 when at least one species is present). Importantly, when group of preys (or predators) are specified by "sp.partition", species are modeled as a function of the composite variable specified by "sp.formula" for each of their prey (or predator) groups.
#'  \item A more flexible option is to specify sp.formula as a list (whose names are species' names) that contains for each species the definition of the biotic part of the model. Notice that, in this case, the function does not check that the model is a DAG. This allow to define any kind of composite variable, or to model interactions between environmental covariates and preys (or predators).
#'}
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # Increase the number of iterations to obtain reliable results.
#' m = trophicSDM(Y,X,G, env.formula, iter = 100,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' print(m)
#' 
#' # Access local models (e.g. species "Y5")
#' m$model$Y5
#' coef(m$model$Y5)
#' # The fitted model can be plotted with `plot(m)`
#' 
#' # Fit a sparse model in the Bayesian framework with the horshoe prior
#' \donttest{
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = "horshoe", 
#'                mode = "prey", method = "stan_glm")
#' }
#' # Fit frequentist glm
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#'                
#' # With elasticnet penalty   
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = "elasticnet", 
#'                mode = "prey", method = "glm")
#'
#' #### Composite variables
#' # See vignette 'Composite variables' for a complete introduction to the use of composite variables
#' # Model species as a function of a quadratic polynomial of prey richness
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                sp.formula = "richness + I(richness^2)",
#'                mode = "prey", method = "glm")
#' m$form.all
#' # Notice that for predators that feed on a single prey (with presence-absence data),
#' # their richness and the square of their richness is exactly the same variable
#' # In this case, `trophicSDM()` removes the redundant variable but prints a warning message
#'
#' # Model species as a function of a dummy variable saying whether they have at leaste one prey
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                sp.formula = "I(richness>0)",
#'                mode = "prey", method = "glm")
#' m$form.all
#'
#' # Define group of preys and model species as a function of the richness (with a quadratic term)
#' # of these groups of preys separately
#' 
#' # Species Y1 and Y2 belong to the same group, species Y3 and Y4 are both alone in their group and 
#' # species Y5 and Y6 form another group
#' sp.partition = list(c("Y1","Y2"),c("Y3"),c("Y4"), c("Y5","Y6"))
#' 
#' m = trophicSDM(Y,X,G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                sp.partition = sp.partition,
#'                sp.formula = "richness + I(richness^2)",
#'                mode = "prey", method = "glm")
#' m$form.all
#' 
#' @importFrom parallel mclapply detectCores
#' @importFrom bayesplot rhat neff_ratio
#' @importFrom  igraph is_igraph V decompose vcount
#' @export

trophicSDM = function(Y, X, G,
                      env.formula = NULL, sp.formula = NULL, sp.partition = NULL,
                      penal = NULL, mode = "prey", method = "stan_glm",
                      family, iter=500, chains = 2, run.parallel = FALSE, verbose = FALSE){

  ################################
  # checks & errors
  if(is.null(colnames(Y))) stop("Please provide species names as the column names of Y")
  
  if(!is_igraph(G)) stop("G is not an igraph object")
  
  if(is.null(env.formula)){
    env.formula = as.list(rep("~ 1", ncol(Y)))
    names(env.formula) = colnames(Y)
  }else{
    if(length(env.formula) == 1){
      env.formula = as.list(rep(env.formula, ncol(Y)))
      names(env.formula) = colnames(Y)
    }
  }
  
  if(!method %in% c("glm","stan_glm")) stop("the selected method has not yet been implemented or it has been misspelled")
  
  if(!(is.null(penal) || penal %in% c("horshoe","elasticnet","coeff.signs"))) stop("the selected penalisation has not yet been implemented or it has been misspelled")

  if(!is.null(sp.partition) &
     (!all(colnames(Y) %in% unlist(sp.partition)) |
      !all(unlist(sp.partition) %in% colnames(Y)))){
    stop("all species must be included in sp.partition")
  }
  if(!is.null(sp.formula) & length(sp.formula)>1){
    #did the user specified a custom sp.formula?
    custom.formula = T
    message("We don't check that G and the graph induced by the sp.formula specified by the user match, nor that the latter is a graph. Please be careful about their consistency.")
    if(!identical(names(sp.formula),colnames(Y))) stop("sp.formula has to be either NULL, richness, or a list whose name equals species names (i.e. colnames(Y))")
  }else{custom.formula=F}
  if(!(mode %in% c("prey","predator"))){stop("mode must be either 'prey' or 'predator'")}
  
  if(is.null(colnames(X))) warning("columns of X must have names in order to match the env.formula argument")

  if(is.character(family) | is.function(family)){stop("If you want to model Gaussian data, please provide family = gaussian(), eventually specifying a link.")}
  ################################


  ################################
  # Laplacian sorting of the graph (component by component)

  if(mode == "prey"){
    sortedV = igraph::V(G)[order(unlist(lapply(igraph::decompose(G), compute_TL_laplacian)), decreasing=T)]
  }else{
    sortedV = igraph::V(G)[order(unlist(lapply(igraph::decompose(G), compute_TL_laplacian)), decreasing=F)]
  }

  # initialize empty lists of models
  m = form.all = as.list(vector(length=igraph::vcount(G)))
  names(m) = names(form.all) = names(sortedV)

  # core part: loop on the species to fit their distribution
  if(run.parallel){
    m=mclapply(1:vcount(G),FUN = function(j){
      if(verbose) print(paste("--- Species", names(sortedV[j]), "---"))

      if(custom.formula){
        sp.form = sp.formula[[names(sortedV[j])]]
      }else{
        sp.form = sp.formula
      }

      # call a function that does a SDM with j as focal species, automatically finds covariates from G and formula.foc
      temp.mod = SDMfit(focal=names(sortedV[j]), Y, X, G,
                        sp.formula = sp.form, sp.partition, mode = mode,
                        formula.foc=paste(as.character(env.formula[[names(sortedV[j])]]),collapse=" "),
                        method=method, penal=penal, family=family, iter=iter,verbose=verbose,chains=chains)

      # assign models and env.formula (that now includes biotic variables too)
      temp.mod
    },mc.cores=detectCores()-1)

    form.all=lapply(m, `[[`, 2)
    m=lapply(m, `[[`, 1)
    names(m) = names(form.all) = names(sortedV)


  }else{


    for(j in 1:vcount(G)){
      if(verbose) print(paste("--- Species", names(sortedV[j]), "---"))

      if(custom.formula){
        sp.form = sp.formula[[names(sortedV[j])]]
      }else{
        sp.form = sp.formula
      }

      #call a function that does a SDM with j as focal species, automatically finds covariates from G and formula.foc
      temp.mod = SDMfit(focal = names(sortedV[j]), Y, X, G, sp.formula = sp.form, sp.partition, mode = mode,
                        formula.foc = paste(as.character(env.formula[[names(sortedV[j])]]),collapse=" "),
                        method = method, penal = penal, family = family,
                        chains = chains,iter = iter,verbose = verbose)

      # assign models and env.formula (that now includes biotic variables too)
      m[[names(sortedV[j])]] = temp.mod$m
      form.all[[names(sortedV[j])]] = temp.mod$form.all
    }

  }
  # return values
  trophicSDMfit = list(model = m, form.all = form.all,
                       data = list(X = X, Y = Y, G = G,n = nrow(X),
                                   p = ncol(X), S = ncol(Y), sp.name = colnames(Y)),
                       model.call = list(form.env = env.formula,
                                         mode = mode, sp.formula = sp.formula,
                                         sp.partition = sp.partition,
                                         method = method, penal = penal,
                                         iter = iter, family = family))

  if(!is.null(penal)){ if(penal == "coeff.signs"){
  trophicSDMfit$coef = lapply(trophicSDMfit$model,function(x) x$coef)
  }}
  
  trophicSDMfit$AIC = do.call(sum, lapply(trophicSDMfit$model, function(x) x$AIC))

  trophicSDMfit$log.lik = do.call(sum, lapply(trophicSDMfit$model, function(x) x$log.lik))

  if(method == "stan_glm"){
    mcmc.diag = data.frame(rhat = unlist(lapply(
      trophicSDMfit$model, function(x) rhat(x$model))),
      neff.ratio = unlist(lapply(
        trophicSDMfit$model, function(x) neff_ratio(x$model))))

    mcmc.diag$species  = sub("\\..*", "",rownames(mcmc.diag))
    mcmc.diag$coef= sub(".*\\.", "",rownames(mcmc.diag))
    
    trophicSDMfit$mcmc.diag = mcmc.diag
  
  }
  class(trophicSDMfit) = "trophicSDMfit"

  return(trophicSDMfit)
}
