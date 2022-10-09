#' Fitting a trophic Species distribution model
#'
#' trophicSDM is used to fit a trophic species distribution model. Requires the species distribution data Y (the sites x species matrix), explanatory variables X and a directed acyclic graph G containing species interactions (i.e., the metaweb, with links going from predators to prey). The function fits the distribution of each species as a function of their preys (with mode = "prey", by default) or predators (if set mode = "predator").
#' @param Y The sites x species matrix containing observed species distribution (e.g. presence-absence).
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param G The species interaction network (aka metaweb). Needs to be an igraph object. Links must go from predator to preys. It needs to be a directed acyclic graph.
#' @param env.formula A list that contains for each species the formula that describes the abiotic part of the model. The names of the list must coincide with the names of the species. The details of model specification are given under ‘Details’.
#' @param method which SDM method to use. For now the available choises are: "glm" (frequentist) or "stan_glm" (full bayesian MCMC, default). Notice that using "glm" does not allow error propagation when predicting.
#' @param mode  "prey" if bottom-up control (default), "predators" otherwise. Notice that G needs to be such that links point from predators to prey.
#' @param family the family parameter of the glm function (see glm). family=gaussian for gaussian data or family=binomial(link = "logit") or binomial(link = "probit") for presence-absence data.
#' @param iter (for method="stan_glm" only) Number of iterations for each MCMC chain if stan_glm is used
#' @param chains (for method="stan_glm" only) Number of MCMC chains (default to 2)
#' @param penal (optional, default to NULL) Penalisation method to shrink regression coefficients.If NULL (default), the model does not penalise the regression coefficient. For now, available penalisation method are "horshoe" for stan_glm, "elasticnet" for glm and  "coeff.signs" (prey coefficients are set to posite and predator coefficients to negative) for glm and stan_glm.
#' @param sp.formula (optional) It allows to specify a particular definition of the biotic part of the model, e.g., using composite variables (e.g., richness), or an interaction of the biotic and abitic component. More details in 'Details'.
#' @param sp.partition (optional) a list to specify groups of species that are used to compute composite variables, e.g., a species can be modelled as a function of the richness of each group of preys. It has to be a list, each element is a vector containing the names of species in the group.
#' @param run.parallel Whether species models are fitted in parallel (can speed computational up time). Default to FALSE.
#' @param verbose Whether to print algorithm progresses
#' @return A "trophicSDMfit" object, containing:
#'    \item{model}{A list containing the local models (i.e. a SDM for each species). Each local model is an object of class "SDMfit". See ?SDMfit for more informations.}
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
#' \item A string defining a formula as function of "richness". E.g., sp.formula="richness+I(richness)^2" (species are modelled as a function of a quadratic polyome of their prey richness), "I(richness>0)" (species are modelled as a function of a dummy variable that is equal to 1 when at least one species is present). Importantly, when group of preys (or predators) are specified by "sp.partition", species are modeled as a function of the composite variable specified by "sp.formula" for each of their prey groups.
#'  \item A more flexible option is to specify sp.formula as a list (whose names are species' names) that contains for each species the definition biotic part of the model. Notice that, in this case, the function does not check that the model is a DAG. This allow to define any kind of composite variable, or to model interactions between environmental covariates and preys (or predators).
#'}
#' @author Giovanni Poggiato and Jérémy Andréletti
#' @examples
#'
#' @export

trophicSDM = function(Y, X, G,
                      env.formula = NULL, sp.formula = NULL, sp.partition = NULL,
                      penal = NULL, mode = "prey", method = "stan_glm",
                      family, iter=1000, run.parallel=TRUE, chains = 2, verbose = F){

  ################################
  # checks & errors
  if(!is_igraph(G)) stop("G is not an igraph object")
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

  ################################


  ################################
  # Laplacian sorting of the graph (component by component)

  if(mode == "prey"){
    sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=T)]
  }else{
    sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=F)]
  }

  # initialize empty lists of models
  m = form.all = as.list(vector(length=vcount(G)))
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
    },mc.cores=detectCores())

    form.all=lapply(m, `[[`, 2)
    m=lapply(m, `[[`, 1)
    names(m) = names(form.all) = names(sortedV)


  }else{


    for(j in 1:vcount(G)){
      print(paste("--- Species", names(sortedV[j]), "---"))

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

  trophicSDMfit$coef = lapply(trophicSDMfit$model,function(x) x$coef)

  trophicSDMfit$AIC = do.call(sum, lapply(trophicSDMfit$model, function(x) x$AIC))

  trophicSDMfit$log.lik = do.call(sum, lapply(trophicSDMfit$model, function(x) x$log.lik))

  if(method == "stan_glm"){
    trophicSDMfit$mcmc.diag = data.frame(rhat = unlist(lapply(
      trophicSDMfit$model, function(x) rhat(x$model))),
      neff.ratio = unlist(lapply(
        trophicSDMfit$model, function(x) neff_ratio(x$model))))

  }
  class(trophicSDMfit) = "trophicSDMfit"

  return(trophicSDMfit)
}
