# Trophic SDM
## Parameters:
# Y : Community data must be a sites x species matrix
# X : Environmental data must be a sites x predictor matrix
# G : Trophic network, igraph
# env.formula : list (names of each element are the name of vertices) of one-hand formulas whose names are the same as in X.
# sp.formula : formula for composite variable. It HAS TO BE a string but written as the rhs of a formula. It can include quadratic terms as formula (also poly works). Be careful ONLY "richness". For example, you could set sp.formula="richness" or sp.formula="richness+richness^2" or sp.formula="sin(richness^2)" (no, you won't)
# sp.partition : list where the edges (i.e. their name as string) are grouped. Each element of the list is a group. E.g. sp.partition = list(c("Y_1","Y_2"),"Y_3",c("Y_4","Y_5","Y_6"))
# penal: how to penalise the likelihood to reduce the dimension. Each SDM method has its own penalisation type, some SDM method have no available penalisation approach. For now, we have "horshoe" for stan_glm, "elasticnet" for glm and  "coeff.signs" (positive prey interactions and negative predator interactions) for glm and stan_glm
# method: which SDM method to use (available for now: glm (frequentist), stan_glm (full bayesian MCMC))
# family: the family parameter of the glm function (see glm). family=gaussian for gaussian data or family=binomial(logit) or binomial(probit) for presence-absence data
# mode : "prey" if bottom-up control, "predators" otherwise
# iter : the number of MCMC samples (if bayesian approach)


#' Fitting a trophic Species distribution model
#'
#' trophicSDM is used to fit a trophic species distribution model. Requires the response variable Y (the sites x traits matrix), explanatory variables X and a directed acyclic graph G containing species interactions (i.e. the metaweb). The function fits the distribution of each species as a function of their preys (default) or predators (set mode = "predator").
#' @param Y The sites x species matrix containing observed species distribution (e.g. presence-absence).
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param G The species interaction network (aka metaweb). Needs to be an igraph object. Links must go from predator to preys. It needs to be a directed acyclic graph
#' @param env.formula A list that contains for each species the formula that describes the abiotic part of the model. The names of the list must coincide with the names of the species. The details of model specification are given under ‘Details’
#' @param method which SDM method to use. For now the choises are: "glm" (frequentist) or "stan_glm" (full bayesian MCMC, default)). Notice that using "glm" does not allow error propagation when predicting.
#' @param mode  "prey" if bottom-up control (default), "predators" otherwise. Notice that G needs to be such that links point from predators to prey.
#' @param family the family parameter of the glm function (see glm). family=gaussian for gaussian data or family=binomial(link = "logit") or binomial(link = "probit") for presence-absence data.
#' @param iter (for method="stan_glm" only) Number of iterations for each MCMC chain if stan_glm is used
#' @param chains (for method="stan_glm" only) Number of MCMC chains (default to 2)
#' @param penal (optional default to NULL) How to penalise the likelihood to reduce the dimension. Each SDM method has its own penalisation type, some SDM method have no available penalisation approach. For now, we have "horshoe" for stan_glm, "elasticnet" for glm and  "coeff.signs" (positive prey interactions and negative predator interactions) for glm and stan_glm
#' @param sp.formula (optional) can be specified to specify a particular definition of the biotic part of the model, e.g., using composite variables (e.g., richness), or an interaction of the biotic and abitic component. More details in 'Details'. It order to specify composite variables, it has to be a string written as the rhs of a formula as a function of "richness". For example, you could set sp.formula="richness" or sp.formula="richness+I(richness)^2" or "I(richness>0)". Richness can be computed for different groups of species, when the parameter sp.partition is specified. A more flexible option is to specify sp.formula as a list (whose names are species' names) that contain for each species the biotic part of the model. Notice that, in this case, the function does not check that the model is a DAG.
#' @param sp.partition (optional) a list specifying groups of species, in order to compute composite variables based separetly for each group of prey. It has to be a list where the species (i.e. their name as string) are grouped. Each element of the list is a group. E.g. sp.partition = list(c("Y_1","Y_2"),"Y_3",c("Y_4","Y_5","Y_6")).
#' @param run.parallel Whether species models are fitted in parallel (can speed computational up time). Default to FALSE.
#' @param verbose Whether to print informations on model convergence.
#' @export

#' @return A list containing:
#'    \item{model}{ An object of class 'runjags' containing the fitted model.}
#'    \item{Y}{A numeric vector of standard errors on parameters}
#'
#'    \item{X_raw}{The design matrix specified as input}
#'
#'    \item{X}{The design matrix transformed as specified in formula}
#'
#'    \item{formula}{The formula specified as input}
#' }
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See formula for more details of allowed formulae.
#' @examples
#' data(Y)
#' data(X)
#' # Short MCMC to obtain a fast example: results are unreliable !
#'m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,
#'         burnin = 100,
#'         sample = 100)
#'
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

  if(method == "stan_glm"){
    trophicSDMfit$mcmc.diag = data.frame(rhat = unlist(lapply(
      trophicSDMfit$model, function(x) rhat(x$model))),
      neff.ratio = unlist(lapply(
        trophicSDMfit$model, function(x) neff_ratio(x$model))))

  }
  class(trophicSDMfit) = "trophicSDMfit"

  return(trophicSDMfit)
}
