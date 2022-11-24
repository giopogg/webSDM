#' Compute K-fold cross-validation predicted values from a fitted trophicSDM model
#' 
#' Once the CV predicted values are obtained, their quality can be evaluated with \code{evaluateModelFit()}.
#' @param tSDM A trophicSDMfit object obtained with trophicSDM()
#' @param K The number of folds for the K-fold cross validation
#' @param partition Optional parameter. A partition vector to specify a partition in K fold for cross validation
#' @param prob.cov Parameter to predict with trophicSDM with presence-absence data. Whether to use predicted probability of presence (prob.cov = T) or the transformed presence-absences (default, prov.cov = F) to predict species distribution.
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param iter For method = "stan_glm": number of iterations of each MCMC chains to fit the trophicSDM model. Default to the number of iterations used to fit the provided trophicSDMfit object
#' @param chains For method = "stan_glm": number of MCMC chains to fit the trophicSDM model. Default to the number of iterations used to fit the provided trophicSDMfit object
#' @param run.parallel Whether to use parallelise code when possible. Default to TRUE. Can speed up computation time
#' @param verbose Whether to print advances of the algorithm
#' @return A list containing:
#' \item{meanPred}{a sites x species matrix of predicted occurrences of species for each site (e.g. probability of presence). With stan_glm the posterior predictive mean is return}
#' \item{Pred975,Pred025}{Only for method = "stan_glm", the 97.5% and 2.5% quantiles of the predictive posterior distribution}
#' \item{partition}{the partition vector used to compute the K fold cross-validation}
#' @author Giovanni Poggiato
#' @importFrom parallel mclapply detectCores
#' @importFrom abind abind
#' @importFrom stats quantile
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' \donttest{
#' m = trophicSDM(Y, X, G, env.formula, iter = 50, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#'
#' # Run a 3-fold (K=3) cross validation. Predictions is done using presence-absences of preys
#' # (prob.cov = FALSE, see ?predict.trophicSDM) with 50 draws from the posterior distribution
#' # (pred_samples = 50)
#' CV = trophicSDM_CV(m, K = 3, prob.cov = FALSE, pred_samples = 10, run.parallel = FALSE)
#' # Use predicted values to evaluate model goodness of fit in cross validation
#' Ypred = CV$meanPred[,colnames(Y)]
#' 
#' evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
#' }
#' # Now with K = 2 and by specifying the partition of site
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' partition = c(rep(1,500),rep(2,500))
#' CV = trophicSDM_CV(m, K = 2, partition = partition, prob.cov = FALSE,
#'                    pred_samples = 10, run.parallel = FALSE)
#' Ypred = CV$meanPred[,colnames(Y)]
#' evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
#' @export
trophicSDM_CV = function(tSDM, K, partition = NULL, prob.cov = FALSE,
                         pred_samples = NULL,
                         iter = NULL, chains = NULL, run.parallel = FALSE, verbose = FALSE){

  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM needs to be a trophicSDMfit object")

  # set MCMC parameters
  if(tSDM$model.call$method == "glm"){ iter = 1; pred_samples = 1
  }else{
    if(is.null(iter)) iter = tSDM$model.call$iter
    if(is.null(pred_samples)) pred_samples = round(iter/10)
  }

  n = tSDM$data$n
  S = tSDM$data$S
  G = tSDM$data$G

  # Create partition (if needed)
  if(!is.null(partition)){
    if(length(partition) != n) stop("partition must be a vector of length n (the number of sites)")
  }else{
    partition <- sample(1:K,size=n,replace=TRUE,prob=rep(n/K,K))
  }


  preds = array(dim=c(n,pred_samples,S))

  for(i in 1:K){

    train = which(partition != i)
    test = which(partition == i)
    Y = tSDM$data$Y[train,]
    X = tSDM$data$X[train,]

    m_K = trophicSDM(Y = Y, X = X, G = tSDM$data$G, env.formula = tSDM$model.call$form.env,
                     penal = tSDM$model.call$penal, method = tSDM$model.call$method, mode = tSDM$model.call$mode,
                     family = tSDM$model.call$family, iter = iter,
                     run.parallel = run.parallel, verbose = verbose,
                     chains=2)

    pred_K = predict(m_K, Xnew = tSDM$data$X[test,],
                     prob.cov = prob.cov, pred_samples = pred_samples, run.parallel = run.parallel)

    if(tSDM$model.call$method == "stan_glm"){
    preds[test,,] = abind(pred_K,along=3)
    }else{
      preds[test,,] = do.call(cbind, pred_K)
    }

    print(paste0("Fold ", i, " out of ", K," \n"))

  }

  if(tSDM$model.call$method == "glm"){

    meanPred = preds[,1,]

    colnames( meanPred ) = names(pred_K)
    
    Pred975 = Pred025 = NULL

  }
  
  if(tSDM$model.call$method == "stan_glm"){

  meanPred = apply(preds,mean,MARGIN = c(1,3))
  Pred975 = apply(preds, quantile, MARGIN=c(1,3),0.975)
  Pred025 = apply(preds, quantile, MARGIN=c(1,3),0.025)

  colnames( meanPred ) = colnames(Pred975) = colnames(Pred025) = names(pred_K)
  }


  list(meanPred = meanPred, Pred975 = Pred975, Pred025 = Pred025, partition = partition)

}
