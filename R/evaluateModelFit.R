#' Evaluates prediction goodness of fit
#' 
#' Evaluate goodness of fit by comparing a true versus a predicted dataset of species distribution. Ypredicted is typically predicted using a prediction method of trophicSDM (in cross-validation if \code{trophicSDM_CV()} is used).
#' @param tSDM A trophicSDMfit object obtained with \code{trophicSDM()}.
#' @param Ynew A sites x species matrix containing the true species occurrences state. If set to NULL (default), it is set to the species distribution data Y on which the model is fitted.
#' @param Ypredicted A sites x species matrix containing the predicted species occurrences state. If set to NULL (default), it is set to the fitted values, i.e. predictions on the dataset used to train the model.
#' @return A table specifying the goodness of fit metrics for each species. For presence-absence data, the model computes TSS and AUC. For Gaussian data, the R2.
#' @author Giovanni Poggiato
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., and Byrnes, J. E. K.. 2018. Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
#' @importFrom dismo evaluate
#' @importFrom stats cor
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#' # Evaluate the quality of model predictions on the training
#' # Predict (fullPost = FALSE) as we used stan_glm to fit the model
#' # but here we are only intested in the posterior mean
#' Ypred = predict(m, fullPost = FALSE)
#' # format predictions to obtain a sites x species dataset whose
#' # columns are ordered as Ynew
#' Ypred = do.call(cbind,
#'                 lapply(Ypred, function(x) x$predictions.mean))
#'                 
#' Ypred = Ypred[,colnames(Y)]
#' evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
#' 
#' # Note that this is equivalent to `evaluateModelFit(m)`
#' # If we fitted the model using "glm"
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' Ypred = predict(m, fullPost = FALSE)
#' # format predictions to obtain a sites x species dataset whose
#' # columns are ordered as Ynew
#' Ypred = do.call(cbind, Ypred)
#' Ypred = Ypred[,colnames(Y)]
#' 
#' evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
#' # Note that this is equivalent to:
#' \donttest{
#' evaluateModelFit(m)
#' }
#' @export

evaluateModelFit = function(tSDM, Ynew = NULL, Ypredicted = NULL){

  if(!inherits(tSDM, "trophicSDMfit")) stop("tSDM is not an object of class trophicSDMfit" )

  if(is.null(Ynew)){

    message("You did not provide Ynew, the observed species distribution Y is used as default.")
    Ynew = tSDM$data$Y
  }

  if(is.null(Ypredicted)) {

    message("You did not provide Ypredicted, species predictions are obtained using predict()")

    if(tSDM$model.call$method == "glm") pred_samples =1
    if(tSDM$model.call$method == "stan_glm") pred_samples = tSDM$model.call$iter/10


    Ypredicted = predict(object = tSDM, pred_samples = pred_samples, fullPost = F)

    if(tSDM$model.call$family$family == "binomial"){
      if(tSDM$model.call$method == "glm") {
        Ypredicted = do.call(cbind, Ypredicted)
      }
      if(tSDM$model.call$method == "stan_glm"){
        Ypredicted = do.call(cbind, lapply(Ypredicted, function(x) x$predictions.mean))
      }

    }
    if(tSDM$model.call$family$family == "gaussian"){
      if(tSDM$model.call$method == "glm") {
        Ypredicted = do.call(cbind, lapply(Ypredicted))
      }
      if(tSDM$model.call$method == "stan_glm"){
        Ypredicted = do.call(cbind, lapply(Ypredicted, function(x) x$predictions.mean))
      }
    }

    Ypredicted = Ypredicted[,colnames(Ynew)]
  }

  if(!all(colnames(Ypredicted) %in% tSDM$data$sp.name) | !all(colnames(Ynew) %in% tSDM$data$sp.name)){
    stop("colnames of Ynew and Ypredicted must be the same of tSDM$data$sp.name)")
  }

  if(!all(colnames(Ypredicted) == colnames(Ynew))){
    stop("colnames of Ynew and Ypredicted must coincide")
  }

  S = tSDM$data$S

  # Compute Joint TSS and AUC
  if(tSDM$model.call$family$family == "binomial"){
    auc = tss = vector(length = S)

    eval = lapply(1:S,function(x){
      eval = dismo::evaluate(p = Ypredicted[which(Ynew[,colnames(Ypredicted)[x]]==1),x],
                             a = Ypredicted[which(Ynew[,colnames(Ypredicted)[x]]==0),x] )

      data.frame(auc = eval@auc, tss = max(eval@TPR+eval@TNR-1))
    })

    metrics = do.call(rbind,eval)
    metrics$species = tSDM$data$sp.name
  }

  if(tSDM$model.call$family$family == "gaussian"){

    R2 = vector(length = S)

    eval = sapply(1:S, function(s) cor(Ynew[,s], Ypredicted[,s])^2)
    metrics = data.frame(R2 = eval, species = tSDM$data$sp.name)
  }
  return(metrics)

}
