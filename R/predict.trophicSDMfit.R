#' Computes predicted values from the fitted trophicSDMfit model
#' 
#' Computes predicted values from the fitted trophicSDMfit model at environmental conditions specified by \code{Xnew}. Once predictions have been obtained, their quality can eventually be evaluated with \code{evaluateModelFit()}.
#' @param object A trophicSDMfit object obtained with trophicSDM()
#' @param Xnew a matrix specifying the environmental covariates for the predictions to be made. If NULL (default), predictions are done on the training dataset (e.g. by setting Xnew = tSDM$data$X).
#' @param prob.cov Parameter to predict with trophicSDM with presence-absence data. Whether to use predicted probability of presence (prob.cov = T) or the transformed presence-absences (default, prov.cov = F) to predict species distribution.
#' @param pred_samples Number of samples to draw from species posterior predictive distribution when method = "stan_glm". If NULL, set by the default to the number of iterations/10.
#' @param run.parallel Whether to use parallelise code when possible. Can speed up computation time.
#' @param verbose Whether to print advances of the algorithm
#' @param fullPost Optional parameter for stan_glm only. Whether to give back the full posterior predictive distribution (default, fullPost = TRUE) or just the posterior mean, and 2.5% and 97.5% quantiles,
#' @param filter.table Optional, default to NULL, should be provided only if the users wants to filter some species predictions. A sites x species matrix of zeros and ones.
#' @param ... 	additional arguments
#' @return A list containing for each species the predicted value at each sites. If method = "stan_glm", then each element of the list is a sites x pred_samples matrix containing the posterior predictive distribution of the species at each sites.
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @importFrom igraph V decompose neighbors vcount
#' @importFrom parallel mclapply detectCores
#' @importFrom stats quantile
#' @method predict trophicSDMfit
#' @export
#' @examples
#' data(Y, X, G)
#' # define abiotic part of the model
#' env.formula = "~ X_1 + X_2"
#' # Run the model with bottom-up control using stan_glm as fitting method and no penalisation
#' # (set iter = 1000 to obtain reliable results)
#' m = trophicSDM(Y, X, G, env.formula, iter = 50,
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "stan_glm")
#'# We can now evaluate species probabilities of presence for the environmental conditions c(0.5, 0.5)
#' predict(m, Xnew = data.frame(X_1 = 0.5, X_2 = 0.5))
#' # Obtain 50 draws from the posterior predictive distribution of species (pred_samples = 10)
#' # using predicted presence-absences of species to predict their predators (prob.cov = TRUE)
#' # Since we don't specify Xnew, the function sets Xnew = X by default
#' Ypred = predict(m, fullPost = TRUE, pred_samples = 10, prob.cov = FALSE)
#' # We can ask the function to only give back posterior mean and 95% credible intervals with
#' # fullPost = F
#' \donttest{
#' Ypred = predict(m, fullPost = TRUE, pred_samples = 30, prob.cov = FALSE)
#' }
#' # If we fit the model using in a frequentist  way (e.g. glm)
#' m = trophicSDM(Y, X, G, env.formula, 
#'                family = binomial(link = "logit"), penal = NULL, 
#'                mode = "prey", method = "glm")
#' # We are obliged to set pred_samples = 1 
#' # (this is done by default if pred_samples is not provided)
#' # In the frequentist case, fullPost is useless.
#'  Ypred = predict(m, pred_samples = 1, prob.cov = FALSE)
predict.trophicSDMfit = function(object, Xnew = NULL, prob.cov = FALSE,
                                 pred_samples = NULL, run.parallel = FALSE, verbose = FALSE,
                                 fullPost = TRUE, filter.table = NULL, ...){

  tSDM = object
  ############################################################
  # checks & errors
  if(!inherits(tSDM, "trophicSDMfit")) stop("You must provide a trophicSDMfit object")

  if(is.null(pred_samples)){

    if(tSDM$model.call$method=="glm") pred_samples = 1
    if(tSDM$model.call$method=="stan_glm") pred_samples = round(tSDM$model.call$iter/10)

  }

  if(tSDM$model.call$method=="glm" & pred_samples != 1 ){stop("glm requires pred_sample  = 1")}

  if(is.null(Xnew)) Xnew = tSDM$data$X

  # checks for filter.table ( must be a list where each element is of length nrow(Xnew))
  if(!is.null(filter.table) &
     (!all(unlist(lapply(filter.table,function(x) length(x) == nrow(Xnew) ))) |
      length(filter.table) != tSDM$data$S)) {
    stop("filter.table must be a list where each element is of length nrow(Xnew) ")
  }
  ############################################################

  family = tSDM$model.call$family
  n = nrow(Xnew)
  S = tSDM$data$S
  G = tSDM$data$G
  mode = tSDM$model.call$mode
  method = tSDM$model.call$method
  
  # Sort species
  sp.prediction = as.list(vector(length=vcount(G)))
  
  if(mode == "prey"){
    #sortedV = igraph::V(G)[order(unlist(lapply(igraph::decompose(G), compute_TL_laplacian)), decreasing=T)]
    sortedV = topological.sort(G, mode = "in")
  }else{
    #sortedV = igraph::V(G)[order(unlist(lapply(igraph::decompose(G), compute_TL_laplacian)), decreasing=F)]
    sortedV = topological.sort(G, mode = "out")
    
  }
  
  names(sp.prediction) = sortedV$name

  ############################################################
  # presence/absence case

  if(family$family %in% c("bernoulli", "binomial")){

    Neighb = lapply(sortedV, function(sV) neighbors(G, sV, ifelse(mode == "prey", "out", "in")))

    # core loop on species (as in trophicSDM)
    for(j in 1:vcount(G)){
      # print (if verbose)
      if(verbose){
        print(paste("--- Species", names(sortedV[j]), "---"));
        print(tSDM$model[[names(sortedV[j])]]$form.all)
      }

      # neighbor species (hat have already been predicted)
      neighb.sp = Neighb[[sortedV$name[j]]]

      # create data to predict with a call to predict.SDMfit
      newdata = array(data = NA, dim = c(n, (ncol(Xnew)+length(neighb.sp)), pred_samples))
      colnames(newdata) = c(colnames(Xnew), names(neighb.sp))

      # fill the abiotic variables
      newdata[,1:ncol(Xnew),] = as.matrix(Xnew)

      #### fill the biotic part of new data

      ## species that have already been predicted
      if (length(neighb.sp)>0){ for(k in 1:length(neighb.sp)){
        if(prob.cov){
          # if prob.cov=TRUE, use the species predicted probabilities of presence
          newdata[, ncol(Xnew)+k,] = sp.prediction[[names(neighb.sp[k])]]$predictions.prob
        }else{
          newdata[, ncol(Xnew)+k,] = sp.prediction[[names(neighb.sp[k])]]$predictions.bin
        }
      }
      }
      # apply the function SDMpredict to each layer of the array (MARGIN=3)
      if(run.parallel){

        pred.array = mclapply(1:dim(newdata)[3],
                              FUN = function(x){
                                predict(object = tSDM$model[[names(sortedV[j])]],
                                               newdata = newdata[,,x],
                                               pred_samples=1,
                                               prob.cov=prob.cov)},
                              mc.cores = detectCores() - 1)

      }else{

        pred.array = apply(newdata,
                           MARGIN = 3,
                           FUN = function(x){
                             predict(object = tSDM$model[[names(sortedV[j])]],
                                            newdata = x,
                                            pred_samples = 1,
                                            prob.cov = prob.cov)})

      }

      # unlist and format
      sp.prediction[[names(sortedV[j])]] = list()

      sp.prediction[[names(sortedV[j])]]$predictions.prob =
        do.call(cbind,
                lapply(pred.array,
                       FUN=function(x) x$predictions.prob))

      if(!is.null(filter.table)){

        sp.prediction[[names(sortedV[j])]]$predictions.prob.unfiltered =
          sp.prediction[[names(sortedV[j])]]$predictions.prob

        sp.prediction[[names(sortedV[j])]]$predictions.prob =
          filter.table[[names(sortedV[j])]] * do.call(cbind,
                                                      lapply(pred.array,
                                                             FUN=function(x) x$predictions.prob))

      }


      if(!prob.cov){# Here we don't care to give back $predictions.bin.unfiltered

        if(!is.null(filter.table)){
          sp.prediction[[names(sortedV[j])]]$predictions.bin =
            filter.table[[names(sortedV[j])]] *
            do.call(cbind,
                    lapply(pred.array,
                           FUN=function(x) x$predictions.bin))
        }else{
          sp.prediction[[names(sortedV[j])]]$predictions.bin =
            do.call(cbind,
                    lapply(pred.array,
                           FUN = function(x) x$predictions.bin))
        }
      }

    } # End loop on species


    ###### Wrap up results and resume (if needed)

    if(pred_samples == 1){

      for(j in 1:length(sp.prediction)){

        sp.prediction[[j]]$predictions.prob = sp.prediction[[j]]$predictions.prob[,1]

        if(!is.null(filter.table)){

          sp.prediction[[j]]$predictions.prob.unfiltered = sp.prediction[[j]]$predictions.prob.unfiltered[,1]

        }

      }

    }else{

      if(!fullPost){

        # predictions.prob
        for(j in 1:length(sp.prediction)){
          predictions.mean = rowMeans(sp.prediction[[j]]$predictions.prob)
          predictions.q975 = apply(sp.prediction[[j]]$predictions.prob, 1, quantile, 0.975)
          predictions.q025 = apply(sp.prediction[[j]]$predictions.prob, 1, quantile, 0.025)

          sp.prediction[[j]]$predictions.prob = list(predictions.mean = predictions.mean,
                                                     predictions.q025 = predictions.q025,
                                                     predictions.q975 = predictions.q975
          )
        }

        # predictions.prob.unfiltered (if needed)
        if(!is.null(filter.table)){

          for(j in 1:length(sp.prediction)){
            predictions.mean = rowMeans(sp.prediction[[j]]$predictions.prob.unfiltered)
            predictions.q975 = apply(sp.prediction[[j]]$predictions.prob.unfiltered,1,quantile,0.975)
            predictions.q025 = apply(sp.prediction[[j]]$predictions.prob.unfiltered,1,quantile,0.025)

            sp.prediction[[j]]$predictions.prob.unfiltered = list(predictions.mean = predictions.mean,
                                                                  predictions.q025 = predictions.q025,
                                                                  predictions.q975 = predictions.q975
            )
          }
        }

      }
    }

    if(is.null(filter.table)){
      sp.prediction = lapply(sp.prediction, function(x) x$predictions.prob)
    }

    return(sp.prediction = sp.prediction)


  }

  # Just the same of above but predictions are directly given (no need to choose between prob or not)
  if(family$family=="gaussian"){

    for(j in 1:vcount(G)){ # eventually modify with apply

      neigh.sp = neighbors(G,v=sortedV[j],mode="out")
      newdata=array(data=NA,dim=c(n,(ncol(Xnew)+length(neigh.sp)),pred_samples))

      colnames(newdata)=c(colnames(Xnew),names(neigh.sp))

      # fill the abiotic variables
      newdata[,1:ncol(Xnew),]=as.matrix(Xnew)

      if(length((neigh.sp)>0)){
        for(k in 1:length(neigh.sp)){

          newdata[,ncol(Xnew)+k,] = sp.prediction[[names(neigh.sp[k])]]# here select the random samples!

        }
      }
      if(run.parallel){

        pred.array = mclapply(1:dim(newdata)[3],
                              FUN = function(x){
                                predict(object = tSDM$model[[names(sortedV[j])]],
                                           newdata = newdata[,,x],
                                           pred_samples = 1,
                                           prob.cov = prob.cov)},
                              mc.cores = detectCores())

      }else{

        pred.array = apply(newdata,
                           MARGIN=3,
                           FUN = function(x){
                             predict(object = tSDM$model[[names(sortedV[j])]],
                                        newdata = x,
                                        pred_samples = 1,
                                        prob.cov = prob.cov)})

      }


      if(!is.null(filter.table)){ # eventually we could change this to give back both filtered and unfilted
        sp.prediction[[names(sortedV[j])]] = filter.table[[names(sortedV[j])]] * pred.array
      }else{
        sp.prediction[[names(sortedV[j])]] = pred.array
      }


    } # End loop on species

    if(pred_samples == 1){
      for(j in 1:length(sp.prediction)){
        sp.prediction[[j]] = sp.prediction[[j]][,1]
      }
    }else{

      
      if(!fullPost & method == "stan_glm"){
        # if !fullPost, only take the mean and intervals
        for(j in 1:length(sp.prediction)){
          predictions.mean = rowMeans(sp.prediction[[j]])
          predictions.q975 = apply(sp.prediction[[j]],1,quantile,0.975)
          predictions.q025 = apply(sp.prediction[[j]],1,quantile,0.025)

          sp.prediction[[j]] = list(predictions.mean = predictions.mean,
                                    predictions.q025 = predictions.q025,
                                    predictions.q975 = predictions.q975
          )
        }
      }
    }

    return(sp.prediction)

  }

}
