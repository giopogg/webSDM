
###################################################################################################
### trophicSDM_predict
#trophicSDM_predict : the function that predicts the output of the trophicSDM model, eventually
#depending on a new set of environmental covariates. Several parameters are available depending
#on the users' needs.

## Parameters:
# m : the fitted model
# Xnew : Environmental data must be a sites x predictor matrix. If NULL, the fitted environmental predictors are used (i.e. in-sample prediction)
# prob.cov : TRUE if the predicted probabilities of presence of the preys are used to predict species. FALSE if the randomly (depending of predicted probabilities of course) generated presence-absence are used
# mode : type of graph neighbors used for prediction, "out" for preys only, "in" for predators only, "all" for both
# pred.sample : the total number of sample from the predictive distribution of each species at each site
# filter.table: ONLY FOR EUROPE. is of a list (of length ncol(Y)) whose elements are vector of list nrow(Xnew)
# fullPost: do you want the model to give back all the samples from the posterior? it can be TRUE, or, if not the confidence level of the credible intervals that are given back (needs to be between 0 and 1, default choice 0.95 for 95% confidence intervals). If !fullPost, only the mean and the confidence intervals of the probabilities are given.
predict.trophicSDMfit = function(tSDM, Xnew = NULL, prob.cov = F,
                                 pred_samples = NULL, run.parallel = T, verbose = F, filter.table = NULL,
                                 fullPost = T){

  ############################################################
  # checks & errors
  if(class(tSDM) != "trophicSDMfit") stop("You must provide a trophicSDMfit object")

  if(is.null(pred_samples)){

    if(tSDM$model.call$method=="glm") pred_samples = 1
    if(tSDM$model.call$method=="stan_glm") pred_samples = tSDM$model.call$iter

  }

  if(tSDM$model.call$method=="glm" & pred_samples != 1 ){stop("glm requires pred_sample  = 1")}

  if(is.null(Xnew)) Xnew = tSDM$data$X

  # checks for filter.table ( must be a list where each element is of length nrow(Xnew))
  if(!is.null(filter.table) &
     (!all(unlist(lapply(filter.table,function(x) length(x) == nrow(Xnew) ))) |
      length(filter.table) != ncol(Y))) {
    stop("filter.table must be a list where each element is of length nrow(Xnew) ")
  }
  ############################################################

  family = tSDM$model.call$family
  n = nrow(Xnew)
  S = tSDM$data$S
  G = tSDM$data$G
  mode = tSDM$model.call$mode

  # Sort species
  sp.prediction = as.list(vector(length=vcount(G)))

  if(mode=="prey"){
    sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=T)]
  }else{
    sortedV = V(G)[order(unlist(lapply(decompose(G), compute_TL_laplacian)), decreasing=F)]
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
                                predict.SDMfit(SDMfit = tSDM$model[[names(sortedV[j])]],
                                               newdata = newdata[,,x],
                                               pred_samples=1,
                                               prob.cov=prob.cov)},
                              mc.cores = detectCores()-1)

      }else{

        pred.array = apply(newdata,
                           MARGIN = 3,
                           FUN = function(x){
                             predict.SDMfit(SDMfit = tSDM$model[[names(sortedV[j])]],
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
                                SDMpredict(SDMfit = tSDM$model[[names(sortedV[j])]],
                                           newdata = newdata[,,x],
                                           pred_samples = 1,
                                           prob.cov = prob.cov)},
                              mc.cores = detectCores()-1)

      }else{

        pred.array = apply(newdata,
                           MARGIN=3,
                           FUN = function(x){
                             SDMpredict(SDMfit = tSDM$model[[names(sortedV[j])]],
                                        newdata = x,
                                        pred_samples = 1,
                                        prob.cov = prob.cov)})

      }


      if(!is.null(filter.table)){ # eventually we could change this to give back both filtered and unfilted
        sp.prediction[[names(sortedV[j])]] = filter.table[[names(sortedV[j])]] * do.call(cbind,pred.array)
      }else{
        sp.prediction[[names(sortedV[j])]] = do.call(cbind,pred.array)
      }


    } # End loop on species

    if(pred_samples == 1){
      for(j in 1:length(sp.prediction)){
        sp.prediction[[j]] = sp.prediction[,1]
      }
    }else{

      if(!fullPost){
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
