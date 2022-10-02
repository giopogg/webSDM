
## Parameters:
# m : the model
# focal : the focal species to be predicted
# newdata: the corresponding predictors (formatted in trophicSDM_predict)
# pred_samples : the number of samples to be predicted. Notice that it corresponds to the third dimension of the newdata array.
# for the other parameters see trophicSDM_predict


predict.SDMfit=function(SDMfit,newdata,pred_samples,prob.cov){

  if(class(SDMfit) != "SDMfit") stop("you need to provide a SDMfit object")
  family=SDMfit$family$family
  method=SDMfit$method

  # If presence/absence data

  if(family %in% c("bernoulli", "binomial")){

    if(method=="stan_glm"){

      # retrieve the expected predicted probability of presence
      ## ! for the particular case of the stan_glm model with a constraint on the coefficients
      ## ! signs another stan-for-R package is used (brsm) and the code has to be adapted
      if (class(SDMfit$model)[1]=='brmsfit'){  # brms package
        predictions.prob=t(posterior_epred(SDMfit$model, newdata=as.data.frame(newdata), nsamples=pred_samples))
      }else{  # rstan package
        predictions.prob=t(posterior_epred(SDMfit$model, newdata=as.data.frame(newdata), draws=pred_samples))
      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov) {
        if (class(SDMfit$model)[1]=='brmsfit'){  # brms package
          predictions=t(posterior_predict(SDMfit$model, newdata=as.data.frame(newdata), nsamples=pred_samples))
        }else{  # rstan package
          predictions=t(posterior_predict(SDMfit$model, newdata=as.data.frame(newdata), draws=pred_samples))
        }
      }else{
        predictions = NULL
        }
    }


    if(method=="glm")  {
      if(pred_samples!=1) stop("pred_sample must be 1 if method=glm!")

      # retrieve the expected predicted probability of presence
      if(is.null(SDMfit$penal)){
      predictions.prob=as.matrix(predict(SDMfit$model,type = "response",newdata=as.data.frame(newdata)))
      }else{

        newx = model.frame(gsub(".*y","",SDMfit$form.all), newdata)

        predictions.prob=as.matrix(predict(SDMfit$model,type = "response",
                                           newx= as.matrix(newx)))

      }
      # retrieve a unique presence/absence prediction
      if(!prob.cov){
        predictions=as.matrix(sapply(as.vector(predictions.prob),FUN=function(x) rbinom(prob=x,size=1,n=1)))
      }else{
        predictions = NULL
        }
    }

    return(list(predictions.prob = predictions.prob, predictions.bin = predictions))

  }

  # If gaussian data

  if(family=="gaussian"){

    if(method=="stan_glm"){
      predictions=t(posterior_predict(SDMfit$model, newdata=as.data.frame(newdata),draws=pred_samples))
    }

    if(method=="glm")  {

      if(pred_samples!=1) {

        stop("pred_sample must be 1 if method=glm!")}

      predictions=as.matrix(predict(SDMfit$model,type = "response",newdata=as.data.frame(newdata)))

    }

    return(predictions)
  }
}
