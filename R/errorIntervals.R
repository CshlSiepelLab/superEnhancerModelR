methods::setGeneric("errorIntervals", function(x,activity,quantiles) {
  standardGeneric("errorIntervals")
})

#' Computes error bands for model
#'
#' Computes error bands for model
#' @param x enhancerDataObject
#' @param activity A range of activity values
#' @param quantiles The quantile bands to plot
#' @name errorIntervals
#' @include enhancerDataObject-class.R

methods::setMethod("errorIntervals", signature(x = "enhancerDataObject"), function(x,activity,quantiles){
  if(any(quantiles>0.5)){
    stop("Quantiles 0<quantiles<0.5")
  }
  err=list()
  express=predictExpression(x,activity)
  if(x@errorModel$type=="lognormal"){
    for(q in quantiles){
      err[[as.character(q)]]=data.frame(x=c(activity,rev(activity)),
                                        y=c(qlnorm(q, log(express), x@errorModel$value[1]),
                                            rev(qlnorm(1-q, log(express), x@errorModel$value[1]))),
                                        Quantile=paste0(q,"-",1-q))
    }
  } else if(x@errorModel$type=="gaussian"){
    for(q in quantiles){
      err[[as.character(q)]]=data.frame(x=c(activity,rev(activity)),
                                        y=c(qnorm(q, express, x@errorModel$value[1]),
                                            rev(qnorm(1-q, express, x@errorModel$value[1]))),
                                        Quantile=paste0(q,"-",1-q))
    }
  } else {
    stop("Unsupported error model")
  }
  return(do.call("rbind",err))
})
