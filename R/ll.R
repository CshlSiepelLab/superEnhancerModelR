methods::setGeneric("ll", function(x) {
  standardGeneric("ll")
})

#' Computes log-likelihood
#'
#' Computes log-likelihood of the enhancerData model
#' @param x enhancerDataObject
#' @name ll
#' @include enhancerDataObject-class.R
#' @export
methods::setMethod("ll", signature(x = "enhancerDataObject"), function(x) {
  ## Predict expression levels from the model
  expression=predictExpression(x)
  if(x@errorModel$type=="gaussian"){
    ll=apply(cbind(x@expressionData,expression),1, function(z) dnorm(z[1],z[2],x@errorModel$value[1],log=TRUE))
  } else if(x@errorModel$type=="lognormal"){
    ll=apply(cbind(x@expressionData,myLog(expression)),1, function(z) dlnorm(z[1],z[2],x@errorModel$value[1],log=TRUE))
    } else {
    stop("Invalid error model.")
  }
  return(sum(ll))
})
