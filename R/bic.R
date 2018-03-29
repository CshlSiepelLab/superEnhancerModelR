methods::setGeneric("bic", function(x) {
  standardGeneric("bic")
})

#' Computes error bands for model
#'
#' Computes error bands for model
#' @param x enhancerDataObject
#' @name bic
#' @include enhancerDataObject-class.R
#' @export
methods::setMethod("bic", signature(x = "enhancerDataObject"), function(x){
  return(-2*ll(x)+(length(x@linkFunction$value)+length(x@errorModel$value))*log(length(x@expressionData)))
})
