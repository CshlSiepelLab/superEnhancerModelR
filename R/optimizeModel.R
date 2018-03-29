methods::setGeneric("optimizeModel", function(params,object) {
  standardGeneric("optimizeModel")
})

#' Predicts expression
#'
#' Predicts expression from activity values and link function
#' @param x enhancerDataObject
#' @name optimizeModel
#' @include enhancerDataObject-class.R
methods::setMethod("optimizeModel", signature(params="numeric",object = "enhancerDataObject"), function(params,object){
  object@linkFunction$value=params[1:length(object@linkFunction$value)]
  object@errorModel$value=params[(length(object@linkFunction$value)+1):length(params)]
  names(object@linkFunction$value)=names(object@linkFunction$constraints$lower)
  names(object@errorModel$value)=names(object@errorModel$constraints$lower)
  return(-ll(object))
})
