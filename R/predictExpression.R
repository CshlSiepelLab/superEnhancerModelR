methods::setGeneric("predictExpression", function(x,activity) {
  standardGeneric("predictExpression")
})

#' Predicts expression
#'
#' Predicts expression from object activity values and link function
#' @param x enhancerDataObject
#' @name predictExpression
#' @include enhancerDataObject-class.R
#' @export
methods::setMethod("predictExpression", signature(x = "enhancerDataObject"), function(x) {
  if(x@linkFunction$type=="additive"){
    expression=computeActivity(x)
  } else if(x@linkFunction$type=="exponential"){
    expression=exp(computeActivity(x))
  } else if(x@linkFunction$type=="logistic"){
    expression=x@linkFunction$value[length(x@linkFunction$value)]/(1+exp(-computeActivity(x)))
  } else {
    stop("Invalid link function")
  }
  return(expression)
})


#' Predicts expression
#'
#' Predicts expression from supplied activity values and link function
#' @param x enhancerDataObject
#' @param activity numeric vector of activity values
#' @name predictExpression
#' @include enhancerDataObject-class.R
#' @export
methods::setMethod("predictExpression", signature = signature(x = "enhancerDataObject", activity="numeric"), function(x,activity) {
  if(x@linkFunction$type=="additive"){
    expression=activity
  } else if(x@linkFunction$type=="exponential"){
    expression=exp(activity)
  } else if(x@linkFunction$type=="logistic"){
    expression=x@linkFunction$value[length(x@linkFunction$value)]/(1+exp(-activity))
  } else {
    stop("Invalid link function")
  }
  return(expression)
})
