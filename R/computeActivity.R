methods::setGeneric("computeActivity", function(x,...) {
  standardGeneric("computeActivity")
})

#' Computes the value of the activity term
#'
#' Computes the value of the activity function for each row of the designMatrix
#' @param x enhancerDataObject
#' @name computeActivity
#' @include enhancerDataObject-class.R
#' @examples
#'
#' @export
methods::setMethod("computeActivity", signature(x = "enhancerDataObject"), function(x) {
    return(as.numeric(x@designMatrix %*% as.matrix(x@linkFunction$value[1:ncol(x@designMatrix)])))
})
