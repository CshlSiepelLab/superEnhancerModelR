methods::setGeneric("optimDE", function(x,...) {
  standardGeneric("optimDE")
})

#' Optimizes enhancer model with evolutionary algorithm
#'
#' Runs differential evolution optimization method on enhancer model
#' @param x enhancerDataObject
#' @name optimMod
#' @include enhancerDataObject-class.R
#' @examples
#'
#' @export
methods::setMethod("optimDE", signature(x = "enhancerDataObject"), function(x,maxit=100,refine=FALSE) {
  opt=list(fn = optimizeModel,object=x)
  opt[["lower"]]=c(x@linkFunction$constraints$lower,x@errorModel$constraints$lower)
  opt[["upper"]]=c(x@linkFunction$constraints$upper,x@errorModel$constraints$upper)
  opt$control=list(itermax=maxit,NP=20*length(opt$lower),trace=100)

  ## Run diff. evol. optimization
  best=do.call(DEoptim::DEoptim,args = opt)
  ## Set object parameters equal to parameters from best run
  x@linkFunction$value[names(x@linkFunction$value)]=best$optim$bestmem[names(x@linkFunction$value)]
  x@errorModel$value[names(x@errorModel$value)]=best$optim$bestmem[names(x@errorModel$value)]

  if(refine){
    best=optimGD(x,restarts=1,refine=FALSE)
    x@linkFunction$value[names(x@linkFunction$value)]=best$par[names(x@linkFunction$value)]
    x@errorModel$value[names(x@errorModel$value)]=best$par[names(x@errorModel$value)]
  }
  return(x)
})
