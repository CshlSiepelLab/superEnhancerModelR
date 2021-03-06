methods::setGeneric("optimDE", function(x,maxit=100,refine=TRUE,threads=1,control=DEoptim::DEoptim.control()) {
  standardGeneric("optimDE")
})

#' Optimizes enhancer model with evolutionary algorithm
#'
#' Runs differential evolution optimization method on enhancer model
#' @param x enhancerDataObject
#' @param maxit the number of generations run (overrides the value set by control)
#' @param refine if TRUE, run gradient descent on the best solution from the evolutionary algorithm
#' @param threads integer value for number of threads to use (overrides the value set by control)
#' @param control A list of control parameters that can be passed to \code{\link[DEoptim]{DEoptim}} (also see \code{\link[DEoptim]{DEoptim.control}})
#' @name optimDE
#' @include enhancerDataObject-class.R
#' @examples
#' ## Create a test design data.frame
#' design=expand.grid(E1=c(0,1),E2=c(0,1),E3=c(0,1))
#' ## Create fake expression data
#' expression=c(0.001,0.2,0.3,0.7,0.4,0.8,0.6,1)*100
#' ## Create activity function
#' actFun=formula(~E1+E2+E3+E1:E2)
#' edo = enhancerDataObject(expression,design,actFun)
#' edo=optimDE(edo,maxit=500,refine=TRUE)
#'
#' @export
methods::setMethod("optimDE", signature(x = "enhancerDataObject"), function(x,maxit=100,refine=TRUE,threads=1,control=DEoptim::DEoptim.control()) {
  opt=list(fn = optimizeModel,object=x)
  opt[["lower"]]=c(x@linkFunction$constraints$lower,x@errorModel$constraints$lower)
  opt[["upper"]]=c(x@linkFunction$constraints$upper,x@errorModel$constraints$upper)
  ## Setup control options
  opt$control=control
  opt$control$itermax=maxit
  if(is.null(opt$control$NP) || is.na(opt$control$NP)){opt$control$NP=20*length(opt$lower)}
  ## Setup parallel cluster
  if (threads > 1) {
    if(!requireNamespace("parallel", quietly = TRUE)){
      stop("parallel not installed")
    } else {
      cl <- parallel::makeCluster(threads)
      opt[["control"]]$cluster=cl
    }
  }

  ## Run diff. evol. optimization
  best=do.call(DEoptim::DEoptim,args = opt)
  ## Set object parameters equal to parameters from best run
  x@linkFunction$value[names(x@linkFunction$value)]=best$optim$bestmem[names(x@linkFunction$value)]
  x@errorModel$value[names(x@errorModel$value)]=best$optim$bestmem[names(x@errorModel$value)]

  if(!is.null(opt[["control"]]$cluster))
    parallel::stopCluster(cl)

  if(refine){
    best=optimGD(x,restarts=1,refine=FALSE)
    x@linkFunction$value=best@linkFunction$value
    x@errorModel$value=best@errorModel$value
  }
  return(x)
})
