methods::setGeneric("optimDE", function(x,maxit,refine=TRUE,threads=1) {
  standardGeneric("optimDE")
})

#' Optimizes enhancer model with evolutionary algorithm
#'
#' Runs differential evolution optimization method on enhancer model
#' @param x enhancerDataObject
#' @param maxit the number of generations run
#' @param refine if TRUE, run gradient descent on the best solution from the evolutionary algorithm
#' @param threads integer value for number of threads to use
#' @name optimDE
#' @include enhancerDataObject-class.R
#' @examples
#' ## Create a test design matrix
#' design = matrix(c(0,1,0,1,0,0,1,1),nrow=4)
#' colnames(design) = c("E1","E2")
#' design = as.data.frame(design)
#' ## Create fake expression data
#' expression = c(0,0.2,0.3,0.9)
#' ## Create activity function
#' actFun = formula(~E1+E2+E1:E2)
#' edo = enhancerDataObject(expression,design,actFun)
#' edo=optimDE(edo,refine=TRUE)
#'
#' @export
methods::setMethod("optimDE", signature(x = "enhancerDataObject"), function(x,maxit=100,refine=TRUE,threads=1) {
  opt=list(fn = optimizeModel,object=x)
  opt[["lower"]]=c(x@linkFunction$constraints$lower,x@errorModel$constraints$lower)
  opt[["upper"]]=c(x@linkFunction$constraints$upper,x@errorModel$constraints$upper)
  opt$control=list(itermax=maxit,NP=20*length(opt$lower),trace=100)
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
