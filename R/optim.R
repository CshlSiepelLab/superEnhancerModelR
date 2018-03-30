methods::setGeneric("optimGD", function(x,...) {
  standardGeneric("optimGD")
})

#' Optimizes enhancer model
#'
#' Runs optimization method on enhancer model
#' @param x enhancerDataObject
#' @name optimGD
#' @include enhancerDataObject-class.R
#' @examples
#' ## Create a test design matrix
#' design = matrix(c(0,1,0,1,0,0,1,1),nrow=4)
#' colnames(design) = c("E1","E2")
#' design = as.data.frame(design)
#' ## Create fake expression data
#' expression = c(0,0.2,0.3,0.9)
#' ## Create activity function
#' actFun = formula(~E1+E2+E1*E2)
#' edo = enhancerDataObject(expression,design,actFun)
#' edo=optimGD(edo,restarts=50)
#'
#' @export
methods::setMethod("optimGD", signature(x = "enhancerDataObject"), function(x,restarts=200,method=NULL,maxit=NULL,refine=FALSE) {
  ## Set up options to be passed to optimizer
  opt=list(par=c(x@linkFunction$value,x@errorModel$value),fn = optimizeModel,object=x)
  if (!is.null(method)){
    opt[["method"]]=method
  }  else {
    opt[["method"]]="L-BFGS-B"
  }
  if(opt[["method"]]=="L-BFGS-B"){
    opt[["lower"]]=c(x@linkFunction$constraints$lower,x@errorModel$constraints$lower)
    opt[["upper"]]=c(x@linkFunction$constraints$upper,x@errorModel$constraints$upper)
  }

  if(!is.null(maxit)){
    opt[["control"]]=list(maxit=maxit)
  }
  ## Run the optimizer of choice
  best=do.call("optim",args = opt)
  x@linkFunction$value[names(x@linkFunction$value)]=best$par[names(x@linkFunction$value)]
  x@errorModel$value[names(x@errorModel$value)]=best$par[names(x@errorModel$value)]

  if(restarts>1){
    for(j in 2:restarts){
      lpar=apply(cbind(x@linkFunction$constraints$lower,x@linkFunction$constraints$upper),1,function(x) rngLinkParameters(1,x))
      epar=apply(cbind(x@errorModel$constraints$lower,x@errorModel$constraints$upper),1,function(x) rngErrorParameters(1,x))
      opt[["par"]]=c(lpar,epar)
      mod=do.call(optim,args = opt)
      if(mod$value<best$value){
        best=mod
      }
    }
  }

  x@linkFunction$value[names(x@linkFunction$value)]=best$par[names(x@linkFunction$value)]
  x@errorModel$value[names(x@errorModel$value)]=best$par[names(x@errorModel$value)]
  ## If refine, run one additional gradient descent
  if(refine){
    opt.ref=opt
    opt.ref$par=best$par
    best=do.call(optim,args = opt.ref)
  }
  ## Return best parameter values to object
  x@linkFunction$value[names(x@linkFunction$value)]=best$par[names(x@linkFunction$value)]
  x@errorModel$value[names(x@errorModel$value)]=best$par[names(x@errorModel$value)]
  return(x)
})
