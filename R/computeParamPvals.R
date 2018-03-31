methods::setGeneric("computeParamPvals", function(x,threads=1,...) {
  standardGeneric("computeParamPvals")
})

#' Computes per parameter p-value
#'
#' Computes a p-value for the intercept and each enhancer coefficient in the link function by performing
#' likelihood ratio tests between the full model and a reduced model without that coefficient.
#' Can be parallelized if doParallel is installed. NOTE. Make sure your full model is already fit before
#' running this function!!!
#' @param x enhancerDataObject
#' @param threads number of threads to use for refiting, requires doParallel
#' @param ... all aditional parameters to be passed to optimizer function \code{\link{optimDE}}
#' @name computeParamPvals
#' @include enhancerDataObject-class.R
#'
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
#' edo=optimDE(edo,maxit=300)
#' param.stats=computeParamPvals(edo,maxit=300)
#' @export
methods::setMethod("computeParamPvals", signature(x = "enhancerDataObject"), function(x,threads=1,...){
  full.ll=ll(x)
  full.formula=x@activityFunction
  ## Check for doParallel
  if (threads > 1 && requireNamespace("doParallel", quietly = TRUE)) {
    stop("Parallelization not yet implemented")
  } else {
    ## A list to hold all the reduced models
    rModelList=list()
    if(threads>1)
      warning("Multiple threads were specified but doParallel is not installed. Running in single threaded mode.")
    for(i in 1:ncol(x@designMatrix)){
      write(paste("Fitting reduced model without term",colnames(x@designMatrix)[i]),stdout())
      reduced.model=x
      if(colnames(x@designMatrix)[i]=="(Intercept)"){
        reduced.formula=update.formula(full.formula, formula(paste0(". ~ . -",colnames(x@designMatrix)[i])))
      } else {
        reduced.formula=update.formula(full.formula, formula(". ~ . +0"))
      }
      ## Replace function and drop column from design matrix
      reduced.model@activityFunction=reduced.formula
      reduced.model@designMatrix=reduced.model@designMatrix[,-i]
      ## Refit model and save
      rModelList[[colnames(x@designMatrix)[i]]]=optimDE(x,...)
    }
  }
  ## Compute log-likelihoods
  model.ll=unlist(lapply(rModelList,ll))
  ## Create output table, D statistic is computed where the reduced model is the null
  out=data.frame(parameter=names(rModelList),ll=as.numeric(model.ll),D=2*(full.ll-model.ll))
  out$p.val=dchisq(out$D,df=1)
  row.names(out)=NULL
  return(list(pVals=out,reducedModels=rModelList))
})