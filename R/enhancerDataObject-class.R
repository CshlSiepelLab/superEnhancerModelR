## Validity checker for enhancerDataObject
validityCheck <- function(object){
  errors=c()
  if(length(object@expressionData)!=nrow(object@designMatrix))
    errors=c("There must be the same number of expression values as row in the design matrix",errors)
  if(any(object@expressionData<=0))
    errors=c(errors,"Zero values in expression data are not permitted, add a very small value")
  if (length(errors) == 0) TRUE else errors
}

#' Class enhancerDataObject
#'
#' Class \code{enhancerDataObject} holds knockout and expression data along with model parameters
#'
#' @name enhancerDataObject-class
#' @rdname enhancerDataObject-class
#' @exportClass enhancerDataObject
methods::setClass(Class = "enhancerDataObject",
                  representation = representation(expressionData = "numeric", designMatrix = "matrix",activityFunction="formula",
                                                    errorModel="list",linkFunction="list"),
                  validity = validityCheck)

#' enhancerDataObject Constructor
#'
#' Contructs a simple kinetic model object from synthesis and degredation rates
#' @param expressionData A numeric vector with expression data, one entry for each row in the designInfo matrix
#' @param designInfo A 0/1 (absence/presence) coded matrix with one column for each enhancer and one row for each experiment
#' @param activityFunction the formula for the activity function (e.g. ~ E1 + E2 + E1*E2)
#' @param errorModel Desired error model for the regression function. Currently only lognormal or gaussian.
#' @param linkFunction Function that converts activity to gene expression
#' @param activityParameterBounds sets upper and lower bounds for coefficients of activity function - c(lower,upper)
#' @param errorParameterBounds sets upper and lower bounds for the error parameter - c(lower,upper)
#' @param scaleParameterBounds sets upper and lower bounds for the saturation parameter (only in logistic model) - c(lower,upper)
#' @name enhancerDataObject
#' @examples
#' ## Create a test design data.frame
#' design=expand.grid(E1=c(0,1),E2=c(0,1),E3=c(0,1))
#' ## Create fake expression data
#' expression=c(0.001,0.2,0.3,0.7,0.4,0.8,0.6,1)*100
#' ## Create activity function
#' actFun=formula(~E1+E2+E3+E1:E2)
#' edo = enhancerDataObject(expression,design,actFun)
#'
#' @export
enhancerDataObject <- function(expressionData,designInfo,activityFunction,errorModel=c("lognormal","gaussian"),
                               linkFunction=c("additive","exponential","logistic"),activityParameterBounds=c(-10^3,10^3),
                               errorParameterBounds=c(10^-3,10^3),scaleParameterBounds=c(10^-3,10^3)){
  ## Makes sure that activity function is a string
  if(!inherits(activityFunction,"formula"))
    stop("activityFormula must be a formula (e.g. \"~E1+E2+E1:E2\")")
  if(!is.data.frame(designInfo))
    stop("designInfo must be a data.frame")
  ## Analyze the formula object to ensure that all variable names are column names
  if(!all(all.vars(activityFunction) %in% colnames(designInfo)))
     stop("activityFunction formula contains variables not in design matrix")

  ## Create design matrix
  designMatrix=model.matrix(activityFunction,data=designInfo)
  nActVars=ncol(designMatrix) ## get number of activity parameters in model

  ## Handle error model
  emTypes=c("lognormal","gaussian")
  if(length(errorModel)>1)
    errorModel=errorModel[1]
  if(!errorModel %in% emTypes)
    stop(paste("Invalid error model type. Valid types are:",paste(emTypes,collapse = ", ")))
  if(length(errorParameterBounds)!=2)
    stop("errorParameterBounds must be a vector of length 2 with form c(lower,upper)")
  ## Build error model list
  emList=list(type=errorModel, value=c(sd=rngErrorParameters(1, errorParameterBounds)),
              constraints=list(lower=c(sd=errorParameterBounds[1]),upper=c(sd=errorParameterBounds[2])))

  ## Handle link function
  lfTypes=c("additive","exponential","logistic")
  if(length(linkFunction)>1)
    linkFunction=linkFunction[1]
  if(!linkFunction %in% lfTypes)
    stop(paste("Invalid link function type. Valid types are:",paste(lfTypes,collapse = ", ")))
  if(length(activityParameterBounds)!=2)
    stop("activityParameterBounds must be a vector of length 2 with form c(lower,upper)")
  ## Throw a warning if the the model selected is linear and the error model is log-normal and
  ## set all the lower bounds of all parameters to 10^-3
  if(linkFunction=="additive" && errorModel=="lognormal"){
    warning("The additive link function can produce negative expression values which are are outside the support of the log-normal function.
            Setting the lower bound of all coefficients to 10^-3 to prevent this.")
    activityParameterBounds[1]=10^-3
  }
  ## Initialize the parameters, biased towards the center of the range
  lfList=list(type=linkFunction, value=rngLinkParameters(nActVars,activityParameterBounds),
              constraints=list(lower=rep(activityParameterBounds[1],nActVars),upper=rep(activityParameterBounds[2],nActVars)))
  names(lfList$value)=colnames(designMatrix)
  names(lfList$constraints$lower)=colnames(designMatrix)
  names(lfList$constraints$upper)=colnames(designMatrix)

  ## Add parameter for logistic model
  if(linkFunction=="logistic"){
    lfList$value=c(lfList$value,scale=runif(1,scaleParameterBounds[1],scaleParameterBounds[2]))
    lfList$constraints$lower=c(lfList$constraints$lower,scale=scaleParameterBounds[1])
    lfList$constraints$upper=c(lfList$constraints$upper,scale=scaleParameterBounds[2])
  }

  ## Check for that there are sufficient parameters to fit the model of interest
  if(length(lfList$value)>nrow(unique(designInfo))){
    warning("The specified model has more parameters as experimental conditions! The model will not be identifiable.\n")
  }

  new("enhancerDataObject",expressionData=expressionData,designMatrix=designMatrix,activityFunction=activityFunction,
      errorModel=emList,linkFunction=lfList)
}
