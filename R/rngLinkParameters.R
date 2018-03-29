rngLinkParameters <- function(nActVars,activityParameterBounds){
  lMean=mean(activityParameterBounds)
  lSd=(activityParameterBounds[2]-activityParameterBounds[1])/40
  vals=rnorm(nActVars,mean = lMean,sd = lSd)
  vals[vals>activityParameterBounds[2]]=activityParameterBounds[2]
  vals[vals<activityParameterBounds[1]]=activityParameterBounds[1]
  return(vals)
}
