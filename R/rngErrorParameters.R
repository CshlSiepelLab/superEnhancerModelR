rngErrorParameters <- function(nVars,errorParameterBounds){
  eMean=errorParameterBounds[1]+(errorParameterBounds[2]-errorParameterBounds[1])*0.05
  eVar=(errorParameterBounds[2]-errorParameterBounds[1])/100
  shape=(eMean^2)/eVar
  scale=eMean/shape
  vals=rgamma(nVars,shape=shape,scale=scale)
  vals[vals>errorParameterBounds[2]]=errorParameterBounds[2]
  vals[vals<errorParameterBounds[1]]=errorParameterBounds[1]
  return(vals)
}
