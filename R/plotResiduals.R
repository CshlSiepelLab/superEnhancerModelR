methods::setGeneric("plotResiduals", function(x,...) {
  standardGeneric("plotResiduals")
})

#' Plots REsidual
#'
#' Plots model residuals
#' @param x enhancerDataObject
#' @name optimMod
#' @include enhancerDataObject-class.R
#' @examples
#'
#' @export
methods::setMethod("plotResiduals", signature(x = "enhancerDataObject"),plotResiduals <- function(x){
  ## Get names of active enhancers for each experiment
  indiv.enhancers=grep(x = colnames(x@designMatrix),pattern = ")|:",invert = TRUE)
  rname=apply(x@designMatrix,1, function(z) paste(colnames(x@designMatrix)[indiv.enhancers][as.logical(z[indiv.enhancers])],collapse = "/"))
  rname[rname==""]="None"

  ## Get factor ordering of rows
  temp=data.frame(rname=rname,nenh=rowSums(x@designMatrix),stringsAsFactors = FALSE)
  lvls=unique(with(temp,temp[order(nenh,rname),])$rname)

  ## Get predicted expression
  predE=predictExpression(x)

  ## Get activity values for actual oberservations
  act=computeActivity(x)
  real=data.frame(activity=act,residual=x@expressionData-predE,enhancers=rname,stringsAsFactors = FALSE)
  real$enhancers=factor(real$enhancers,levels=lvls)

  ## Compute expression values for intermediate activity to get smooth curve
  sim.act=seq(min(real$activity)-abs(min(real$activity))*0.1,max(real$activity)+abs(max(real$activity))*0.1,by=0.01)
  out=data.frame(activity=sim.act,expression=predictExpression(x,sim.act))

  err=errorIntervals(x,sim.act,quantiles = 0.1)
  err$y=err$y-c(out$expression,rev(out$expression))

  cc <- scales::seq_gradient_pal("light blue", "blue", "Lab")(seq(0,0.5,length.out=length(unique(err$Quantile))))

  g=ggplot2::ggplot()+
    ggplot2::geom_polygon(data=err,ggplot2::aes(x=x,y=y,fill=Quantile),alpha=1)+
    ggplot2::scale_fill_manual(values=cc)+
    ggplot2::theme_bw(base_size = 28)+
    ggplot2::xlab("Activity/-Energy")+
    ggplot2::ylab("Expression")+
    ggplot2::geom_hline(yintercept = 0,linetype=3)+
    ggplot2::geom_point(data=real,ggplot2::aes(activity,residual,color=enhancers))


  return(g)
})
