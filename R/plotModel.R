methods::setGeneric("plotModel", function(x,...) {
  standardGeneric("plotModel")
})

#' Plots Model
#'
#' Plots data and model prediction
#' @param x enhancerDataObject
#' @name optimMod
#' @include enhancerDataObject-class.R
#' @examples
#'
#' @export
methods::setMethod("plotModel", signature(x = "enhancerDataObject"),plotModel <- function(x){
  ## Get activity values for actual oberservations
  act=computeActivity(x)
  real=data.frame(activity=act,observed=x@expressionData)

  sim.act=seq(min(real$activity)-abs(min(real$activity))*0.1,max(real$activity)+abs(max(real$activity))*0.1,by=0.01)
  out=data.frame(activity=sim.act,expression=predictExpression(x,sim.act))

  err=errorIntervals(x,sim.act,quantiles = 0.1)

  cc <- scales::seq_gradient_pal("light blue", "blue", "Lab")(seq(0,0.5,length.out=length(unique(err$Quantile))))

  g=ggplot2::ggplot()+
    ggplot2::geom_polygon(data=err,ggplot2::aes(x=x,y=y,fill=Quantile),alpha=1)+
    ggplot2::scale_fill_manual(values=cc)+
    ggplot2::geom_line(data=out,ggplot2::aes(activity,expression),color="black",size=2)+
    ggplot2::geom_point(data=real,ggplot2::aes(activity,observed),color="red")+
    ggplot2::theme_bw(base_size = 28)+
    ggplot2::xlab("Activity/-Energy")+
    ggplot2::ylab("Expression")+
    ggplot2::theme(legend.position = c(0, 1),legend.justification = c(0, 1))

  return(g)
})
