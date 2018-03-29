myLog<-function(x){
  x[x<0]=-Inf
  x[x==0]=-300
  x[x>0]=log(x[x>0])
  return(x)
}
