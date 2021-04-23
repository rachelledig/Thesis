normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
}

centered <-function(y) {
  x<-y[!is.na(y)]
  x<-(x - mean(x)) / (sd(x)*2)
  y[!is.na(y)]<-x
  return(y)
}