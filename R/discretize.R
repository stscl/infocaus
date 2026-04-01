
discretize = \(x, n = 5, method = 'natural', large = 3000, prop = 0.15, 
               seed = 123456789, thr = 0.4, iter = 100, bps = NULL, right_closed = TRUE){
  if (any(inherits(x,'factor'),inherits(x,'character'))){
    return(as.integer(as.factor(x)))
  }

  return(RcppDisc(x,n,method,large,prop,seed,thr,iter,bps,right_closed))
}
