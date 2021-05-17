
#' @title RRE_RE_BIC
#' @description Xia et al(2015), Luo et al(2009), Zhu et al(2006,2010)
#' @param x, predictors
#' @param y, response
#' @param method, mtthod choices
#' @param nslices, used in SIR
#' @export

RRE_RE_BIC <- function(x=x, y=y, method=method, nslices=nslices){
  if(method=='AFM'){
    p <- nrow(y)
    n <- ncol(y)
  }else{
    p <- ncol(x)
    n <- nrow(x)}

  # ?˴???r?Ǹ???ladel?????ж????ģ???????Zhu(2020)??û?п??????????壬???ǽ?r=p
  # if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
  r <- p-1

  eigvalue <- matrix(0, p, 1)

  c <- log(n)/(10*sqrt(n))
  if(method == 'AFM'){
    c <- log(n)/(10*n)
    mhat <- mhat_AFM(y)}

  if(method=='DEE'){
    mhat <- mhat_DEE(x, y)

  }
  if(method=='sir'){
    if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
    z <- stand(x)
    obs<-cbind(y,z)
    sdr<-dr(y~z,method=method,nslices=nslices,numdir=r)
    mhat<-sdr$evectors%*%diag(sdr$evalues)%*%t(sdr$evectors)
    r<-min(r,nslices-1)
  }
  if (method=="dr")
  {
    if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
    z <- stand(x)
    mhat<-mhat_dr(z,y,nslices)}
  eigval<-eigen(mhat)$values[1:(r+1)]
  lambda <- eigval
  # lambdanew <- eigval/(1+eigval)
  lambdanew <- eigval
  d1 <- length(lambda)
  # if (p > 10) {dmax<-as.integer(p/log(p))} else {dmax <- p-1}
  lambda1 <- (lambdanew[2:d1]+c)/(lambdanew[1:(d1-1)]+c)

  # dmax <- min(d1, 10)
  # RE
  dmax <- d1
  if(method=='AFM'){dmax <- min(p, n)/2}

  # lambda2 <- (lambda[1:(dmax-1)])/(lambda[2:dmax])
  lambda2 <- (lambda[2:dmax])/(lambda[1:(dmax-1)])
  dRRE <- which.min(lambda1)
  # dRE <- which.max(lambda2)
  dRE <- which.min(lambda2)
  #------------BIC--------------
  L <- log(eigval+1) - eigval
  denominator <- 2*sum(L)
  G <- matrix(0, p, 1)
  for (k in 1:p) {
    if(method!='AFM'){G[k] <- n*sum(L[1:k])/denominator -sqrt(n)*k*(k+1)/p}
    else{G[k] <- n*sum(L[1:k])/denominator -log(n)*k*(k+1)/p}

  }
  dBIC <- which.max(G)

  resu <- list(0)
  resu$RRE_lam = lambda1
  resu$RE_lam = lambda2
  resu$BIC_lam = G
  resu$RRE <- dRRE
  resu$RE <- dRE
  resu$BIC <- dBIC
  return(resu)
}
