

####################################################
########## zhu et.al(2020) #########################
############## TDRR ################################

#' @title TDRR
#' @description zhu et.al(2020)
#' @param x, predictors
#' @param y, response
#' @param method, mtthod choices
#' @param nslices, used in SIR
#' @param tao, tunning parameters
#' @export

TDRR <- function(x=x, y=y, method = method, nslices=nslices, tao){
  if(missing(tao)){tao<-0.5}
  if(method=='AFM'){
    n <- nrow(y)
    p <- ncol(y)
    m <- min(p, n)
    c1 <- log(m)/(10*sqrt(m))
    c2 <- log(m)/(5*sqrt(m))
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
    c1 <- 0.1*log(n)/sqrt(n)
    c2 <- 0.2*log(n)/sqrt(n)
    # c1 <- 2/(sqrt(n)*log(n))
    # c2 <- 0.3/(n^(1/4))
  }

  structure1 <- matrix(0, p-1, 1)
  if(method=='AFM'){
    mhat <- mhat_AFM(y)
  }
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

  eigval<-eigen(mhat)$values
  lambda11 <- eigval/(1+eigval)
  if(method!='AFM'){lambda11 <- lambda11^2}
  d1 <- length(lambda11)
  lambda12=(lambda11[1:(d1-1)]+c1)/(lambda11[2:d1]+c1)-1;
  lambda13=(lambda12[2:(d1-1)]+c2)/(lambda12[1:(d1-2)]+c2);
  res <- c()
  for(j in seq(1, length(tao), 1)){
    lambda14=lambda13<tao[j]
    if(sum(lambda14) == 0){res[j] <- 0}else{
      for (i in seq(length(lambda14), 1, -1)) {
        if(lambda14[i] == 1){
          res[j] <- i
          break}
      }
    }
  }


  return(res)
}
