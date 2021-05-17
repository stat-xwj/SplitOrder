#' @title ladle
#' @description The bootstrap ladle method by Ye
#' @param x, predictors
#' @param y, response
#' @param nslices, parameter needed in the SIR methods
#' @param nboot, bootstrap 次数
#' @param method, AFM, pca, cca, ica, sir, save, dr, DEE, kernel sir.
#' @param order, parameter needed in the SIR methods
#' @export
ladle<-function(x=x,y=y,nslices=nslices,nboot=nboot,
                method=method,order=order)
{
  if(method=='AFM'){x<-y}
  p<-ncol(x)
  n<-nrow(x)
  z<-stand(x)
  if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
  if (missing(y)) {obs<-x} else {obs<-cbind(y,x)}
  if (!(missing(nslices))) {H <- nslices}
  if (missing(nboot)) {nboot <- as.integer(n/2)}

  if (method=="pca")
  {
    if (missing(order)) {nx <- x} else {
      nx<-x
      if (order>1)
      {for (j in 2:order) {nx<-cbind(nx,x^j)}}}
    mhat<-var(nx)
    obs<-nx
    p<-ncol(nx)
  }
  if (method=="cca")
  {
    if (is.vector(y)==TRUE) {y<-matrix(y,,1)}
    q<-ncol(y)
    if (q>10) {r<-min(r,as.integer(q/log(q)))}
    mhat<-mhat_cca(x,y)
  }
  if (method=="ica")
  {mhat<-mhat_ica(x)}
  if (method=="sir")
  {
    obs<-cbind(y,z)
    sdr<-dr(y~z,method=method,nslices=H,numdir=r)
    mhat<-sdr$evectors%*%diag(sdr$evalues)%*%t(sdr$evectors)
    r<-min(r,H-1)
  }
  if (method=="save")
  {
    obs<-cbind(y,z)
    sdr<-dr(y~z,method=method,nslices=H,numdir=r)
    mhat<-sdr$evectors%*%diag(sdr$evalues)%*%t(sdr$evectors)
  }
  if (method=="dr")
  {mhat<-mhat_dr(z,y,H)}
  if (method=="DEE")
  {mhat<-mhat_DEE(z,y)
  Bhat<-eigen(mhat)$vectors[,1:r]
  lam<-eigen(mhat)$values[1:(r+1)]
  }
  if (method=='AFM'){
    mhat<-mhat_AFM(y)
  }
  ################################################
  ## adjustment for NSDR to save computing time ##
  ################################################
  if (method!="kernel sir")
  {
    Bhat<-eigen(mhat)$vectors[,1:r]
    lam<-eigen(mhat)$values[1:(r+1)]
  } else {
    m <- mhat_ksir(x=x,y=y,nslices=H)
    r <- min(r,H-2) #only consider nonzero sample eigenvalues
    #    nx <- m$nx
    #    obs <- cbind(y,nx)
    #    Bhat <- m$evectors[,1:r]
    lam <- m$evalues[1:(r+1)]
    nx<-(m$nx)%*%m$evectors[,1:p]
    obs<-cbind(y,nx)
    Bhat<-diag(1,p)
  }
  fn0<-rep(0,r+1)
  for (j in 1:nboot)
  {
    u<-round(runif(n,min=-0.5,max=n+0.5))
    bs<-obs[u,]
    if (method=="pca")
    {mstar<-var(bs)}
    if (method=="cca")
    {mstar<-mhat_cca(bs[,-(1:q)],bs[,1:q])}
    if (method=="ica")
    {mstar<-mhat_ica(bs)}
    if ((method=="sir")|(method=="save"))
    {
      bsdr<-dr(bs[,1]~bs[,-1],method=method,nslices=H)
      vsdr<-qr.Q(qr(bsdr$evectors))
      mstar<-vsdr%*%diag(bsdr$evalues)%*%t(vsdr)
    }
    if (method=="dr")
    {mstar<-mhat_dr(stand(bs[,-1]),bs[,1],H)}
    if (method=="DEE")
    {mstar<-mhat_DEE(stand(bs[,-1]),bs[,1])}
    if (method=='AFM'){
      mstar<-mhat_AFM(y)
    }
    if (method!="kernel sir")
    {Bstar<-eigen(mstar)$vectors[,1:r]} else {
      Bstar<-mhat_ksir(phi=bs[,-1],y=bs[,1],nslices=H)$evectors[,1:r]}
    for (i in 1:r)
    {
      fn0[i+1]<-fn0[i+1]+1-abs(det(t(Bstar[,1:i])%*%Bhat[,1:i]))
    }
  }
  fn0<-fn0/nboot
  res<-list(0)
  res$fn0 <- fn0
  res$fn <-fn0 / (1+sum(fn0))
  res$lam <- lam
  res$phin <- lam/(1+sum(lam))
  res$gn <- res$fn + res$phin
  res$d <- which.min(res$gn) - 1
  return(res)
}





#' @title yeweiss
#' @description The original bootstrap ladle consider eigenvectors only
#' @param x, predictors
#' @param y, response
#' @param nslices, parameter needed in the SIR methods
#' @param nboot, bootstrap 次数
#' @param method, AFM, pca, cca, ica, sir, save, dr, DEE, kernel sir.
#' @param order, parameter needed in the SIR methods
#' @param delta is the scalar to select the threshold for "small" bootstrap variability
#' @export
yeweiss<-function(x=x,y=y,nslices=nslices,nboot=nboot,method=method,
                  order=order, delta=delta)
{
  ladel_ans <- ladle(x=x,y=y,nslices=nslices,nboot=nboot,method,order)
  a<-ladel_ans$fn0
  m<-max(a)
  res <- c()
  for (i in 1:length(delta)){
    thd<-delta[i]*m
    if (sum(a<thd)==0) {res <- which.min(a)} else {
      res[i] <- which((a[a<thd])[sum(a<thd)]==a)-1}
  }
  ans <- list(0)
  ans$ladel <- ladel_ans
  ans$yeweiss <- res
  # ans <- list(0)
  # ans$ladel <- ladel_ans$d
  # a <- a[-1]
  # ans$YW <- which.min(a)
  return(ans)
}
