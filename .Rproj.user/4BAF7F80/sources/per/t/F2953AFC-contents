#' @title ladle_weight
#' @description ladle weight methods
#' @param x, predictors
#' @param y, response
#' @param nslices, parameter needed in the SIR methods
#' @param nsplit, nsplit means we want to split the whole samples into nsplit parts
#' @param method, AFM, pca, cca, ica, sir, save, dr, DEE, kernel sir.
#' @param order, used in SIR
#' @param wt, the given weight type 'normal','bernoulli'
#' @export

ladle_weight<-function(x=x,y=y,nslices=nslices,s=s, method=method,order=order, wt = wt)
{
  if(method=='AFM'){x<-y}
  p<-ncol(x)
  n<-nrow(x)
  z<-stand(x)
  r <- p-1
  # if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
  if (missing(y)) {obs<-x} else {obs<-cbind(y,x)}
  if (!(missing(nslices))) {H <- nslices}
  if (missing(s)) {s=10}

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
    # if (q>10) {r<-min(r,as.integer(q/log(q)))}
    mhat<-mhat_cca(x,y)
  }
  if (method=="ica")
  {mhat<-mhat_ica(x)}
  if (method=="sir")
  {
    obs<-cbind(y,z)
    sdr<-dr(y~z,method=method,nslices=H,numdir=r)
    mhat<-sdr$evectors%*%diag(sdr$evalues)%*%t(sdr$evectors)
    # r<-min(r,H-1)
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
    # r <- min(r,H-2) #only consider nonzero sample eigenvalues
    #    nx <- m$nx
    #    obs <- cbind(y,nx)
    #    Bhat <- m$evectors[,1:r]
    lam <- m$evalues[1:(r+1)]
    nx<-(m$nx)%*%m$evectors[,1:p]
    obs<-cbind(y,nx)
    Bhat<-diag(1,p)
  }
  fn0<-rep(0,r+1)
  a = (1-sqrt(5))/2
  pa = (1+sqrt(5))/(2*sqrt(5))
  for (j in 1:s)
  {
    if(wt == 'normal'){
      weight_vector = rnorm(nrow(obs))
    } else if(wt == 'bernoulli'){
      weight_vector = c()
      for (i in 1:nrow(obs)) {
        weight_vector[i] = sample(c(a,1-a), size = 1, prob = c(pa, 1-pa))
      }
    }

    bs = obs * weight_vector

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

  fn0<-fn0/s
  res<-list(0)

  Fn0 <- fn0

  ####################
  phinR <- rep(0, r+1)
  for (i in 1:(length(lam)-1)) {
    phinR[i+1]=lam[i+1]/lam[i]
  }

  if (p > 10) {r0<-as.integer(p/log(p))} else {r0 <- p-1}

  res<-list(0)
  res$fn0 <- Fn0[1:(r0+1)]
  res$fn <-res$fn0 / (1+sum(res$fn0))
  res$lam <- lam
  res$phin = lam[1:(r0+1)]/(1+sum(lam[1:(r0+1)]))
  res$phinR <- phinR[1:(r0+1)]

  res$gn <- res$fn + res$phin
  res$rn = cu(res$fn) + res$phinR  # # 有截断
  res$an = cu(Fn0 / (1+sum(Fn0)))+res$lam/(sum(res$lam)+1)

  res$SSLE = which.min(res$gn) - 1
  res$SSRE = which.min(res$rn[-1])
  res$SSAE = which.min(res$an)-1

  res$r = r0


  return(res)
}
