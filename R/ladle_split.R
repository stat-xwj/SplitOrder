#' @title ladle_split
#' @description The main function in splitLadle method
#' @import MASS
#' @import mnormt
#' @import dr
#' @import mixAK
#' @import LassoSIR
#' @import elasticnet
#' @import stats
#' @param x
#' @param y
#' @param nslices, parameter needed in the SIR methods
#' @param nsplit, nsplit means we want to split the whole samples into nsplit parts
#' @param method, AFM, pca, cca, ica, sir, save, dr, DEE, kernel sir.
#' @param order, parameter needed in the SIR methods
#' @param criterion, criterion = 1 means the first B_k combination methods, means the total B_k dots each split part of B_ik, and sum
#' criterion = 2 means each split parts dots (even split is needed)
#' criterion = 3 means criterion1 + criterion2 and then average them
#' @param n0, n0 means to repeat split method n0 times
#' @export

ladle_split<-function(x=x,y=y,nslices=nslices,nsplit=nsplit,
                      method=method,order=order, criterion=criterion, n0=n0)


{
  if(method=='AFM'){x<-y}
  p<-ncol(x)
  n<-nrow(x)
  z<-stand(x)
  r = p-1
  # if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
  if (missing(y)) {obs<-x} else {obs<-cbind(y,x)}
  if (!(missing(nslices))) {H <- nslices}
  if (missing(nsplit)) {nsplit <- as.integer(n/2)}

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
  Fn0 <- rep(0,r+1)
  for (kn0 in 1:n0) {
    fn0<-rep(0,r+1)
    # Split <- createFolds(x[, 1], nboot)
    set.seed(kn0)
    U<-sample(1:n, n, replace = FALSE)
    block = floor(n/nsplit) # ????ȡ????????ȱʧֵ
    if(criterion==1){
      for (j in 0:(nsplit-1))
      {
        u<-U[(j*block+1):(j*block+block)]
        # u<-round(runif(n,min=-0.5,max=n+0.5))
        bs<-obs[u,]

        if (method!="kernel sir")
        {mstar <- Fmstar(bs=bs, method=method, q=q, H=H, y=y)
        Bstar<-eigen(mstar)$vectors[,1:r]} else {
          Bstar<-mhat_ksir(phi=bs[,-1],y=bs[,1],nslices=H)$evectors[,1:r]}
        for (i in 1:r)
        {
          fn0[i+1]<-fn0[i+1]+1-abs(det(t(Bstar[,1:i])%*%Bhat[,1:i]))
        }
      }
      fn0<-fn0/nsplit
    }else if(criterion==2){
      DBstar <- array(0,dim=c(p, r, nsplit))
      for (j in 0:(nsplit-1))
      {
        u<-U[(j*block+1):(j*block+block)]
        # u<-round(runif(n,min=-0.5,max=n+0.5))
        bs<-obs[u,]

        if (method!="kernel sir")
        {mstar <- Fmstar(bs=bs, method=method, q=q, H=H, y=y)
        DBstar[,,j+1]<-eigen(mstar)$vectors[,1:r]} else {
          DBstar[,,j+1]<-mhat_ksir(phi=bs[,-1],y=bs[,1],nslices=H)$evectors[,1:r]}
      }
      for (i in 1:r)
      {
        # MDBstar <- diag(i)
        # for (sp in seq(1,nsplit,2)) {
        #   MDBstar= t(DBstar[,1:i,sp]) %*% DBstar[,1:i,sp+1] %*% MDBstar
        # }
        # fn0[i+1] = 1-abs(det(MDBstar))
        flag = 0
        sp1=1
        while (sp1<=nsplit) {
          sp2 = sp1 + 1
          while(sp2<=nsplit){
            fn0[i+1] = fn0[i+1] + 1-abs(det(t(DBstar[,1:i,sp1]) %*% DBstar[,1:i,sp2]))
            flag = flag + 1
            sp2 = sp2 + 1
          }
          sp1 = sp1 + 1
        }
      }
      fn0<-fn0/(flag)
    }else if(criterion == 3){
      DBstar <- array(0,dim=c(p, r, nsplit+1))
      for (j in 0:(nsplit-1))
      {
        u<-U[(j*block+1):(j*block+block)]
        # u<-round(runif(n,min=-0.5,max=n+0.5))
        bs<-obs[u,]

        if (method!="kernel sir")
        {mstar <- Fmstar(bs=bs, method=method, q=q, H=H, y=y)
        DBstar[,,j+1]<-eigen(mstar)$vectors[,1:r]} else {
          DBstar[,,j+1]<-mhat_ksir(phi=bs[,-1],y=bs[,1],nslices=H)$evectors[,1:r]}
      }
      DBstar[,,nsplit+1] <- Bhat[,1:r]
      for (i in 1:r)
      {
        # MDBstar <- diag(i)
        # for (sp in seq(1,nsplit,2)) {
        #   MDBstar= t(DBstar[,1:i,sp]) %*% DBstar[,1:i,sp+1] %*% MDBstar
        # }
        # fn0[i+1] = 1-abs(det(MDBstar))
        flag = 0
        sp1=1
        while (sp1<=nsplit+1) {
          sp2 = sp1 + 1
          while(sp2<=nsplit+1){
            fn0[i+1] = fn0[i+1] + 1-abs(det(t(DBstar[,1:i,sp1]) %*% DBstar[,1:i,sp2]))
            flag = flag + 1
            sp2 = sp2 + 1
          }
          sp1 = sp1 + 1
        }
      }
      fn0<-fn0/(flag)

    }

    Fn0 <- Fn0+fn0
  }
  Fn0 <- Fn0/n0

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

  # res$rn = cu(Fn0/(sum(Fn0))) + phinR  # 没有截断
  # res$rn = res$fn + res$phinR

  res$rn = cu(res$fn) + res$phinR  # # 有截断


  res$an = cu(Fn0 / (1+sum(Fn0)))+res$lam/(sum(res$lam)+1)

  res$SSLE = which.min(res$gn) - 1
  res$SSRE = which.min(res$rn[-1])
  res$SSAE = which.min(res$an)-1

  res$r = r0



  return(res)
}
