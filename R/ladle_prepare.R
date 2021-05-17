

#' @title stand
#' @description standardize x
#' @param x, a matrix
#' @return standardized x (z)
#' @export
stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/(n-1)
  eva <- eigen(sigma)$values
  # 由于精度原因，会造成一些很小的值取到负号，在这里直接取反。
  eva[eva < 0] = 0.00000001
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}





#' @title matpower
#' @description  power of a matrix
#' @param a, a matrix
#' @param alpha, a power
#' @export
matpower = function(a,alpha)
{
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))
}




#' @title mppower
#' @description  Moore-Penrose type power
#' @param matrix, a matrix
#' @param power, a number
#' @param ignore, criterion=ignore
#' @export
mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power,m)%*%
    t(evec[,1:m])
  return(tmp)
}





#' @title center
#' @description  center X (n*p matrix)
#' @param x, a matrix
#' @export
center = function(x){
  return(t(t(x)-apply(x,2,mean)))}







#' @title mhat_cca
#' @description  Mhat for CCA
#' @param x, a matrix
#' @param y, a matrix
#' @export
mhat_cca<-function(x,y)
{
  n<-nrow(x)
  p<-ncol(x)
  q<-ncol(y)
  nsx<-matpower(cov(x),-0.5)
  nvy<-matpower(cov(y),-1)
  cxy<-cov(x,y)
  M<-nsx%*%cxy%*%nvy%*%t(cxy)%*%nsx
  return(M)
}







#' @title mhat_ica
#' @description  Mhat for ICA
#' @param x, a matrix
#' @export
mhat_ica<-function(x)
{
  n<-nrow(x)
  p<-ncol(x)
  z<-stand(x)
  w<-apply(z^2,1,sum)
  m<-t(z)%*%diag(w)%*%z/n - (p+2)*diag(p)
  #m<-m%*%m
  return(m%*%m)
}






#' @title diffin
#' @description  vector transformation needed in ST for CCA
#' @param v, a matrix
#' @export
diffin<-function(v)
{
  s<-length(v)
  if (s==1) {return(1)} else {m<-(v%*%t(rep(1,s))-rep(1,s)%*%t(v))
  r<-numeric(0)
  for (i in 1:(s-1)) {r<-c(r,m[i,(i+1):s])}
  return(r)}
}




#' @title diffbw
#' @description some transformation
#' @param u, vector
#' @param v, scaler
#' @export
diffbw<-function(u,v)
{
  s<-length(u)
  r<-numeric(0)
  for (i in 1:s)
  {r<-c(r,u[i]-v)}
  return(r)
}






#' @title slicing
#' @description discretizing y into H slices
#' @param y, a vector or matrix
#' @param H, a slice
#' @export
slicing<-function(y,H)
{
  n = length(y)
  if (length(levels(as.factor(y)))>H)
  {
    ytilde<-rep(0,H+1)
    ytilde[1]<-min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-quantile(y,h/H)
    }
  }
  if (length(levels(as.factor(y)))<=H)
  {
    H <- length(levels(as.factor(y)))
    ytilde<-rep(0,H+1)
    ytilde[1]=min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-min(y[y>ytilde[h]])
    }
  }
  ytilde[H+1]=max(y)+1
  prop<-rep(1,H)
  for (i in 1:H)
  {
    prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  return(res)
}




#' @title mhat_dr
#' @description mhat for directional regression
#' @param x, a matrix
#' @param y, a vector or matrix
#' @param H, a slice
#' @export
mhat_dr<-function(x,y,H){
  n<-nrow(x)
  p<-ncol(x)
  z<-stand(x)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  ind<-matrix(0,n,H)
  zbar<-matrix(0,p,H)
  for (j in 1:(H-1))
  {
    ind[,j]<-((y >= ytilde[j])&(y < ytilde[j+1]))
    zbar[,j]<- (t(z)%*%(ind[,j]))/sum(ind[,j])
  }
  ind[,H]<-(y >= ytilde[H])
  zbar[,H]<- (t(z)%*%(ind[,H]))/sum(ind[,H])
  A<-matrix(0,p,p)
  B<-matrix(0,p,p)
  C<-0
  for (q in 1:H)
  {
    Z<-(t(z))[,ind[,q]==1]-zbar[,q]
    A<-A + prop[q]*((Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q]))%*%
                      (Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q])) - diag(1,p))
    B<-B + sqrt(prop[j])*(zbar[,q]%*%t(zbar[,q]))
    C<-C + sqrt(prop[j])*(t(zbar[,q])%*%zbar[,q])
  }
  C<-as.vector(C)
  M<-2*A + 2*(B%*%B) + 2*B*C
  return(M)
}


#########################################
###  mhat for DEE  ###
#########################################

# mhat_DEE <- function(x,y){
# # x the observed predictor matrix
# # y the observed response vector
# # K the dimension of central subspace
# # This code is applicable when the response y is multivariate.
#  size=dim(x)
#  n=size[1]
#  p=size[2]
# xscore=x-matrix(1,n,1)%*%apply(x,2,mean)
# Indexy=order(y)
# Mn=matrix(0,p,p)
# for(i in 2:(n-2)){
#   p1=i/n;
#   p0=(n-i)/n;
#   mean1=apply(xscore[Indexy[1:i],],2,mean)
#   mean0=apply(xscore[Indexy[(i+1):n],],2,mean)
#   Mn=Mn+p1*mean1%*%t(mean1)+p0*mean0%*%t(mean0)
# }
# Mn=Mn/(n-2)
# return(Mn)
# }





#' @title mhat_DEE
#' @description   mhat for DEE
#' @param x, a matrix
#' @param y, a vector or matrix
#' @param K, a slice
#' @export
mhat_DEE <- function(x, y, K){
  # x the observed predictor matrix
  # y the observed response vector
  # K the dimension of central subspace
  # This code is applicable when the response y is multivariate.


  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  if(missing(K)){K <- p}
  x_orig <- x
  x <- x-matrix(1,n,1)%*%apply(x,2,mean)
  Lambda <- array(0, dim = c(p, p, n))
  for (ii in 1:n) {
    address <- as.matrix(as.integer((y[, ]<y[ii, ])))
    # address <- matrix(0, n, q)
    # for (kk in 1:q) {
    #   address[, kk] <-(y[, kk]<y[ii, kk])
    # }
    unit <- address[!duplicated(address),]
    if(!is.matrix(unit)){
      unit <- as.matrix(unit)
    }
    x_mean <- matrix(0, nrow(unit), p)
    for(jj in 1:nrow(unit)){
      temp <- x[diag((address-matrix(1, n, 1)*unit[jj, ]) %*% t(address-matrix(1, n, 1)*unit[jj, ]))==0,]
      if (is.vector(temp)){
        x_mean[jj, ]<-matrix(0, 1, p)}else{
          x_mean[jj, ] <- apply(temp, 2, mean)
        }
    }
    x_mean <- x_mean[apply(x_mean^2, 1, sum)!=0,]
    x_mean <- as.matrix(x_mean)
    Lambda[,,ii] <- cov(x_mean)
  }

  mean_Lambda <- matrix(0, p, p)
  for(i in 1:n){
    mean_Lambda = mean_Lambda + Lambda[,,i]
  }
  mean_Lambda = mean_Lambda/n
  return(mean_Lambda)
  # Eigen_vector<-eigen(mean_Lambda)
  #
  # res <- list(0)
  # res$mhat <- mean_Lambda
  # res$lambda <- Eigen_vector$values
  # res$vectors <- Eigen_vector$vectors
  #
  # res$direction <- M_12(cov(x_orig))%*% Eigen_vector$vectors[,1:K]
  # return(res)
}



#' @title mhat_AFM
#' @description   mhat for Approximate factor models
#' @param Y, a matrix
#' @export
mhat_AFM <- function(Y){
  p <- ncol(Y)
  n <- nrow(Y)
  return(Y%*%t(Y)/(n*p))
}




################################
# caulate the M-1/2 ######
#' @title M_12
#' @description  caulate the M-1/2
#' @param Sigma, a matrix
#' @export
M_12 <-function(Sigma){
  lamda <- solve(eigen(Sigma)$vectors)%*%Sigma%*%(eigen(Sigma)$vectors)
  lamda_sqrt <- matrix(0,nrow = nrow(lamda),ncol = nrow(lamda))
  diag(lamda_sqrt) <- sqrt(diag(lamda))
  Sigma_sqrt <- (eigen(Sigma)$vectors)%*%lamda_sqrt%*%solve(eigen(Sigma)$vectors)
  return(solve(Sigma_sqrt))
}






#' @title mhat_ksir
#' @description mhat for kernel SIR
#' @param x, the original predictor
#' @param y, the original response
#' @param phi, the transformed predictor -in kernel space
#' @param Sigma, a matrix
#' @param nslices, slices number
#' @export
mhat_ksir<-function(x=x,phi=phi,y=y,nslices=nslices)
{
  n<-length(y)
  H<-nslices
  if ((missing(phi)))
  {
    p<-ncol(x)
    phi=cbind(1,x)
    for(i in 1:p)
    {phi = cbind(phi,x[,i]*x[,1:i])}
  }
  oy = order(y)   # ordered y
  ophi = phi[oy,] # ordered phi
  bslice = H        # number of slices; b means beween
  wslice = n/bslice  # number of observations in each slice: w means within
  ephiy=numeric()    # E(Y|slice)
  for(i in 1:bslice){
    ephiy=rbind(ephiy,apply(ophi[((i-1)*wslice+1):(i*wslice),],2,mean))}
  ephi = apply(phi,2,mean)
  m<-var(t(t(ephiy)-ephi))
  vd<-svd(t(t(ephiy)-ephi))
  v<-vd$v[,order((vd$d)^2,decreasing=TRUE)]
  vc<-svd(diag(1,ncol(phi))-v%*%t(v))$u
  lam<-(vd$d)^2/(H-1)
  res<-list()
  res$mhat<-m
  res$nx<-phi
  res$evectors<-cbind(v,vc)
  res$evalues<-lam
  return(res)
}





#' @title pick
#' @description   pick d from p-values of ST
#' @param v, vector
#' @param alpha, confidence terms
#' @export
pick<-function(v,alpha){
  p<-length(v)
  m=0
  while ((m<p)&(v[m+1]<alpha))
  {m=m+1}
  return(m)
}





#' @title Fmstar
#' @description function used in YW abnd ladel
#' @param bs,
#' @param method,
#' @param q,
#' @param H,
#' @param y,
#' @export
Fmstar <-function(bs=bs, method=method, q=q, H=H, y=y){
  p = ncol(bs)-1
  if (method=="pca")
  {mstar<-var(bs)}
  if (method=="cca")
  {mstar<-mhat_cca(bs[,-(1:q)],bs[,1:q])}
  if (method=="ica")
  {mstar<-mhat_ica(bs)}
  if ((method=="sir")|(method=="save"))
  {
    # 这里修改添加了一个 stand()
    bsdr<-dr(bs[,1]~bs[, -1],method=method,nslices=H)
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
  return(mstar)
}

###################
###     end     ###
###################





