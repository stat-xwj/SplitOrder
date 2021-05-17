# The following function is given bu Luo Wei.

#' @title L-infinity norm
#' @description L-infinity norm
#' @param v, vector
#' @export
li<-function(v){
  return(max(abs(v)))
}




#' @title L-2 norm
#' @description L-2 norm
#' @param v, vector
#' @export
#'
lt<-function(v){
  return(sqrt(sum(v^2)))
}





#' @title cu
#' @description to accumulate a vector
#' @param v, vector
#' @export
cu<-function(v)
{
  p<-length(v)
  u<-v
  for (i in 1:p)
  {u[i]<-sum(v[1:i])}
  return(u)
}




#' @title dcu
#' @description to de-accumulate a vector
#' @param u, vector
#' @export
dcu<-function(u)
{
  p<-length(u)
  v<-u
  for (i in 2:p)
  {v[i]<-u[i]-u[i-1]}
  return(v)
}



#' @title mstand
#' @description standardize X componentwise
#' @param x, matrix
#' @export
mstand<-function(x)
{
  p<-ncol(x)
  z<-x
  for (i in 1:p) {z[,i]<-(x[,i]-mean(x[,i]))/sqrt(var(x[,i]))}
  return(z)
}





#' @title runis
#' @description generate RV from a uniform distribution
#' on a hyper-sphere centered at the origin and with radius sqrt p
#' @param n, items
#' @param p, dimension
#' @export
runis<-function(n,p)
{
  x<-matrix(rnorm(n*p),n,)
  xc<-x-rep(1,n)%*%t(apply(x,2,mean))
  s<-sqrt(c(apply(xc^2,1,sum)))
  xunis<-(xc/s)*sqrt(p)
  sunib<-(runif(n))^(1/p)*sqrt((p+2)/p)
  x<-xunis*sunib
  return(x)
}



#' @title agv
#' @description main function
#' @param x is n*p matrix
#' @param y is n-dim vector
#' @param k is dimension of S (augmentation predictor), default p/10
#' @param rep is number of repitations of augmentation, default 10
#' @param H is number of slices (for inverse regression)
#' @param method is the working matrix function
#' @param dist is the distribution of S, normal(default) or runis
#' @export
agv<-function(x,y,k,rep,H,method,dist)
{
  p<-ncol(x)
  n<-nrow(x)



  if (p<n) {z<-stand(x)} else {z<-x}
  if (missing(k)==TRUE) {k=as.integer(p/10)+1}
  if (missing(rep)==TRUE) {rep=10}
  if (missing(dist)==TRUE) {dist="normal"}

  val<-rep(0,p+1)
  vec<-rep(0,p+1)

  if(method == 'kernel sir'){
    r = min(p-1, H-2)
    val<-rep(0,r+1)
    vec<-rep(0,r+1)
  }

  if ((method=="pca")|(method=="sparse pca"))
  {
    sigm<-median(diag(var(x)))
  }
  for (j in 1:rep)
  {
    set.seed(j)
    if (dist=="normal") {ax<-matrix(rnorm(n*k),n,k)}
    if (dist=="unis")   {ax<-runis(n,k)}
    xs<-cbind(z,ax)
    if ((method=="sir")|(method=="save"))
    {
      ve<-dr(y~xs,method=method,nslices=H,numdir=p)$evectors[(p+1):(p+k),1:p]
    }
    if (method=="dr")
    {
      ve<-eigen(mhat_dr(xs,y,H))$vectors[(p+1):(p+k),1:p]
    }
    if (method=="sparse sir")
    {
      q<-min(H-1,p)
      vec<-vec[1:(q+1)]
      z<-mstand(x)
      xs<-cbind(z,ax)
      v<-LassoSIR(X=xs,Y=y,H=H,no.dim=q)
      vv<-v$eigen.value[1:q]
      vv<- vv - min(vv)
      act<-(apply((v$beta[,1:q])^2,1,sum)>0)
      pr<-sum(act[1:p])
      pf<-sum(act)
      rxs<-xs[,act]
      u<-eigen(dr(y~rxs,method="sir",nslices=H)$M)
      ve<-u$vectors[-(1:pr),1:q]
      ve<-matrix(ve,,q)
    }
    if (method=="ica")
    {
      ve<-eigen(mhat_ica(xs))$vectors[(p+1):(p+k),1:p]
    }
    if (method=="cca")
    {
      if (is.vector(y)==TRUE) {y<-matrix(y,,1)}
      ve<-eigen(mhat_cca(xs,y))$vectors[(p+1):(p+k),1:p]
    }
    if (method=="pca")
    {
      xs<-cbind(x,ax*sqrt(sigm))
      ve<-eigen(var(xs))$vectors[(p+1):(p+k),1:p]
    }
    if (method=="sparse pca")
    {
      dm<-as.integer(p^(1/3))
      vec<-vec[1:(dm+1)]
      val<-val[1:(dm+1)]
      xs<-cbind(x,ax*sqrt(sigm))
      ve<-spca(xs,K=dm,para=rep(0.1,dm),type="predictor",
               sparse="penalty")$loadings[(p+1):(p+k),]
    }
    # modified by xiawj
    if (method=='kernel sir')
    {
      ve <- mhat_ksir(x=xs,y=y,nslices=H)$evectors[(p+1):(p+k),1:p]
      # m <- mhat_ksir(x=xs,y=y,nslices=H)[(p+1):(p+k),]
      # r <- min(r,H-2) #only consider nonzero sample eigenvalues
      #    nx <- m$nx
      #    obs <- cbind(y,nx)
      #    Bhat <- m$evectors[,1:r]

    }
    if(method=='DEE'){
      ve = eigen(mhat_DEE(xs,y))$vectors[(p+1):(p+k),1:p]
    }
    ##
    if (method!="sparse sir") {ve<-matrix(ve,k,)}
    vec[-1]<-vec[-1] + (apply(ve^2,2,sum))
  }
  vec<-vec/rep

  if ((method=="sir")|(method=="save"))
  {val[1:p]<-dr(y~x,method=method,nslices=H,numdir=p)$evalues}
  if (method=="dr")
  {val[1:p]<-eigen(mhat_dr(x,y,H))$values}
  if (method=="sparse sir")
  {
    val=c(vv,0)
  }
  if (method=="ica")
  {val[1:p]<-sqrt(eigen(mhat_ica(x))$values)}
  if (method=="cca")
  {
    if (is.vector(y)==TRUE) {y<-matrix(y,,1)}
    val[1:p]<-eigen(mhat_cca(x,y))$values
  }
  if (method=="pca")
  {val[1:p]<-eigen(var(x))$values}
  if (method=="sparse pca")
  {
    dm<-as.integer(p^(1/3))
    val[1:dm]<-spca(x,K=dm,para=rep(0.1,dm),type="predictor",
                    sparse="penalty")$pev
    h<-eigen(var(x))$values[1]
    val[1:dm]<-val[1:dm]*h/val[1]
  }
  if (method=='kernel sir'){
    val[1:r] = mhat_ksir(x=x,y=y,nslices=H)$evalues[1:p]
  }
  if(method=='DEE'){
    val[1:p] = eigen(mhat_DEE(x,y))$values
  }
  res<-list()
  res$vectors<-vec
  res$values<-val
  res$d<-which.min(cu(vec)+val/(1+cu(val)))-1

  # if (method=='kernel sir'){
  #   res$d<-which.min((cu(vec)+val/(1+cu(val)))[1:r])-1
  # }

  res$k<-k
  return(res)
}

