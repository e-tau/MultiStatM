#' Hermite_Coeff
#' n the Hermite polynomial Hen(x) of variance 1, 
#' the absolute value of the #'coefficient of \eqn{x^k} 
#' is the number of (unordered) partitions of an n-element
#' set into \eqn{k} singletons and \eqn{(n ??? k)/2}   (unordered) pairs.
#' coefficient of \eqn{x^k} in  \eqn{H_N(x)}
#' 
#' @family Hermite
#' @export
Hermite_Coeff<- function(N){ 
  
PTA<-Partition_Type_All(N)
el_j<-PTA$eL_r
S_m_j<-PTA$S_r_j

Hcoeff<-rep(0,ceiling(N/2+1-(N%%2)))
kk=0
  for (k in 0:N) {
      if (N%%2== k%%2) {
        el=c(k,(N-k)/2,rep(0,N-2))
        loc_type_el<-c(0,0)
        for (m in 1:length(el_j)){
          if (is.vector(el_j[[m]])){
            if (prod((el==el_j[[m]]))) {loc_type_el<-c(m,1)}
          }
         else {
           for (mm in 1:dim(el_j[[m]])[1])
             if (prod((el==el_j[[m]][mm,]))) {loc_type_el<-c(m,mm)}
         } 
        }
      kk=kk+1;
      Hcoeff[kk] = (-1)^((N-k)/2)*S_m_j[[loc_type_el[1]]][loc_type_el[2]]
  }
 }
Hcoeff=Hcoeff[ceiling((N/2)+1-(N%%2)):1]
return(Hcoeff)
}


#' Hermite_Poly_HN
#' % (4.24) p.157
#' given x and sigma2=sigma^2
#' output values of  Hermite polynomials from 1 to N 
#' @family Hermite
#' @export
Hermite_Poly_HN<-function(x,N,sigma2=1){
H_N_x<-rep(0,N)
H_N_x[1]<-x
for (n in 2:N){
  Hcoeff<-Hermite_Coeff(n)
  nH<-length(Hcoeff)
  powersX<-seq(n,0,by=-2)
  powersSigma2<-(n-powersX)/2
  x_powers=rep(x,nH)^powersX
  sigma2_powers<-rep(sigma2,nH)^powersSigma2
  H_N_x[n]<-sum(Hcoeff*x_powers*sigma2_powers)
}
return(H_N_x)
}

#' Hermite_Poly_HN_Inv
#' % (4.23) p.157
#' given Hermite polynomials  H_N_x from 1 to N and sigma2=sigma^2
#' output powers of x_val: x^n, n=1:N
#' 
#' @family Hermite



Hermite_Poly_NH_Inv<-function(H_N_x,sigma2=1){
N<-length(H_N_x)
x_val<-rep(0,N)
x_val[1]=H_N_x[1]
for (n in 2:N) {
  Hcoeff<-Hermite_Coeff(n)
  nH<-length(Hcoeff)
  powersX<-seq(n,0,by=-2)
  if ((n%%2==0)) {H_N_x1 =rev(c(1,H_N_x[seq(2,n,by=2)]))}
  else {H_N_x1 =rev(H_N_x[seq(1,n,by=2)])}

powersSigma2<-(n-powersX)/2
signCoeff<-rep(-1,nH)^powersSigma2
Xcoeff<-Hcoeff*signCoeff
sigma2_powers<-rep(sigma2,nH)^powersSigma2
x_val[n]<-sum(Xcoeff*H_N_x1*sigma2_powers)
}
return(x_val)
}

#' Multivariate Hermite coefficients
#' 
#' 
#' @family Hermite
#' @export
Hermite_CoeffMulti<-function(N,d){
PTA<-Partition_Type_All(N)
el_j<-PTA$eL_r
HcoeffMatrix<-vector(mode = "list", length = ceiling(N/2+1-(N%%2)))
kk=0
for (k in 0:N) {
  if (N%%2== k%%2) {
    el=c(k,(N-k)/2,rep(0,N-2))
    loc_type_el<-c(0,0)
    for (m in 1:N) {
      if (is.vector(el_j[[m]])){
        if (prod((el==el_j[[m]]))) {loc_type_el<-c(m,1)}
      }
      else {
        for (mm in 1:dim(el_j[[m]])[1])
          if (prod((el==el_j[[m]][mm,]))) {loc_type_el<-c(m,mm)}
      } 
    
    }
    kk=kk+1;
    HcoeffMatrix[[kk]]<- (-1)^((N-k)/2)*t(Commutator_Moment_L(el,loc_type_el[2],loc_type_el[1],d))
  }
}
HcoeffMatrix<-HcoeffMatrix[seq(ceiling(N/2+1-(N%%2)),1,by=-1)]
return(HcoeffMatrix)
}

#' Multivariate Hermite Polynomials
#' 
#' 
#' @family Hermite
#' @examples 
#' x<-c(1,3)  
#' N<-3
#' Sig2<-matrix(c(1,0,0,1),2,2,byrow = T)
#' Hermite_Poly_HN_Multi(x,N,Sig2)
#' @export
Hermite_Poly_HN_Multi<-function(x,Sig2,N){

d=length(x)
H_N_x<-vector(mode = "list", length = N)
x_ad<- vector(mode="list",length=N)
x_ad[[1]]=x
if (N>1){
for (n in 2:N){
  x_ad[[n]]<-kronecker(x_ad[[n-1]],x)
  }
}

vSig2<-c(Sig2)
vSig2_ad<-vector(mode="list",length=ceiling(N/2+1-(N%%2)))
vSig2_ad[[1]]<-1
if (N>1){
  for (n in 2:ceiling(N/2+1-(N%%2))){
    vSig2_ad[[n]]<-kronecker(vSig2_ad[[n-1]],vSig2)
  }
}

H_N_x[1]<-x_ad[1]
if (N>1) {
  for (n in 2:N){
    HcoeffMatrix<-Hermite_CoeffMulti(n,d)
    nH=length(HcoeffMatrix)
    X_powers<-vector(mode="list",length=nH)
    
    if ((n%%2)==0){
      X_powers[1:(nH-1)]<-x_ad[seq(n,1,by=-2)]
      X_powers[[nH]]<-1
    }
    else {X_powers<-x_ad[seq(n,1,by=-2)]}
    
    Sigma2_powers<-vSig2_ad[1:nH]
    H_N_x0=0
    for (k in 1:nH){
      H_N_x0<-H_N_x0+HcoeffMatrix[[k]]%*%kronecker(Sigma2_powers[[k]],X_powers[[k]])
    }
    H_N_x[[n]]=H_N_x0  
  }
}

return(H_N_x)
}

#' Inverse d-variate Hermite Polynomial
#' 
#' @family Hermite
#' 
#' @examples
#' Sig2=matrix(c(1,0,0,1),2,2,byrow=T)
#' x<-c(1,3)
#' N<-4
#' H_N_X<-Hermite_Poly_HN_Multi(x,sigma2,N)
#' Hermite_Poly_NH_Multi_Inv(H_N_X,Sig2,N)

Hermite_Poly_NH_Multi_Inv<-function(H_N_X,Sig2,N) {
  
d<-length(H_N_X[[1]])  
x_ad_val<-vector(mode="list",length=N)
vSig2<-c(Sig2)
vSig2_ad<-vector(mode="list",length = ceiling(N/2+1-(N%%2)))
vSig2_ad[[1]]<-1
if (N>1) {
  for (n in 2:ceiling(N/2+1-(N%%2))){
    vSig2_ad[[n]]<-kronecker(vSig2_ad[[n-1]],vSig2)
  }
}
x_ad_val[1]<-H_N_X[1]

if (N>1){
  for (n in 2:N) {
    HcoeffMatrix<-Hermite_CoeffMulti(n,d)
    nH<-length(HcoeffMatrix)
    H_N_Xp<-vector(mode="list",length=nH)
    if ((n%%2==0)){
      H_N_Xp[1:(nH-1)]<-H_N_X[seq(n,1,by=-2)]
      H_N_Xp[[nH]]<-1
    }
    else {H_N_Xp <- H_N_X[seq(n,1,by=-2)]}
    Sigma2_powers<-vSig2_ad[1:nH]
    X0<-0
    for (k in 1:nH) {
      X0<-X0+(-1)^(k-1)*HcoeffMatrix[[k]]%*%kronecker(Sigma2_powers[[k]],H_N_Xp[[k]])
    }
    x_ad_val[[n]]<-X0
  }
}
return(x_ad_val)

}
