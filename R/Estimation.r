########### Estimation
## 0. Center_X - not visible in the help
## 1. Stand_Multi
## 2. Esti_Skew_Mardia
## 3. Esti_Kurt_Mardia
## 4.  Esti_SkewVec_Mori
## 5.  Esti_SkewInd_Mori
## 6.  Esti_Skew_Mori_2
## 7.  Esti_Kurt_Mori_2
## 8. Esti_H3
## 9. Esti_H4
## 10. Hermite_Nth
## 11. Skew_Mardia
## 12. Kurt_Mardia
## 13. Estimates_MMom_MCum
## 14 Esti_MVSK
## 15 Variance_of_Esti_Skew
## 16 Esti_Variance_MVSK
## !!!!!!!!!!!!!!!
##
## Variance_of .....
## !!!!!!!!!!!!!!!!!

##  SkewEsti
##  SkewKron
##  SkewKronT
##  Variance_of
##################

# Centering a sample of vector variates,
#  @param x matrix of  sample, rows are observations of a d variate,
#  sample size is the number of rows
#  @return data matrix with rows centered by the sample means of columns
#  @family Standard
#
# @export
Center_X <- function(x) {
  apply(x, 2, function(y) y - mean(y))
  # this is: scale(y, center = TRUE, scale = FALSE)
}


#' Standardize multivariate data
#'
#' For data formed by d-variate vectors x with sample covariance S and sample mean M,
#' it computes the values
#' \eqn{z=S^{-1/2}(x-M)}
#'
#' @param x a multivariate data matrix, sample size is the number of rows
#' @return z multivarate data with null mean vector and
#' identity sample covariance matrix
#' @examples
#' x<-rmvnorm(1000,mean=c(0,0,1,3))
#' z<-Stand_Multi(x)
#' mu_z<- apply(z,2,mean)
#' cov_z<- cov(z)
#'
#' @family Standard
#' @export
Stand_Multi<-function(x){
  # x is a multivariate vector of data
  z<-Center_X(x)
  cx<-cov(x)
  svdx<-svd(cx)
  sm12<-svdx$u%*%diag(1/sqrt(svdx$d))%*%t(svdx$u)
  z1<-t(sm12%*%t(z))
  return(z1)
}


#' Mardia Skewness index
#'
#' Compute the multivariate Mardia's skewness index and
#' provides the p-value for the hypothesis of zero symmetry under the
#' Gaussian assumption
#' @param x A vector of multivariate data
#' @return The skewness index
#' @return The p-value
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.1
#' @family Estimation
#' @export
Esti_Skew_Mardia<-function(x){
  z<-Stand_Multi(x)
  n=dim(z)[1]
  d=dim(z)[2]
  MSkew<-sum((z%*%t(z))^3)/n^2
  pval<-pchisq(n*MSkew/6,choose(d+2,3),lower.tail = FALSE)
  return(list("Mardia.Skewness"=MSkew,"p.value"=pval))
}

#' Mardia Kurtosis Index
#'
#' @param x A vector of multivariate data
#' @return The kurtosis index
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.1
Esti_Kurt_Mardia<-function(x){
  z<-Stand_Multi(x)
  n=dim(z)[1]
  MK<-sum(t(z^2)%*%z^2)/n
  return(list("Mardia.Kurtosis"=MK))
}




#' Mori's skewness vector
#'
#' @param x A vector of multivariate data
#' @return The skewness vector
#'
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.2
#'
#' @export
Esti_SkewVec_Mori<-function(x){
  z<-Stand_Multi(x)
  z2<-apply(z^2,1,sum)
  MS<-apply(z2*z,2,mean)
  return(list("Mori.Skewness.Vector"=MS))
}

#' Mori's skewness index
#'
#' @param x A vector of multivariate data
#' @return The skewness index
#' @return The p-value for the hypothesis of zero symmetry under the
#' Gaussian assumption
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.2
#'
#' @export
Esti_SkewInd_Mori<-function(x){
  z<-Stand_Multi(x)
  n=dim(z)[1]
  d=dim(z)[2]
  z2<-apply(z^2,1,sum)
  MS<-apply(z2*z,2,mean)
  MS<-sum(MS^2)
  pval<-pchisq(n*MS/(2*(d+2)),d,lower.tail = FALSE)
  return(list("Mori.Skewness.Index"=MS,"p.value"=pval))
}

#' Mori's kurtosis matrix
#'
#' @param x A vector of multivariate data
#' @return The kurtosis matrix
#' @family Indexes
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Example 6.9
#'
#' @export
Esti_Kurt_Mori<-function(x){
  z<-Stand_Multi(x)
  n=dim(z)[1]
  d<-dim(z)[2]
  MM<-array(0,c(d,d))
  for (i in 1: n){
    temp<-z[i,]%*% t(z[i,])%*%z[i,]%*% t(z[i,])
    MM<-MM+temp
  }
  MK<-MM/n -(d+2)*diag(d)

  return(list("Mori.Kurtosis"=MK))
}


#' Mori's skewness vector
#'
#' Computes Mori's skewness vector given the skewness vector.
#' It can be used to compute estimated or population values.
#'
#' @param H3 the skewness vector.
#' @return Mori's skewness vector
#' @examples
#' #Compute Mori's skewness vector for the multivariate Skew-Normal
#' alpha<-c(10,5)
#' omega<-diag(rep(1,2))
#' H3SN<-MVSK_SkewNorm_Th(omega, alpha)$SkewX
#' Skew_Mori(H3SN)
#' @family Indexes
#' @references  S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644, Example 2.

Skew_Mori<-function(H3){
  d<-(length(H3))^(1/3)
  MoS<-kronecker(t(as.matrix(c(diag(d)))),diag(d))%*%H3
  MoS<-as.vector(MoS)
  return(MoS)
}

#' Mori's kurtosis matrix
#'
#'  Computes Mori's kurtosis matrix given the kurtosis vector.
#' It can be used to compute estimated or population values.
#'
#' @param H4 the d-variate kurtosis vector.
#' @return Mori's kurtosis matrix
#' @family Indexes
#' @examples
#' #Compute Mori's kurtosis matrix for the multivariate Skew-Normal
#' alpha<-c(10,5)
#' omega<-diag(rep(1,2))
#' H4SN<-MVSK_SkewNorm_Th(omega, alpha)$KurtX
#' Kurt_Mori(H4SN)
#' @references  S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644, Example 9.
Kurt_Mori<-function(H4){
  d<-(length(H4))^(1/4)
  MoK<-kronecker(diag(d^2),t(as.matrix(c(diag(d)))))%*%H4
  MoK<-matrix(MoK,d,d)
  return(MoK)
}



#' Estimate the third d-variate polynomial
#'
#' The vector x is standardized and the skewness vector is computed
#' @param x a d-variate data vector
#'
#' @examples
#' x<-mvtnorm::rmvnorm(100,rep(0,3))
#' H3<-Esti_H3(x)
#' @export
Esti_H3<-function(x){
  z<-Stand_Multi(x)
  H3<-apply(apply(z,1,Hermite_Third),1,mean)
  return(H3)
}

#' Estimate the fourth d-variate polynomial
#'
#' The vector x is standardized and the kurtosis vector is computed
#'
#' @param x a d-variate data vector
#' @examples
#' x<-mvtnorm::rmvnorm(100,rep(0,3))
#' H4<-Esti_H4(x)
Esti_H4<-function(x){
  z<-Stand_Multi(x)
  H4<-apply(apply(z,1,Hermite_Fourth),1,mean)
  return(H4)
}



#' Mardia Skewness index
#'
#' Computes Mardia's skewness index given the skewness vector.
#' It can be used to compute estimated or population values.
#'
#' @param H3 the skewness vector.
#' @return Mardia's skewness index
#' @family Indexes
#' @examples
#' #Compute Mardia's skewness index for the multivariate Skew-Normal
#' alpha<-c(10,5)
#' omega<-diag(rep(1,2))
#' H3SN<-MVSK_SkewNorm_Th(omega, alpha)$SkewX
#' Skew_Mardia(H3SN)
#'
#'
#' @references  S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644, Example 1.
Skew_Mardia <-function(H3){
  MS<-sum(H3^2)
  return(MS)
}

#' Mardia Kurtosis index
#'
#' Computes Mardia's kurtosis index given the kurtosis vector.
#' It can be used to compute estimated or population values.
#'
#' @param H4 the kurtosis vector.
#' @return Mardia's kurtosis index
#' @family Indexes
#' @examples
#' #Compute Mardia's kurtosis index for the multivariate Skew-Normal
#' alpha<-c(10,5)
#' omega<-diag(rep(1,2))
#' H4SN<-MVSK_SkewNorm_Th(omega, alpha)$KurtX
#' Kurt_Mardia(H4SN)
#' @references  S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644, Example 7.

Kurt_Mardia<-function(H4){
  d<-(length(H4))^(1/4)
  t(c(diag(d^2)))%*%H4+d*(d+2)
}


##############################
#' Estimate for Multivariate Moments and Cumulants
#'
#' Provides estimates of multivariate moments and cumulants up to order k.
#' By default data are standardized, only centering can be used
#'
#' @param X d-vector data
#' @param k The highest moment order
#' @param centr set to 1 (and sc_sigM=0) if only centering is needed
#' @param sc_sigM set to 1 (and centr=0) if standardization of multivariate data is needed
#' @family Estimation
#' @export
Esti_MMom_MCum  <- function(X,k,centr=0,sc_sigM=1){
  if (dim(X)[1]<dim(X)[2]) (stop("X is not proper data matrix"))
  Mu_X <-  colMeans(X)  # column
  Vari_X <-cov(X)
  if (centr==1) X <- Center_X(X)
  if (sc_sigM==1 ) X <- Stand_Multi(X)
  Mu_X_k <- NULL
  Mu_X_k[[1]] <- Mu_X
  for (j in c(2:k)) {
    Xadk <- apply(X, 1, function(y) KronPower(y,j))
    Mu_X_k[[j]] <- rowMeans(Xadk) #apply(Xadk, 2,mean)
  }
  Cum_X_k <- Mom2CumMulti(Mu_X_k)
  EstiMomCum <- list(estMu = Mu_X,estVar=Vari_X,
                     estMu.k=Mu_X_k,estCum.k=Cum_X_k)
  return(EstiMomCum)
}

####################
#' Estimates for multivariate Skewness and Kurtosis
#'
#' Provides estimates of mean, variance, skewness and kurtosis vectors for  d-variate data
#' @param X d-variate data vector
#' @return The list of the estimated quantities
#' @examples
#' x<-mvtnorm::rmvnorm(100,rep(0,d),sigma = 3*diag(rep(1,d)))
#' MVSK<-Esti_MVSK(x)
#' names(MVSK)
#' MVSK$estSkew
#' @family Estimation
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.
#' @export
Esti_MVSK <- function(X){
  kim <- Esti_MMom_MCum(X,4,centr=0,sc_sig=1)
  estiMVSK <- list(kim[[1]],kim[[2]], kim[[4]][3],kim[[4]][4] )
  names(estiMVSK) <- c("estMu" ,   "estVar"  , "estSkew" , "estCurt")
  return(estiMVSK)
}


############
#' Asymptotic Variance for estimated skewness
#'
#' @param cum The theoretical/estimated cumulants up to order 6 in vector form
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. skewness (6.13)
#' @family Estimation
#' @export
Variance_of_Esti_Skew_Th <- function(cum){
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j3
  d <- length(cum[[1]])
  type_2j3 <- c(0,0,2,0,0,0)
  loc_type_2j3  <-  Partition_Type_eL_Location(type_2j3)
  L_2j3 <- Commutator_Moment_L(type_2j3,loc_type_2j3[1],loc_type_2j3[2],d)
  # %%%%%%%%%%%%%%%%%%%%%%%  Commutator  L_1j2_1j1
  # type_1j2_1j1  <- [1,1,0];%
  # loc_type_1j2_1j1   <-  PartitionType_eL_Location(type_1j2_1j1 );
  # L_1j2_1j1   <-  MomentCommutators_L(type_1j2_1j1 ,loc_type_1j2_1j1 (1),loc_type_1j2_1j1 (2),d);
  Id3 <-diag(d^3)
  K2 <-  kronecker( Commutator_Kmn(d,d),diag(d)) #
  K3  <- Commutator_Kmn(d,d^2);# L1_2_1_1 X van elõl
  Sym22 <-  Id3+K2+K3; #  Sym22 - L_1j2_1j1'
  #%%%%%%%%%%%%%%%% Commutator L2_H4
  # L2_H4i  <-  kronecker(L_1j2_1j1,L_1j2_1j1)*PermutingMatrixkroneckerProd([1 3 4 2 5  6 ],d);
  L2_H4i  <-  kronecker(t(Sym22),t(Sym22))%*%Commutator_Kperm(c(1,3,4,2,5,6 ),d)

  #%%%%%%%%%%%%%%%5  Commutator M_3 inverse already!!!!
  d1 <- c(d,d,d)
  M3_m_ni  <- t(Commutator_Mixing( d1,d1))
  #%%%%%%%%%%%%%%%%%%%%%
  Id <- diag(d)
  vec_Var_Skew  <- cum[[ 6]] + t(L_2j3)%*%kronecker(cum[[ 3]],cum[[ 3]]) +
    L2_H4i%*%kronecker(as.vector(Id),cum [[4]])+
    M3_m_ni%*%kronecker(as.vector(Id), kronecker(as.vector(Id),as.vector(Id)))-
    kronecker(cum[[ 3]],cum[[ 3]])
  Var_Skew  <-  matrix(vec_Var_Skew, nrow=d^3)
  return(Var_Skew)
}

####################################

#' Estimated Variance of  skewness and kurtosis
#'
#' Provides the estimated covariance matrices of the data-estimated skewness and kurtosis
#' vectors.
#'
#' @param X d-variate data
#' @return The list ov covaraince matrices of the skewness and kurtosis vectors
#' @examples
#' x<-mvtnorm::rmvnorm(100,rep(0,d),sigma = 3*diag(rep(1,d)))
#' MVSK<-Esti_MVSK(x)
#' EV_MVSK<-Esti_Variance_MVSK(x)
#' names(EV_MVSK)
#' EV_MVSK$Vari_Skew_e
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Chapter 6, formulae (6.13) and (6.22).
#' @family Estimation
#' @export
Esti_Variance_MVSK<- function(X){
  estiMVSK <- Esti_MVSK(X) #
  Z <-  Stand_Multi(X)
  H3Zt<- apply(Z,1,Hermite_Third) # Hermite
  #
  cH3 <-  -  apply(H3Zt,2,function(U) U-estiMVSK$estSkew[[1]])
  # cov(t(H3Zt))
  Vari_Skew_e <- cov(t(cH3))
  #############
  H4Zt <-  apply(Z,1,Hermite_Fourth)
  # Est_H4(Z)
  cH4 <-  -  apply(H4Zt,2,function(U) U-estiMVSK$estCurt[[1]])
  Vari_Kurt_e <- cov(t(cH4))
  esti.var.SK <- list(Vari_Skew_e, Vari_Kurt_e)
  names(esti.var.SK) <- c("Vari_Skew_e" ,   "Vari_Kurt_e")
  return(esti.var.SK)
}

###################

#' Asymptotic Variance of Estimated  kurtosis
#'
#' @param cum The theoretical/estimated cumulants up to the 8th order in vector form
#'
#' @references Gy.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. (6.26)
#' @family Estimation
#' @export

Variance_of_Esti_Kurt_Th <- function(cum){
  ##
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_1j3_1j5
  d <- length(cum[[1]])
  # kron2 <- function(Mm) kronecker(Mm,Mm)
  Id <- diag(d)
  vId <- as.vector(Id)
  #%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j3
  type_2j3 <- c(0,0,2,0,0,0)
  loc_type_2j3  <-  Partition_Type_eL_Location(type_2j3)
  L_2j3 <- Commutator_Moment_L(type_2j3,loc_type_2j3[1],loc_type_2j3[2],d)
  ###########L_1j3_1j5
  type_1j3_1j5 <- c(0,0,1,0,1,0,0,0)
  loc_type_1j3_1j5  <-  Partition_Type_eL_Location(type_1j3_1j5)
  L_1j3_1j5 <- Commutator_Moment_L(type_1j3_1j5,loc_type_1j3_1j5[1],loc_type_1j3_1j5[2],d)
  # %%%%%%%%%%%%%%%%%%%%%%%  Commutator  L_2j4
  type_2j4 <- c(0,0,0,2,0,0,0,0)
  loc_type_2j4  <-  Partition_Type_eL_Location(type_2j4)
  L_2j4 <- Commutator_Moment_L(type_2j4,loc_type_2j4[1],loc_type_2j4[2],d)
  #######
  ############## L2_H6i
  # L_1j1_1j3 available in canonical for only i.e. L_1j3_1j1, hence direct construction
  L_1j1_1j3 <-  diag(d^4) + Commutator_Kperm(c(2, 1, 3, 4),d)+
    Commutator_Kperm(c(3, 1, 2, 4) ,d)+
    Commutator_Kperm(c( 4, 1, 2, 3),d)
  Legy <- kronecker(L_1j1_1j3,L_1j1_1j3)
  L2_H6i=Legy*Commutator_Kperm(c(1, 3, 4,  5, 2, 6, 7, 8),d)
  #%%%%%%%%%%%%%%%5  Commutator M_4_m_n
  d1 <- rep(d,4)
  M4_m_ni  <- t(Commutator_Mixing( d1,d1))
  #%%%%%%%%%%%%%%%%%%%%% L_1j3_1j1
  type_1j3_1j1=c(1,0,1,0)
  loc_type_1j3_1j1  <- Partition_Type_eL_Location(type_1j3_1j1)
  L_1j3_1j1=Commutator_Moment_L(type_1j3_1j1,loc_type_1j3_1j1[1],loc_type_1j3_1j1[2],d)
  ################
  L22_H4m<- L22_H4(d)

  krId2 <- kron2(vId)
  # %% kurtosis (6.26)
  vec_Var_Kurt =cum[[ 8]] +
    t(L_1j3_1j5)%*%kronecker(cum[[5]],cum[[3]])+
    t(L_2j4)%*%kron2(cum[[4]]) - as.vector(kron2(cum[[4]])) +
    L2_H6i%*%kronecker(vId,cum[[ 6]]+
                         t(L_2j3)%*%kron2(cum[[3]]))+
    t(L22_H4m)%*%kronecker(krId2,cum[[4]])+
    M4_m_ni%*%kron2(krId2)+
    kron2(t(L_1j3_1j1))%*%kronecker(kron2(cum[[3]]),vId)

  Var_Kurt=matrix(vec_Var_Kurt, nrow= d^4);
  return(Var_Kurt)
}
