## 1. Rand_Uni_d
## 2. Rand_SkewNorm
## 3. Rand_CFUSN
## 4. Rand_CFUSSD
## 5.
## 6.
## 7. MVSK_UniDistr_Th
## 8. MomCum_UniDistr_Th
## 9. MVSK_UniAbsDistr_Th_
## 10. MVSK_SkewNorm_Th
## 11. MomCum_Zabs_Th
## 12. MomCum_MZabs_Th  (Multi)
## 13. MomCum_CFUSN_Th


#####################
#########################################################
#########################    #############
##########################################################
#' Random Uniform on the sphere
#'
#' Generate random d-vectors from the Uniform distribution on the sphere
#' @param n sample size
#' @param d dimension
#'
#' @return  A random matrix \eqn{n x d}
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Random generation
#' @export

Rand_Uni_d <-  function(n,d){
Z <- matrix( rnorm(d*n), nrow = n)
U3 <- t(apply(Z, 1, function(y) y/sqrt(sum(y^2))))
return(U3)
}
#  rowSums(U3^2)

########################################################
###########      ###################
#######################################################
#' Random Multivariate Skew Normal
#'
#' Generate random d-vectors  from the multivariate Skew Normal distribution
#' @param n sample size
#' @param Om correlation matrix   with d dimension
#' @param alf shape parameter vector of dimension d
#' @return  A random matrix \eqn{n x d}
#'
#' @examples
#' alpha<-c(10,5,0)
#' omega<-diag(rep(1,3))
#' x<-Rand_SkewNorm(20,omega,alpha)
#'
#' @references Azzalini, A. with the collaboration of Capitanio, A. (2014).
#' The Skew-Normal and Related Families. Cambridge University Press,
#'   IMS Monographs series.
#' @references Gy.H.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Section 5.1.2
#'  @family Random generation
#' @export
Rand_SkewNorm <- function(n,Om,alf){

  d <- length(alf)
  delt<- as.vector((t(Om%*%alf)))*(1/as.vector(sqrt(1+t(alf)%*%Om%*%alf)))
  SI <- rbind(c(1,t(delt)),cbind(delt,Om))
  ieg<- eigen(SI)
  V <- ieg$vectors
  Sig <- V %*% diag(sqrt(ieg$values)) %*% t(V)
  Z <- matrix(rnorm((d+1)*n),nrow=n)
  Z <- (apply(Z, 1, function(y) Sig%*%y))
  X3 <- t(Z[2:(d+1),])*
    matrix(kronecker(rep(1,d),Z[1,]>0)-kronecker(rep(1,d),Z[1,]<0),nrow=n)
  return(X3)
}



###################################################################
#########################(Rand_CFUSN) ################
#######################################################################

#' Random multivariate CFUSN
#'
#' Generate random d-vectors  from the multivariate
#' Canonical Fundamental Skew-Normal (CFUSN) distribution
#' @param n The number of variates to be generated
#' @param Delta Correlation matrix, the skewness matrix Delta
#'
#' @examples
#' d <- 2; p <- 3
#' Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
#' ieg<- eigen(diag(p)+t(Lamd)%*%Lamd)
#' V <- ieg$vectors
#' Delta <-Lamd %*% V %*% diag(1/sqrt(ieg$values)) %*% t(V)
#' x<-Rand_CFUSN(20,Delta)
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 (5.5) p.247
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Random generation
#' @export

Rand_CFUSN <- function(n,Delta){
  #
  d <- dim(Delta)[1]
  p <- dim(Delta)[2]
  iegD<- eigen(diag(d)-Delta%*%t(Delta))
  VD <- iegD$vectors

  Delta2 <-VD %*% diag(sqrt(iegD$values)) %*% t(VD)
  Z1 <- abs(matrix( rnorm(p*n), nrow = n))
  Z2 <- abs(matrix( rnorm(d*n), nrow = n))

  X2 <- apply(Z2, 1, function(y) Delta2%*%(y))
  #Delta2%*%Z2[1,]
  X <- Delta%*%t(Z1)+X2
  return(t(X))
}


#####################################################################
#################### Canonical Fundamental Skew-Spherical Distribution CFUSSD
###############################################################


################# #################
#' Random multivariate CFUSSD
#'
#' Generate random d-vectors  from the multivariate Canonical Fundamental
#' Skew-Spherical distribution (CFUSSD)  with Gamma generator
#' @param n sample size
#' @param d  dimension
#' @param p  dimension of the first term
#' @param a  shape  parameter of the Gamma generator
#' @param b  scale parameter of the Gamma generator
#' @param Delta skewness matrix
#'
#' @return  a random matrix \eqn{n x d}
#'
#' @examples
#' d <- 2; p <- 3
#' Lamd <-  matrix(sample(1:50-25, d*p), nrow=d)
#' ieg<- eigen(diag(p)+t(Lamd)%*%Lamd)
#' V <- ieg$vectors
#' Delta <-Lamd %*% V %*% diag(1/sqrt(ieg$values)) %*% t(V)
#' Rand_CFUSSD_G(20,p,d,1,1,Delta)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, (5.36) p. 266, (see p.247 for Delta)
#' @family Random generation
#' @export
Rand_CFUSSD_G <- function(n,d,p,a,b,Delta){
  #
  ieg2 <- eigen(diag(d)-Delta%*%t(Delta))
  V2 <- ieg2$vectors
  Delta2 <-  V2 %*% diag(sqrt(ieg2$values)) %*% t(V2)
  R <- rgamma(n, shape= a, scale = b ) # a > 0 scale = 1/rate
  p1 <-  rbeta(n,p/2,d/2,n)
  # %% symmetric
  U2 <-  Rand_Uni_d(n,d)
  mU1 <- abs(Rand_Uni_d(n,p));
  Z1 <-  mU1*matrix(kronecker(rep(1,p),R*p1),ncol = p)
  Z2 <-  U2*matrix(kronecker(rep(1,d),R*(1-p1)),ncol = d)
  X <- t(Delta%*%t(Z1)+Delta2%*%t(Z2))
  return(X)
}


############################################
#########   ################
###########################################
#####
#####
#' MVSK Uniform on the sphere
#'
#'
#' Computes the theoretical values of the mean vector,  covariance, skewness vector,
#' total skenwness, kurtosis vector and total
#' kurtosis for the Uniform distribution on the sphere.  Note that Skewness is ZERO
#' @param d dimensions
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 Proposition 5.3 p.297
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Theoretical Moments, Cumulants
#' @export
MVSK_UniDistr_Th <- function(d,nCum = 0){
  eL <- c(0,2,0,0)
  loc_type_el <- Partition_Type_eL_Location(eL)
  L2 <- Commutator_Moment_L(eL,2,2,d)
  Idv <- as.vector(diag(d))
  EU1 <- rep(0,d)
  EU2 <- 1/d*Idv
  varU<- matrix(1/d*Idv,nrow = d)
  EU3 <- rep(0,d^3)
  EU4 <- 1/d/(d+2)*L2%*%kronecker(Idv,Idv)
  EU.k <- list(EU1,EU2,EU3,EU4)

  Umvsk <- list(EU1,varU,EU.k)
  names(Umvsk) <- c("EU1","varU","EU.k")
  if (nCum != 0){
    cumU.k <- Mom2CumMulti(EU.k)
    kurtU <- cumU.k[[4]]*d
    Umvsk <- list(EU1,varU,EU.k,cumU.k,kurtU)
    names(Umvsk) <- c("EU1","varU","EU.k","cumU.k","kurtU")
  }
  return(Umvsk)
}


#########################################################################
#' Moments and cumulants Uniform Distribution on the Sphere
#'
#' By default, only moments are provided
#'
#' @param r highest order of moments and cumulants
#' @param d dimension
#' @return The list of moments and cumulants in vector form
#' @examples
#' # The first four moments for d=3
#' MomCum_UniDistr_Th(4,3,nCum=0)
#' # The first four moments and cumulants for d=3
#' MomCum_UniDistr_Th(4,3,nCum=4)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 Proposition 5.3 p.297
#' @family Theoretical Moments, Cumulants
#' @export

MomCum_UniDistr_Th <- function(r,d,nCum=0)
{
  if (r%%2 !=0 )   return(EUM.k <- 0)# (4%%2)*2

  m0 <- floor(r/2)
  # 1d moments
  EU.k0 <- c(rep(0,r))
  EU.k0[seq(2,2*m0,2)] <- cumprod(seq(1,2*m0,2))/cumprod(seq(d,d+2*(m0-1),2))
  Idv <- as.vector(diag(d))
  per <- cumprod(seq(1,2*m0,2))
  Idv.u <- Idv
  EUM.k <- NULL
  EUM.k[[1]]  <- rep(0,d)
  EUM.k[[2]] <- 1/d*Idv
  EUM.k[[3]] <- c(rep(0,d^3))
  varU<- matrix(1/d*Idv,nrow = d)
  #,"varU","EU.k","cumU.k","kurtU")
  # names(EUM.k[1] ) <- c("EU.1")
  for (k in 2:m0) {
    # tic
    eL <- c(0,k,c(rep(0,2*k-2)))
    # loc_type_el <- Partition_Type_eL_Location(eL)
    # L2 <- Commutator_Moment_L(eL,loc_type_el[1],loc_type_el[2],d)
    L2 <- Commutator_Moment_eL(eL,d)
     Idv.u <- kronecker(Idv.u,Idv)
    EUM.k[[2*k]] <- EU.k0[2*k] *L2%*% Idv.u /per[k]
    EUM.k[[2*k+1]] <-c(rep(0,d^(2*k+1)))
    #   print( toc)
  }
  if ( r%%2 == 0 ) EUM.k <- EUM.k[-(2*m0+1)]
  if (nCum != 0) {
    CumU <- Mom2CumMulti(EUM.k )
    MomCumU <- list(EUM.k,CumU )
    names(MomCumU) <- c("EUM" , "CumU")
    return(MomCumU)
  }

  nev <- NULL
  for (k in 1:r) {
    nev[[k]] <- paste("EU.",as.character(k), sep='')
  }
  names(EUM.k) <- nev
  return(EUM.k)
}


#####################################################################
############################# UniAbsDistr_Th ##############
#############################################################
#' Moments of the modulus of the Uniform distribution on the sphere
#'
#' Moments (up to the 4th order) of the modulus of the d-variate Uniform distribution on
#' the sphere on (d-1)
#' @param d vector-dimension
#' @return  The list of the first four moments in vector form
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Lemma 5.12 p.298
#' @family Theoretical Moments, Cumulants
#' @export
MomCum_UniAbsDistr_Th <- function(d){
#
  Id <- diag(d)
EabsU1 <- matrix(sqrt(1/pi)/Gkd(1,d)*rep(1,d),ncol = 1)
EabsU2 <- matrix(as.vector(Id)/d+(rep(1,d^2)-as.vector(Id))/pi/Gkd(2,d),ncol = 1)
###########
e3 <- rep(0,d^3)
for (k in c(1:d)) {
  e3 <- e3 + kron3(Id[,k])
}
e3 <- matrix(e3,ncol = 1)
e21 <- rep(0,d^3)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e21 <- e21 +  kronecker(kron2(Id[,k]),Id[,j])
  }
}
e21 <- Matrix_Symmetry(d,3) %*%e21*3
e111 <- rep(0,d^3)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          e111 <- e111 +  kronecker(Id[,k],kronecker(Id[,j],Id[,m]))
        }
      }
    }
  }
}
e111<-matrix(e111,ncol=1) # 6*Matrix_Symmetry(d,3) %*% e111
EabsU3 <-( e3 + e21/2 + e111/pi)/sqrt(pi)/Gkd(3,d)
#########################################################################
e4 <- rep(0,d^4)
for (k in c(1:d)) {
  e4 <- e4 + kron4(Id[,k])
}
e4 <-matrix(e4,ncol = 1)
##############################################################x
e31 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e31 <- e31 +  kronecker(kron3(Id[,k]),Id[,j])
  }
}
e31 <- Matrix_Symmetry(d,4)%*%e31*4
###########################################
e211 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          e211 <- e211 +  kronecker(kron2(Id[,k]),kronecker(Id[,j],Id[,m]))
        }
      }
    }
  }
}
e211 <- Matrix_Symmetry(d,4)%*%e211*6
##################################
e22 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k)  e22 <- e22 + kronecker(kron2(Id[,k]),kron2(Id[,j]))
  }
}
e22 <- matrix(e22,ncol = 1)
###########################################
e1111 <- rep(0,d^4)
for (k in c(1:d)) {
  for (j in c(1:d) ) {
    if (j != k) {
      for (m in c(1:d)) {
        if (m !=k & m!= j){
          for (n in c(1:d)) {
            if (n !=k & n!= j & n!=m) {
              e1111 <- e1111 + kronecker(kronecker(Id[,k],Id[,j]),kronecker(Id[,m],Id[,n]))}
          }

        }
      }
    }
  }
}
e1111 <- Matrix_Symmetry(d,4)%*%e1111
###########################################################
e.vect4 <-( e31/pi +  e22/4 + e211/sqrt(pi^3)  +  e1111/pi^2)/Gkd(4,d)
EabsU4 <- 3*e4/d/(d+2) +  e.vect4
Eabs <- list(EabsU1,EabsU2,EabsU3,EabsU4 )
names(Eabs) <- c("EabsU1","EabsU2","EabsU3","EabsU4" )
return(Eabs)
}



 #################################
 ###############################################################################
 ##### Mod Z    ##########
 ###################################
 # n <- 9
#' Moments and cumulants Central folded Normal distribution
#'
#' Provides the theoretical moments and cumulants of the univariate Central Folded
#' Normal distribution. By default only moments are provided.
#' @param n The highest moment (cumulant) order
#' @return  The list of moments and cumulants
#' @examples
#' # The first three moments
#' MomCum_Zabs(3, nCum = F)
#' # The first three moments and cumulants
#' MomCum_Zabs(3, nCum = T)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Proposition 5.1 p.242 and  formula: p. 301
#' @export
 MomCum_Zabs <- function(n,nCum=FALSE){
   #
   n0 <- floor((n-2)/2) +(n/2-n%/%2>0)
   mu <- NULL
   mu[1] <- sqrt(2/pi)
   mu[2] <- 1
   double.fact1 <- NULL
   double.fact2 <- NULL
   double.fact1[1] <- 1
   double.fact2[1] <- 1
   #sz <- 1
   if (n0>0){
     for (sz in 1:n0) {
       double.fact1[sz+1] <- (2*sz)*double.fact1[sz]
       mu[2*sz+1] <- mu[1]*double.fact1[sz+1]
       double.fact2[sz+1] <- (2*sz+1)*double.fact2[sz]
       mu[2*sz+2] <- double.fact2[sz+1]
     }
   }
   if (n/2-n%/%2 >0) mu <- mu[-(2*n0+2)]
   if (nCum==1){
     cum <- Mom2Cum(mu )
     MomCumZ <- list(mu,cum )
     names(MomCumZ) <- c("MuZ" , "CumZ")
     return(MomCumZ)}
   return(mu)
 }

##########################

####################################################
#######
########################################################
###########       SkewNormMVSK_Th ###################
#######################################################
 #' MVSK multivariate Skew Normal
 #'
 #' Computes the theoretical values of the mean vector,  covariance, skewness vector,
 #' total skenwness, kurtosis vector and total
 #' kurtosis for the multivariate Skew Normal distribution
 #' @param Om A \eqn{d \times d} correlation matrix
 #' @param alf shape parameter d-vector
 #' @return  A list of theoretical values
 #'
 #' @examples
 #' alpha<-c(10,5,0)
 #' omega<-diag(rep(1,3))
 #' MVSK_SkewNorm_Th(omega,alpha)
 #' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
 #' Springer 2021 (5.5) p.247
 #' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
 #' skewness and kurtosis. Sankhya A, 83(2), 607-644.
 #' @family Theoretical Moments, Cumulants
 #' @export

MVSK_SkewNorm_Th <- function(Om,alf){
  # Lemma 5.1 p. 246
  # Om correlation matrix
  #if (nMu !=0) {nCum=nMu}
  # kron3 <- function(A) kronecker(A,kronecker(A,A))
  oszt <- sqrt(1+t(alf)%*%Om%*%alf)
  dim(oszt) <- NULL
  delt <- Om%*%alf/oszt
  EX <- sqrt(2/pi)*delt
  VarX <- Om-2/pi*delt%*%t(delt)
  ieg <- eigen(VarX)
  V <- ieg$vectors
  VarX.inv <-  V %*% diag(1/sqrt(ieg$values)) %*% t(V)
  VarXid <- VarX.inv %*% delt
  kr3 <- kron3(VarXid)
  SkewX <-(2*sqrt(2/pi)^3-sqrt(2/pi))*kr3
  SkewX.tot <- sum(SkewX*SkewX)
  kz4  <-  -6*4/(pi^2)+4*2/pi
  KurtX <- kz4*kronecker(kr3,VarXid)
  KurtX.tot <- sum(KurtX*KurtX)
  #if (nCum!=0){
  #  deltk <- NULL
  #  deltk[[1]] <- delt
  #  deltk[[2]] <- kronecker(delt,delt)
  #  MC.absZ <- MomCum_Zabs(r,nCum = 1)
  #  cum <- NULL
  #  cum[[1]] <- EX
  #  cum[[2]] <- as.vector(VarX)
  #  for (k in 3:r) {
  #    deltk[[k]] <- kronecker(delt,deltk[[k-1]])
  #    cum[[k]] <- MC.absZ$CumZ[k]* deltk[[k]]
  #  }
  #  if (nMu!=0){
  #    Mu <- Cum2MomMulti(cum)
  #    MVSKTheo <- list(EX,VarX, SkewX, SkewX.tot,
  #                     KurtX, KurtX.tot,Mu,cum )
  #    names(MVSKTheo) <- c("EX","VarX", "SkewX",
  #                         "SkewX.tot","KurtX", "KurtX.tot", "MuX","CumX")
  #    return(MVSKTheo)
  #  }
  #  MVSKTheo <- list(EX,VarX, SkewX, SkewX.tot,KurtX, KurtX.tot,cum )
  #  names(MVSKTheo) <- c("EX","VarX", "SkewX", "SkewX.tot","KurtX", "KurtX.tot","CumX")
  #  return(MVSKTheo)
  #
  # }
  MVSKTheo <- list(EX,VarX, SkewX, SkewX.tot,KurtX, KurtX.tot )
  names(MVSKTheo) <- c("EX","VarX", "SkewX", "SkewX.tot","KurtX", "KurtX.tot")
  return(MVSKTheo)
}

########################################################
###########       MomCum_SkewNorm_Th ###################
#######################################################
#' Moments and cumulants d-variate Skew Normal
#'
#' Computes the theoretical values of moments and cumulants up to the r-th order
#' @param r the highest moment and cumulant order
#' @param Om A \eqn{d \times d} correlation matrix
#' @param alf shape parameter d-vector
#' @return  A list of theoretical moments and cumulants
#'
#' @examples
#' alpha<-c(10,5,0)
#' omega<-diag(rep(1,3))
#' MomCum_SkewNorm_Th(r=4,omega,alpha)
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021 (5.5) p.247
#' @references S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
#' skewness and kurtosis. Sankhya A, 83(2), 607-644.
#' @family Theoretical Moments, Cumulants
#' @export

MomCum_SkewNorm_Th <- function(r=4,Om,alf){
  # Lemma 5.1 p. 246
  # Om correlation matrix
  oszt <- sqrt(1+t(alf)%*%Om%*%alf)
  dim(oszt) <- NULL
  delt <- Om%*%alf/oszt
  EX <- sqrt(2/pi)*delt
  VarX <- Om-2/pi*delt%*%t(delt)
  deltk <- NULL
  deltk[[1]] <- delt
  deltk[[2]] <- kronecker(delt,delt)
  MC.absZ <- MomCum_Zabs(r,nCum = 1)
  cum <- NULL
  cum[[1]] <- EX
  cum[[2]] <- as.vector(VarX)
  for (k in 3:r) {
    deltk[[k]] <- kronecker(delt,deltk[[k-1]])
    cum[[k]] <- MC.absZ$CumZ[k]* deltk[[k]]
  }

  Mu <- Cum2MomMulti(cum)
  MVSKTheo <- list(Mu,cum )
  names(MVSKTheo) <- c("MuX","CumX")
  return(MVSKTheo)

}



################################## Multi  ##############
#' Moments and cumulants multivariate central folded Normal distribution
#'
#' Provides the theoretical moments and cumulants of the multivariate central Folded
#' Normal distribution. By default only moments are provided.
#' @param r The highest moment (cumulant) order
#' @param d dimension
#' @return  The list of moments and cumulants in vector form
#'
#'
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Lemma 5.2 p. 249
#' @family Theoretical Moments, Cumulants
#' @export
MomCum_MZabs_Th <- function(r,d){
  #
  MC.Z <- MomCum_Zabs(r,nCum = 1)
  Id <- diag(d)
  i_kron_k <- NULL
  i_kron_k[[1]] <- rep(1,d)
  e_kron_k <- NULL
  e_kron_k[[1]] <-Id
  CumMZ <- NULL
  for (k in 1:r) {
    i_kron_k[[k+1]] <- rep(0,d^(k+1))
    #e_kron_k[[k+1]] <- NULL
    e_kron_k[[k+1]] <- matrix(rep(0,d^(k+2)),ncol=d)
    for (j in 1:d) {
      e_kron_k[[k+1]][,j] <- kronecker(e_kron_k[[k]][,j],Id[,j])
      i_kron_k[[k+1]] <- i_kron_k[[k+1]]+e_kron_k[[k+1]][,j]
    }
    CumMZ[[k]] <-  MC.Z$CumZ[k]*i_kron_k[[k]]
  }
  return(CumMZ)
}


#############################################
#' Moments and cumulants CFUSN
#'
#' Provides the theoretical cumulants of the multivariate Canonical Fundamental
#' Skew Normal distribution
#' @param r The highest cumulant order
#' @param d The multivariate dimension and number of rows of the skewness matrix Delta
#' @param p The number of cols of the skewness matrix Delta
#' @param Delta The skewness matrix
#' @return  The list of theoretical  cumulants in vector form
#' @references Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
#' Springer 2021, Lemma 5.3 p.251
#' @family Theoretical Moments, Cumulants
#' @export

MomCum_CFUSN_Th <- function(r,d,p,Delta){
  #
  EX <- sqrt(2/pi)*Delta%*%rep(1,p)
  VarX <- diag(d)-2/pi*Delta%*%t(Delta)
  MC.Z <- MomCum_Zabs(r,nCum = 1)
  Ip <- diag(p)
  i_kron_k <- NULL
  i_kron_k[[1]] <- rep(1,p)
  e_kron_k <- NULL
  e_kron_k[[1]] <- Ip
  De_kron_k <- NULL
  De_kron_k[[1]] <- Delta
  CumMX <- NULL
  CumMX[[1]] <- EX
  for (k in 1:r) {
    i_kron_k[[k+1]] <- rep(0,p^(k+1))
    # e_kron_k[[k+1]] <- NULL
    e_kron_k[[k+1]] <- matrix(rep(0,p^(k+2)),ncol=p)
    for (j in 1:p) {
      e_kron_k[[k+1]][,j] <- kronecker(e_kron_k[[k]][,j],Ip[,j])
      i_kron_k[[k+1]] <- i_kron_k[[k+1]]+e_kron_k[[k+1]][,j]
    }
    De_kron_k[[k+1]] <- kronecker(De_kron_k[[k]],Delta)
    CumMX[[k+1]] <-  MC.Z$CumZ[k]*
      De_kron_k[[k+1]]%*%(i_kron_k[[k+1]])
  }
  return(CumMX)
}

