
########### Support
##
## 1. KronProd
## 2. perm

## 3. KronPower
## 4. KronProdCells
## 5. KronProdVects ??????
## 6.
## 7. Scaling_MultiSample
## 8.
## 9.
## 10.

#####################
kron2 <- function(A) kronecker(A,A)
kron3 <- function(A) kronecker(A,kronecker(A,A))
kron4 <- function(A) kronecker(A,kron3(A))
# kron4
# Gkd
Gkd <- function(k,d) gamma((d+k)/2)/gamma(d/2)
##################################
#KronProdList
###############################
#
#Kronecker product of a list
# @examples
# Mc<-list(c(1,2),c(1,2,3,4),c(5,5,5,5))
# KronProd(Mc)

KronProd<-function(Mc){
  m<-length(Mc)
  Kr<-Mc[[1]]
  if (m==1) {return("Kr"=Kr)}
  for (k in 2:m) {Kr<-kronecker(Kr,Mc[[k]])}
  return("Kr"<-Kr)
}

#########################################
#Kronecker power of a vector
#################################################
# x <- c(1:4)
KronPower <- function(x,k){
  Xk <-  rep(list(x),k)
  return(KronProd(Xk))
}



# is.scalar
#
# Check if an element is scalar
# @examples
# a<-c(1,2)
# is.scalar(a)
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

## arrangements::permutations(1:q)
# instead of
# perm <- function(v) {
#   n <- length(v)
#   if (n == 1) v
#   else {
#     X <- NULL
#     for (i in 1:n) X <- rbind(X, cbind(v[i],perm(v[-i]) ))
#     X
#   }
# }

##############
### L22_H4 commutator matrix for Variance_of_Esti_Kurt
L22_H4 <- function(d){
  sz <- 1
  v1 <- 1:4
  v2 <- 5:8
  perm22_H4 <- matrix(rep(0,72*8),nrow = 72,ncol = 8)
  Cperm1  <-  arrangements::combinations(4, 2)#(v1,2);
  Cperm2  <-  arrangements::combinations(v2,2);
  for (k in 1:6){
    for (j in 1:6){
      A <- c(setdiff(1:4,Cperm1[j,]), setdiff(5:8,Cperm2[k,]))
      perm22_H4[sz,]  <-  c(A[1], A[3], A[2], A[4], Cperm1[j,], Cperm2[k,])
      perm22_H4[sz+1,]  <- c(A[1], A[4], A[2], A[3], Cperm1[j,], Cperm2[k,])
      sz <- sz+2;
    }
  }
  L22_H4  <- matrix(rep(0,d^16),nrow=d^8);
  for (k in 1:72)
    L22_H4  <- L22_H4 + Commutator_Kperm(perm22_H4[k,],d )
  end
  return(L22_H4)
}

#
# #' Third d-variate T-Hermite polynomial at standardized vector x
# #' @param x multivariate data of size d
Hermite_Third<-function(x){
  d=length(x)
  H3<-Hermite_Poly_HN_Multi(x,diag(d),3)[[3]]
  return(as.vector(H3))
}

# #' Fourth d-variate T-Hermite polynomial at standardized vector x
# #' @param x multivariate data of size d
Hermite_Fourth<-function(x){
  d=length(x)
  H4<-Hermite_Poly_HN_Multi(x,diag(d),4)[[4]]
  return(as.vector(H4))
}
