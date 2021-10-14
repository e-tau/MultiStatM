###########
## 1. KommMatrix  1.2.3 Commutation Matrix p. 7, (1.12)
## 2. pik - used for PermutingMatrix
## 3. PermutatingMatrix - used for symmetrizer
## 4. Symmetrizing - calculates symmetrizer
## 5. Elimination matrix
## 6. Qplication


### 1.2.3 Commutation Matrix p. 7, (1.12)

#' Commutator_Kmn
#'
#' Provides the commutation matrix
#' @param m Row-dimension
#' @param n Col-dimension
#' @return An commutation matrix matrix of dimension \eqn{mn x mn}
#' @references G.H.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021 (p.7, (1.12)).
#' @export
Commutator_Kmn <-  function(m,n){
  pp <- 1:(m*n)
  p2 <- t(matrix(data=pp,nrow = m, ncol = n))
  pc <- as.vector(p2)
  Ic <- diag(m*n)
  Ic <- Ic[pc,]
  return(Ic)
}

#######################
### 2

pik = function(perm,d,k){
  out = 1 + (perm-1)%*%cumprod(d*array(1,c(k,1)))/d
  return(out)
}

###################
### 3

# Permutation matrix of any order
#
# It does not use sparse matrices
#
#
#PermutingMatrix = function(d,permutation0){
#  k = length(permutation0)
#  egyk = array(1,c(k,1))
#  perm0 = permutation0
#  perm01 = rev(perm0)
#  # sr = sort(perm01)
#  perm02 = order(perm01)
#  permutation0 = rev(perm02)
#  # sr = sort(permutation0)
#  invIk0 = order(permutation0)
#
#  Allind = unique(t(combn(rep(1:d, k), k)))
#  K_perm = array(0, c(d^k,d^k))
#  for(jj in 1:(d^k)){
#    pA = Allind[jj,]
#    pInv = pA[invIk0]
#    K_perm[pik(pInv,d,k), pik(pA,d,k)] = 1
#  }
#
#  return(K_perm)
#}


############################
### 4

#' Matrix_Symmetry
#'
#' It does not use sparse matrices or parallel computations
#'
#' @export
Matrix_Symmetry = function(d,n){

  Symmetr=array(0,c(d^n,d^n))
  Pall = arrangements::permutations(1:n)
  db = factorial(n)
  for(k in 1:db){
    Symmetr = Symmetr + PermutingMatrix(d, Pall[k,])
  }
  Symmetr = Symmetr/factorial(n)
  return(Symmetr)
}


##########################
### 5

#' Matrix_Elimination
#'
#' Eliminated the duplicated elements in the vector
#'
#' @export
Matrix_Elimination<-function(d,q){
  x<- primes::generate_n_primes(1005)
  x<- x[1:d]
  y=x
  for (k in 1:(q-1)){ y=kronecker(x,y)}
  int<-match(unique(y), y)
  I<-diag(d^q)
  return(I[int,])
}

###############################
### 6

#' Matrix_Qplication
#'
#' Restores the duplicated elements in the vector
#'
#' @export
Matrix_Qplication<-function(d,q){
  x<- primes::generate_n_primes(1005)
  x<- x[1:d]
  y=x
  for (k in 1:(q-1)){ y=kronecker(x,y)}
  int<-match(unique(y), y)
  im<-match(y, unique(y))
  esz=length(int)
  De=matrix(0,d^2,esz)
  I<-diag(d^q)

  for (k in 1:esz){
    if (sum((im==k))>1){De[,k]<-apply(I[ ,(im==k)],1,sum)}
    else De[,k]<-I[,(im==k)]
  }
  return(De)
}
