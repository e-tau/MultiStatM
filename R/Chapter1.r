
## 2. Permutation_Inverse
## 4. Matrix_Symmetry - calculates symmetrizer; 1.3.1 Symmetrization, p.14. (1.29)
## 5. Matrix_Elimination 1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
## 6. Matrix_Qplication, p.21, (1.31)

############################
### 4

#' Inverse of a Permutation
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, Remark 1.1, p.2
#'
#' @param permutation0 A permutation of numbers 1:n
#' @return Inverse permutation of  permutation0
#'
#' @export
Permutation_Inverse  <-  function(permutation0)
{ return(order(permutation0))}


#' Symmetrizer Matrix
#'
#' It does not use sparse matrices or parallel computations
#'
#'  @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.Section 1.3.1 Symmetrization, p.14. (1.29)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#' @return Symmetrizer Matrix with order d^q x d^q
#'
#' @family Matrix
#' @export
Matrix_Symmetry  <-  function(d,q){

  Symmetr <- array(0,c(d^q,d^q))
  Pall  <-  arrangements::permutations(1:q)
  db  <-  factorial(q)
  for(k in 1:db){
    Symmetr  <-  Symmetr + Commutator_Kperm(Pall[k,],d )
  }
  Symmetr  <-  Symmetr/factorial(q)
  return(Symmetr)
}


##########################
### 5

#' Elimination Matrix
#'
#' Eliminated the duplicated/q-plicated elements in the vector
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.
#' Section  1.3.2 Multi-Indexing, Elimination, and Duplication, p.21,(1.32)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#' @return Elimination Matrix of order  ŋ\eqn{_{d,q}x d^q}, see (1.30), p.15 # \eqn{_{d,q}}
#'
#' @family Matrix
#' @export
Matrix_Elimination<-function(d,q){
  x<- primes::generate_n_primes(1005)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  int<-match(unique(y), y)
  I<-diag(d^q)
  return(I[int,])
}

###############################
### 6

#' Qplication Matrix
#'
#' Restores the duplicated/q-plicated  elements which are eliminated
#' by  Matrix_Elimination in the vector,
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, p.21, (1.31)
#'
#' @param d dimension of a vector x
#' @param q  power of the Kronecker product
#' @return Qplication Matrix of order \eqn{d^q x ŋ_{d,q}}, see (1.30), p.15
#'
#' @family Matrix
#' @export
Matrix_Qplication<-function(d,q){
  x<- primes::generate_n_primes(1005)
  x<- x[1:d]
  y <- x
  for (k in 1:(q-1)){ y <- kronecker(x,y)}
  int<-match(unique(y), y)
  im<-match(y, unique(y))
  esz <- length(int)
  De <- matrix(0,d^2,esz)
  I<-diag(d^q)

  for (k in 1:esz){
    if (sum((im==k))>1){De[,k]<-apply(I[ ,(im==k)],1,sum)}
    else De[,k]<-I[,(im==k)]
  }
  return(De)
}
