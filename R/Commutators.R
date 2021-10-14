#' Commutator_Kperm
#'
#' Produces any permutation of kronecker products of vectors of any length
#' @param perm vector indicating the permutation of the order in the kronecker product
#' @param dims vector indicating the dimensions of the vectors
#' @return A permutation matrix
#' @references G.H.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.
#' @examples
#' dims=c(2,3,2)
#' perm = c(1,3,2)
#' K.matr(perm,dims)
#' @examples
#' perm =c(3,1,4,2) dims=4 # All vectors with dimension 4
#' @export
Commutator_Kperm <- function(perm,dims){

  n <- length(perm)
  if (length(dims)== 1) {dims <- rep(dims,n)}
  S.perm <-        sort(perm , index.return = TRUE)
  perm <- S.perm$ix # % inverse ( sort is increasing)

  K.matr <- diag(prod(dims))
  if (length(perm)== 1) {return( K.matr )}
  U.before <- NULL
  U.after <- NULL
  for (i in 1 : (n - 1)) {
    # run loop (n-i) times
    for (j in 1 : (n - i)) {
      # compare elements
      if (perm[j] > perm[j + 1]) {
        if ((j-1)>0) {U.before <- diag(prod(dims[1:(j-1)]))}
        else {U.before <-1}

        if ((j+1)==n){ U.after <- 1 }
        else { U.after <- diag(prod(dims[(j+2):n]))}

        K1.matr <-  kronecker(U.before,kronecker(Commutator_Kmn(dims[j+1],dims[j]),U.after))
        # K2.matr <- (K1.matr,)
        K.matr <- K1.matr %*% K.matr
        temp.d <- dims[j]
        dims[j] <- dims[j+1]
        dims[j+1] <- temp.d
        temp.p <- perm[j]
        perm[j] <- perm[j + 1]
        perm[j + 1] <- temp.p
      }
    }
  }
  return(K.matr)
}

#' Commutator_Moment_L
#'
#'
#' @examples
#' n=4
#' r=2
#' m=1
#' d=2
#' PTA<-Partition_Type_All(n)
#' el_rm<-PTA$eL_r[[2]][1,]
#' el_r_m is always a vector

#n=3;r=1;m=1;d=2;
#PTA<-PartitionsType_All(n)
#el_rm<-PTA$eL_r
#el_rm
#el<-el_rm[[2]]
#MC<-MomentCommutators.L(el,m,r,d)
#MC


Commutator_Moment_L<-function(el_rm,m,r,d) {
  N<-length(el_rm)
  PTB<-Partition_Type_All(N)
  part_class<-PTB$Part.class
  S_N_r<-PTB$S_N_r
  S_m_j<-PTB$S_r_j

  sepL<-cumsum(S_N_r)
  sepS_r<-cumsum(S_m_j[[r]])

  if (r==1) {perm_Urk1<- 1:N
  L_eL<-rep(0,d^N)
  L_eL<-L_eL+Commutator_Kperm(perm_Urk1,d)
  return("LeL"=L_eL)
  }
  else {
    if (m==1) {l_ind<-sepL[r-1]+1} else {l_ind<-sepL[r-1]+sepS_r[m-1]+1}
    if (is.scalar(S_m_j[r])) {u_ind<-l_ind+S_m_j[[r]]-1}
    else {u_ind<-l_ind+S_m_j[[r]][m]-1}
    perm_Urk1<-matrix(0,S_m_j[[r]][m],N)
    sz<- 1
    for (k in l_ind:u_ind){
      perm_Urk1[sz,]<-Partition_2Perm(part_class[[k]])
      sz<-sz+1
    }
  }
  L_eL<-rep(0,d^N)
  for (ss in 1:dim(perm_Urk1)[1]) {
    L_eL<-L_eL+Commutator_Kperm(perm_Urk1[ss,],d)

  }
  return("L_eL"=L_eL)
}

#' Commutator_Mixing
#'
#' mixing commutator
#'
#' @param d1 dimension of the first group of vectors
#' @param d2 dimension of the second group of vectors
#' @references G.H.Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021. (4.58) p. 171
#'
#' @examples
#' d1 <- c(2, 3, 2)
#' d2<-  c(3 ,2, 2)
#' Commutator_Mixing(d1,d2)
#' @export
Commutator_Mixing <- function( d1,d2) {
# permutations
# dim(d1)<-dim(d2)
n <- length(d2)
i1<- matrix(data =c(1:n,(1:n)+n) , nrow = 2, ncol = n, byrow =TRUE)
fact_n<-factorial(n)  # number of permutations
B <- matrix(data = rep(c(1:n),fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
Permut <- perm(c(1:n))
q<-cbind(B, (Permut+n))
indUj<- c(i1) # reorder q
q<-q[,indUj] # permutations
# dimensions
d1B <- matrix(rep(d1,fact_n), nrow = fact_n, ncol = n, byrow =TRUE)
d2q<- matrix(rep(0,fact_n*n),nrow = fact_n, ncol = n) # zeros(fact_n,n);
for (k in c(1:fact_n)) {
  d2q[k,] <- d2[Permut[k,]]
}
Bd <- cbind(d1B, d2q)
Bdq <- Bd[,indUj] # dimensions with respect to permutations
# Commutator
M_m_n<-matrix(rep(0,(prod(d1)*prod(d2))^2),nrow = prod(d1)*prod(d2))
for  (kk in c(1:fact_n)) {
  M_m_n<- M_m_n + Commutator_Kperm(q[kk,],Bdq[kk,])
}
return(M_m_n)
}
