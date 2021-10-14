#KronProdList
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

# is.scalar
#
# Check if an element is scalar
# @examples
# a<-c(1,2)
# is.scalar(a)
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

##
perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i],perm(v[-i]) ))
    X
  }
}
