% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation.r
\name{Skew_Mardia}
\alias{Skew_Mardia}
\title{Mardia Skewness index}
\usage{
Skew_Mardia(H3)
}
\arguments{
\item{H3}{the skewness vector.}
}
\value{
Mardia's skewness index
}
\description{
Computes Mardia's skewness index given the skewness vector.
It can be used to compute estimated or population values.
}
\examples{
#Compute Mardia's skewness index for the multivariate Skew-Normal
alpha<-c(10,5)
omega<-diag(rep(1,2))
H3SN<-MVSK_SkewNorm_Th(omega, alpha)$SkewX
Skew_Mardia(H3SN)


}
\references{
S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
skewness and kurtosis. Sankhya A, 83(2), 607-644, Example 1.
}
\seealso{
Other Indexes: 
\code{\link{Esti_Kurt_Mardia}()},
\code{\link{Esti_Kurt_Mori}()},
\code{\link{Esti_SkewInd_Mori}()},
\code{\link{Esti_SkewVec_Mori}()},
\code{\link{Esti_Skew_Mardia}()},
\code{\link{Kurt_Mardia}()},
\code{\link{Kurt_Mori}()},
\code{\link{Skew_Mori}()}
}
\concept{Indexes}
