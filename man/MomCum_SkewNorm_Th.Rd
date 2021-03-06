% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distributions.R
\name{MomCum_SkewNorm_Th}
\alias{MomCum_SkewNorm_Th}
\title{Moments and cumulants d-variate Skew Normal}
\usage{
MomCum_SkewNorm_Th(r = 4, Om, alf)
}
\arguments{
\item{r}{the highest moment and cumulant order}

\item{Om}{A \eqn{d \times d} correlation matrix}

\item{alf}{shape parameter d-vector}
}
\value{
A list of theoretical moments and cumulants
}
\description{
Computes the theoretical values of moments and cumulants up to the r-th order
}
\examples{
alpha<-c(10,5,0)
omega<-diag(rep(1,3))
MomCum_SkewNorm_Th(r=4,omega,alpha)
}
\references{
Gy.Terdik, Multivariate statistical methods - Going beyond the linear,
Springer 2021 (5.5) p.247

S. R. Jammalamadaka, E. Taufer, Gy. Terdik. On multivariate
skewness and kurtosis. Sankhya A, 83(2), 607-644.
}
\seealso{
Other Theoretical Moments, Cumulants: 
\code{\link{MVSK_SkewNorm_Th}()},
\code{\link{MVSK_UniDistr_Th}()},
\code{\link{MomCum_CFUSN_Th}()},
\code{\link{MomCum_MZabs_Th}()},
\code{\link{MomCum_UniAbsDistr_Th}()},
\code{\link{MomCum_UniDistr_Th}()}
}
\concept{Theoretical Moments, Cumulants}
