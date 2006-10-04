\name{globulin}
\alias{globulin}
\docType{data}
\title{ Differences in Globulin Fraction in Two Groups }
\usage{data(globulin)}
\description{
Globulin fraction of plasma (g/l) in two groups of 10 patients. 
}
\format{
  This data frame contains the following variables:
  \describe{
    \item{gfrac}{Globulin fraction of plasma}
    \item{group}{a factor with levels \code{group1} and \code{group2}}
  }
}
\details{
See page 75 of Gardner & Altman (1989).
}
\source{
  M. J. Gardner & D. G. Altman (1989), 
  \emph{Statistics with Confidence}.
  Published by the British Medical Journal.
}
\references{
  Joachim R\"ohmel (1996),
  Precision intervals for estimates of the difference in success 
  rates for binary  random variables based on the permutation principle.
  \emph{Biometrical Journal}, \bold{38}(8), 977--993.  
}
\examples{
data(globulin)
perm.test(gfrac ~ group, data=globulin, conf.int=TRUE)
}
\keyword{datasets}
