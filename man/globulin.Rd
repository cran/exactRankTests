\name{globulin}
\alias{globulin}
\non_function{}
\title{ Differences in Globulin Fraction in Two Groups }
\usage{data(globulin)}
\description{
Globulin fraction of plasma (g/l) in two groups of 10 patients. 
}
\format{
  This data frame contains the following variables:
  \describe{
    \item{group1}{Globulin fraction of plasma for group 1}
    \item{group2}{Globulin fraction of plasma for group 2}
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
  Joachim R\"ohmel (1996).
  Precision intervals for estimates of the difference in success 
  rates for binary  random variables based on the Permutation Principle.
  \emph{Biometrical Journal}, \bold{38}(8), 977--993.  
}
\examples{
data(globulin)
attach(globulin)
pt <- perm.test(group1, group2, conf.int=TRUE)
pt
stopifnot(pt$conf.int == c(-8.50, 1.25))
}
\keyword{datasets}
