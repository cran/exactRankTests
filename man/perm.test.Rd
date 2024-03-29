\name{perm.test}
\alias{perm.test}
\alias{perm.test.default}
\alias{perm.test.formula}
\title{One and Two Sample Permutation Test }
\description{
  Performs the permutation test for the one and two sample problem.
}
\usage{
\method{perm.test}{default}(x, y, paired=FALSE, alternative=c("two.sided", "less", "greater"),
          mu=0, exact=NULL, conf.int=FALSE, conf.level=0.95, tol=NULL, \dots)
\method{perm.test}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{numeric vector of integer data values.}
  \item{y}{numeric vector of integer data values.}
  \item{paired}{a logical indicating whether you want a paired test.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{mu}{a number specifying an optional location parameter.}
  \item{exact}{a logical indicating whether an exact p-value should be
    computed.}
  \item{conf.int}{a logical indicating whether a confidence interval
    should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{tol}{real. real valued scores are mapped into integers by
             multiplication. Make sure that the absolute difference between
             the "true" quantile and the approximated quantile is less
             than \code{tol}. This might not be possible due to memory/time limitations.
             See \code{\link{pperm}}.}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric variable giving the data values and \code{rhs} a factor
    with two levels giving the corresponding groups.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{further arguments to be passed to or from methods.}

}
\details{
  The permutation test is performed for integer valued observations or
  scores. If real values \code{x} or \code{y} are passed to this function
  the following applies: if \code{exact} is true (i.e. the sample size is
  less than 50 observations) and \code{tol} is not given, the scores are
  mapped into \eqn{\{1,\dots,N\}}, see \code{\link{pperm}} for the details.
  Otherwise the p-values are computed using \code{tol}. If the sample size
  exceeds $50$ observations, the usual normal approximation is used.

  P-values are computed according to the StatXact-manual, see
  \code{\link{pperm}}. 

  For (in principle) continuous variables the confidence sets represent the
  "largest shift in location being consistent with the observations". For
  discrete variables with only a few categories they are hard to interpret.
  In the case of binary data (e.g. success / failure) the confidence sets
  can be interpreted as the differences of two success-rates covered by the
  data. For a detailed description see R\"ohmel (1996). 

  Confidence intervals are only available for independent samples. When the
  sample sizes are unbalanced, \code{length(x)} needs to be smaller than
  \code{length(y)}.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic with a name
    describing it.}
  \item{p.value}{the p-value for the test.}
  \item{pointprob}{this gives the probability of observing the test
   statistic itself.}
  \item{null.value}{the location parameter \code{mu}.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{the type of test applied.}
  \item{data.name}{a character string giving the names of the data.}
  \item{conf.int}{a confidence interval for the location parameter.
    (Only present if argument \code{conf.int = TRUE}.)}
}
\references{ 

  Joachim R\"ohmel (1996),
  Precision intervals for estimates of the difference in success
  rates for binary  random variables based on the permutation principle.
  \emph{Biometrical Journal}, \bold{38}(8), 977--993.

  Cyrus R. Mehta & Nitin R. Patel (2001),
  \emph{StatXact-5 for Windows.}
  Manual, Cytel Software Cooperation, Cambridge, USA

}

\note{ Confidence intervals may need some cpu-time ... }

\examples{

# Example from Gardner & Altman (1989), p. 30
# two treatments A and B, 1 means improvement, 0 means no improvement
# confidence sets cf. R\"ohmel (1996)

A <- c(rep(1, 61), rep(0, 19))
B <- c(rep(1, 45), rep(0, 35))
perm.test(A, B, conf.int=TRUE, exact=TRUE)

# one-sample AIDS data (differences only), Methta and Patel (2001),
# Table 8.1 page 181

data(sal)
attach(sal)
ppdiff <- pre - post
detach(sal)

# p-values in StatXact == 0.0011 one-sided, 0.0021 two.sided, page 183

perm.test(ppdiff)
perm.test(ppdiff, alternative="less")
perm.test(ppdiff, exact=FALSE)

}
\keyword{htest}
