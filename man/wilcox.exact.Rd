\name{wilcox.exact}
\alias{wilcox.exact}
\alias{wilcox.exact.default}
\alias{wilcox.exact.formula}
\title{Wilcoxon Rank Sum and Signed Rank Tests}
\description{
  Performs one and two sample Wilcoxon tests on vectors of data for possibly
tied observations.
}
\usage{
\method{wilcox.exact}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"),
             mu = 0, paired = FALSE, exact = NULL,  
             conf.int = FALSE, conf.level = 0.95, \dots)
\method{wilcox.exact}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{numeric vector of data values.}
  \item{y}{an optional numeric vector of data values.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{mu}{a number specifying an optional location parameter.}
  \item{paired}{a logical indicating whether you want a paired test.}
  \item{exact}{a logical indicating whether an exact p-value should be
    computed.}
  \item{conf.int}{a logical indicating whether a confidence interval
    should be computed.}
  \item{conf.level}{confidence level of the interval.}
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
  This version computes exact conditional (on the data) 
  p-values and quantiles using the Shift-Algorithm by Streitberg & 
  R\"ohmel for both tied and untied samples.  

  If only \code{x} is given, or if both \code{x} and \code{y} are given
  and \code{paired} is \code{TRUE}, a Wilcoxon signed rank test of the
  null that the median of \code{x} (in the one sample case) or of
  \code{x-y} (in the paired two sample case) equals \code{mu} is
  performed.

  Otherwise, if both \code{x} and \code{y} are given and \code{paired}
  is \code{FALSE}, a Wilcoxon rank sum test (equivalent to the
  Mann-Whitney test) is carried out.  In this case, the null hypothesis
  is that the location of the distributions of \code{x} and \code{y}
  differ by \code{mu}.

  By default (if \code{exact} is not specified), an exact p-value is
  computed if the samples contain less than 50 finite values and there
  are no ties.  Otherwise, a normal approximation is used.

  Optionally (if argument \code{conf.int} is true), a nonparametric
  confidence interval for the median (one-sample case) or for the
  difference of the location parameters \code{x-y} is computed.  If
  exact p-values are available, an exact confidence interval is obtained
  by the algorithm described in Bauer (1972).  Otherwise, an asymptotic
  confidence interval is returned.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic with a name
    describing it.}
  \item{p.value}{the p-value for the test.}
  \item{pointprob}{this gives the probability of observing the test
   statistic itself (called \code{point-prob}).}
  \item{null.value}{the location parameter \code{mu}.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{the type of test applied.}
  \item{data.name}{a character string giving the names of the data.}
  \item{conf.int}{a confidence interval for the location parameter.
    (Only present if argument \code{conf.int = TRUE}.)}
  \item{estimate}{Hodges-Lehmann estimate of the location parameter.
    (Only present if argument \code{conf.int = TRUE}.)}
}
\references{
  Myles Hollander & Douglas A. Wolfe (1973),
  \emph{Nonparametric statistical inference}.
  New York: John Wiley & Sons.
  Pages 27--33 (one-sample), 68--75 (two-sample).

  David F. Bauer (1972),
  Constructing confidence sets using rank statistics.
  \emph{Journal of the American Statistical Association}
  \bold{67}, 687--690.

  Cyrus R. Mehta & Nitin R. Patel (2001),
  \emph{StatXact-5 for Windows.}
  Manual, Cytel Software Cooperation, Cambridge, USA

}
\seealso{
  \code{\link{perm.test}} for the one and two sample permutation test.
}
\examples{
## One-sample test.
## Hollander & Wolfe (1973), 29f.
## Hamilton depression scale factor measurements in 9 patients with
##  mixed anxiety and depression, taken at the first (x) and second
##  (y) visit after initiation of a therapy (administration of a
##  tranquilizer).
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
wilcox.exact(x, y, paired = TRUE, alternative = "greater")
wilcox.exact(y - x, alternative = "less")    # The same.

## Two-sample test.
## Hollander & Wolfe (1973), 69f.
## Permeability constants of the human chorioamnion (a placental
##  membrane) at term (x) and between 12 to 26 weeks gestational
##  age (y).  The alternative of interest is greater permeability
##  of the human chorioamnion for the term pregnancy.
x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
wilcox.exact(x, y, alternative = "g")        # greater

## Formula interface.
data(airquality)
boxplot(Ozone ~ Month, data = airquality)
wilcox.exact(Ozone ~ Month, data = airquality,
            subset = Month \%in\% c(5, 8))


# Hollander & Wolfe, p. 39, results p. 40 and p. 53

x <- c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

wilcox.exact(y,x, paired=TRUE, conf.int=TRUE)

# Hollander & Wolfe, p. 110, results p. 111 and p. 126

x <- c(0.8, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)

wilcox.exact(y,x, conf.int=TRUE)

}
\keyword{htest}
