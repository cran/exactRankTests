\name{cscores}
\alias{cscores}
\alias{cscores.default}
\alias{cscores.Surv}
\alias{cscores.factor}
\title{ Computation of Scores }
\description{
  This function can be used to compute several scores for a data vector.
}
\usage{
\method{cscores}{default}(y, type=c("Data", "Wilcoxon", "NormalQuantile", 
        "AnsariBradley", "Median", "Savage", "ConSal"), int=FALSE,
        maxs=length(y), \dots )
\method{cscores}{factor}(y, \dots)
\method{cscores}{Surv}(y, type="LogRank", int=FALSE, maxs=nrow(y), \dots)
}
\arguments{
  \item{y}{ a numeric, factor or logical vector or an object of class
           \code{Surv}. }
  \item{type}{a character string which specifies the type of the scores to be
              computed. \code{Data} just returns \code{y} if 
              \code{y} is numeric.}
  \item{int}{a logical, forcing integer valued scores. }
  \item{maxs}{an integer defining the maximal value of the scores
              if \code{int=TRUE}.}
  \item{\dots}{ additional arguments, not passed to anything at the moment. }
}
\details{
  This function will serve as the basis for a 
more general framework of rank and permutation tests in future versions of
this package. Currently, it is only used in the examples.

The logrank scores are computed as given in Hothorn & Lausen (2002).

If integer valued scores are requested (\code{int = TRUE}), the
\code{scores} are mapped into integers by
\code{round(scores*length(scores)/max(scores))}. See \code{\link{dperm}} for
more details.

\code{type} is self descriptive, except for \code{ConSal} which implements
scores suggested by Conover & Salsburg (1988).
}
\value{
  A vector of scores for \code{y} with an attribute \code{scores} indicating
the kind of scores used is returned.
}
\references{

        Torsten Hothorn & Berthold Lausen (2003),
        On the exact distribution of maximally selected rank statistics.
        \emph{Computational Statistics & Data Analysis}, \bold{43}(2),
        121-137.

	William J. Conover & David S. Salsburg (1988), Locally most powerful tests
	for detecting treatment effects when only a subset of patients can
	be expected to "respond" to treatment. \emph{Biometrics}, \bold{44},
	189-196.

}
\examples{

y <- rnorm(50)
# v.d. Waerden scores
nq <- cscores(y, type="Normal", int=TRUE)
# quantile for m=20 observations in the first group
qperm(0.1, nq, 20)

}

\keyword{misc}
