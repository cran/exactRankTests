\name{sal}
\alias{sal}
\docType{data}
\title{ Serum Antigen Level }
\description{
  The response of serum antigen level to AZT in 20 patients suffering AIDS.
}
\usage{data(sal)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{pre}{level pre treatment.}
    \item{post}{level post treatment.}
  }
}
\details{
  The data is given in Metha and Patel (2001), Table 8.1, page 181.
Two-sided p-value for the Wilcoxon-Signed Rank Test: 0.0021 (page 183) or
0.0038 (asymptotically).
}
\references{

        Cyrus R. Mehta & Nitin R. Patel (2001),
        \emph{StatXact-5 for Windows.}
        Manual, Cytel Software Cooperation, Cambridge, USA

}
\examples{
data(sal)
attach(sal)

wilcox.exact(pre, post, paired=TRUE, conf.int=TRUE)
wilcox.exact(pre,post, paired=TRUE, conf.int=TRUE, exact=FALSE)

detach(sal)

}
\keyword{datasets}
