\name{lungcancer}
\alias{lungcancer}
\docType{data}
\title{ Lung Cancer Clinical Trial }
\description{
  Survival times for patients suffering lung cancer for a treatment and
control group.
}
\usage{data(lungcancer)}
\format{
  A data frame with 14 observations on the following 3 variables.
  \describe{
    \item{time}{survival time in days.}
    \item{cens}{censoring indicator: 0 censored, 1 event.}
    \item{group}{a factor with levels \code{control} and \code{newdrug}.}
  }
}
\details{
  The data is given in Table 9.19, page 293, of Metha and Pathel (2001). The
two-sided p-value for the log-rank test is 0.001 (page 295).

}
\references{

        Cyrus R. Mehta & Nitin R. Patel (2001),
        \emph{StatXact-5 for Windows.}
        Manual, Cytel Software Cooperation, Cambridge, USA

}
\examples{
data(lungcancer)
attach(lungcancer)

# round logrank scores
scores <- cscores.Surv(cbind(time, cens))
T <- sum(scores[group=="newdrug"])
mobs <- sum(group=="newdrug")
(prob <- pperm(T, scores, m=mobs, al="le"))
pperm(T, scores, m=mobs, al="tw")
pperm(T, scores, m=mobs, al="tw", simulate=TRUE)

# map into integers, faster
scores <- cscores.Surv(cbind(time, cens), int=TRUE)
T <- sum(scores[group=="newdrug"])
mobs <- sum(group=="newdrug")
(prob <- pperm(T, scores, m=mobs, al="le"))
pperm(T, scores, m=mobs, al="tw")
pperm(T, scores, m=mobs, al="tw", simulate=TRUE)

detach(lungcancer)

}
\keyword{datasets}
