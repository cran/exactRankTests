\name{ocarcinoma}
\alias{ocarcinoma}
\docType{data}
\title{ Ovarian Carcinoma }
\description{
  Survival times of 35 women suffering ovarian carcinoma at stadium II and
IIA.
}
\usage{data(ocarcinoma)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{time}{time in days.}
    \item{cens}{censoring indicator: 0 censored, 1 event.}
    \item{stadium}{a factor at levels \code{II} and \code{IIA}. }
  }
}
\details{
  Data from Fleming et al. (1980, 1984), reanalysed in Schumacher and
Schulgen (2002).
}
\source{
  T. R. Fleming, J. R. O'Fallon, P. C. O'Brien and D. P. Harrington (1980).
Modified Kolmogorov-Smirnov test procedures with applications to arbitrarily
censored data. \emph{Biometrics}, \bold{36}, 607--625.

  T. R. Fleming, S. J. Green and D. P. Harrington (1984).
Considerations of monitoring and evaluating treatment effects in clinical
trials. \emph{Controlled Clinical Trials}, \bold{5}, 55--66.

}
\references{
  M. Schumacher and G. Schulgen (2002). \emph{Methodik klinischer Studien:
methodische Grundlagen der Planung, Durchführung und Auswertung.} Springer,
Heidelberg.
}
\examples{
data(ocarcinoma)
attach(ocarcinoma)
# compute integer valued logrank scores
logrsc <- cscores.Surv(cbind(time, cens), int=TRUE)
# the test statistic
lgT <- sum(logrsc[stadium == "II"])
# p-value
round(pperm(lgT, logrsc, m=sum(stadium=="II"), al="tw"), 4)
}
\keyword{datasets}
