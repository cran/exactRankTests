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
    \item{stadium}{a factor with levels \code{II} and \code{IIA}. }
  }
}
\details{
  Data from Fleming et al. (1980, 1984), reanalysed in Schumacher and
Schulgen (2002).
}
\source{
  
    Thomas R. Fleming, Judith R. O'Fallon, Peter C. O'Brien & David 
    P. Harrington (1980),
    Modified Kolmogorov-Smirnov test procedures with applications to 
    arbitrarily censored data. \emph{Biometrics}, \bold{36}, 607--625.

    Thomas R. Fleming, Stephanie J. Green & David P. Harrington (1984),
    Considerations of monitoring and evaluating treatment effects in clinical
    trials. \emph{Controlled Clinical Trials}, \bold{5}, 55--66.

}
\references{

    Martin Schumacher & Gabi Schulgen (2002), \emph{Methodik klinischer 
    Studien: methodische Grundlagen der Planung, Durchf\"uhrung und
    Auswertung.} Springer, Heidelberg.

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

# compute logrank scores and simulate p-value
logrsc <- cscores.Surv(cbind(time, cens), int=FALSE)
# the test statistic
lgT <- sum(logrsc[stadium == "II"])
# p-value
round(pperm(lgT, logrsc, m=sum(stadium=="II"), al="tw", simulate=TRUE), 4)

}
\keyword{datasets}
