\name{ASAT}
\alias{ASAT}
\docType{data}
\title{ Toxicological Study on Female Wistar Rats }
\description{
  ASAT-Values for a new compound and a control group of 34 female Wistar
rats.
}
\usage{data(ASAT)}
\format{
  A data frame with 34 observations on the following 2 variables.
  \describe{
    \item{asat}{the ASAT-values (a liver enzyme)}
    \item{group}{a factor at levels \code{Compound} and \code{Control}.}
  }
}
\details{
  The aim of this toxicological study is the proof of safety for the new
compound. The data are originally given in Hothorn (1992) and reproduced in
Hauschke et al. (1999).
}
\source{
  L. A. Hothorn (1992), Biometrische Analyse toxikologischer Untersuchungen.
  In: J. Adams (ed.): \emph{Statistisches Know how in der medizinischen 
  Forschung.} Ullstein-Mosby, Berlin, 475--590.
}
\references{
  D. Hauschke, M. Kieser & L. A. Hothorn (1999), Proof of safety in
  toxicology based on the ratio of two means for normally distributed data.
  \emph{Biometrical Journal}, \bold{41}(3), 295--304.
}
\examples{
data(ASAT)
# does not really look symmetric
plot(asat ~ group, data=ASAT)
# proof-of-safety based on ration of medians
pos <- wilcox.exact(I(log(asat)) ~ group, data = ASAT, alternative = "less", 
             conf.int=TRUE)
# one-sided confidence set. Safety cannot be concluded since the effect of
# the compound exceeds 20\% of the control median
exp(pos$conf.int)
}
\keyword{datasets}
