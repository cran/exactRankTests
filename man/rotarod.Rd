\name{rotarod}
\alias{rotarod}
\docType{data}
\title{ Rotating Rats Data}
\description{
  The endurance time of 24 rats in two groups in a rotating cylinder.
}
\usage{data(rotarod)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{time}{the endurance time}
    \item{group}{a factor at levels \code{control} and \code{treatment}.}}
}
\details{
  The rats were randomly assigned to the control and treatment group. The
treatment group received an oral dose of a centrally ating muscle relaxant
and the control group a saline solvent. This dataset serves as the basis of
an comparision of the results of the 
Wilcoxon-Mann-Whitney test computed by 11 statistical packages in Bergmann
et al. (2000). The exact conditional p-value is $0.0373$ (two-sided) and
$0.0186$ (one-sided). The asymptotic two-sided 
p-value (corrected for ties) is reported as $0.0147$.
}
\source{
  R. Bergmann, J. Ludbrook and W. P.J.M. Spooren (2000). Different outcomes
of the Wilcoxon-Mann-Whitney test from different statistics packages.
\emph{The American Statistician}, \bold{54}(1), 72--77.
}
\examples{
data(rotarod)
wilcox.exact(time ~ group, data=rotarod, alternative="g")
wilcox.exact(time ~ group, data=rotarod, conf.int=TRUE)
wilcox.exact(time ~ group, data=rotarod, exact=FALSE)
# the permutation test
perm.test(time ~ group, data=rotarod)
perm.test(time ~ group, data=rotarod, exact=FALSE)
}
\keyword{datasets}
