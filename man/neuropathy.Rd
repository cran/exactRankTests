\name{neuropathy}
\alias{neuropathy}
\docType{data}
\title{ Acute Painful Diabetic Neuropathy }
\description{
  The logarithm of the ratio of pain scores at baseline and after four weeks
for a control and treatment group.
}
\usage{data(neuropathy)}
\format{
  A data frame with 58 observations on the following 2 variables.
  \describe{
    \item{pain}{Pain scores: ln(baseline/final).}
    \item{group}{a factor with levels \code{control} and \code{treat}.} }
}
\details{
  Data from Table 1 of Conover & Salsburg (1988).
}
\source{

        William J. Conover and David S. Salsburg (1988), Locally most powerful tests
        for detecting treatment effects when only a subset of patients can  
        be expected to "respond" to treatment. \emph{Biometrics}, \bold{44},
        189--196.

}
\examples{
data(neuropathy)
# compare with Table 2 of Conover & Salsburg (1988)
wilcox.exact(pain ~ group, data=neuropathy, alternative="less")
css <- cscores(neuropathy$pain, type="ConSal")
pperm(sum(css[neuropathy$group=="control"]),css,
      m=sum(neuropathy$group=="control"))

}
\keyword{datasets}
