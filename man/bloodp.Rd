\name{bloodp}
\alias{bloodp}
\docType{data}
\title{ Diastolic Blood Pressure }
\description{
  Diastolic blood pressure for a two groups of patients.
}
\usage{data(bloodp)}
\format{
  A data frame with 15 observations on the following 2 variables.
  \describe{
    \item{bp}{the diastolic blood pressure.}
    \item{group}{a factor with levels \code{group1} and \code{group2}.}
  }
}
\details{
  The data is given in Table 9.6, page 227, of Metha and Pathel (2001). Note
that there are some tied observations. The permutation test using the raw
blood pressure values does not lead to a rejection of the null hypothesis of
exchangeability: p-value = 0.1040 (two-sided) and p-value = 0.0564
(one-sided). The asymptotic two-sided p-value is 0.1070.

For the Wilcoxon-Mann-Whitney test,
the one-sided p-value is 0.0542 and the two-sided one is 0.0989 (Metha &
Patel, 2001, page 229).

The one-sided p-value for the v.d.Waeren test is 0.0462 (Metha &
Patel, 2001, page 241) and the two-sided p-value is 0.0799.

}
\references{

        Cyrus R. Mehta & Nitin R. Patel (2001),
        \emph{StatXact-5 for Windows.}
        Manual, Cytel Software Cooperation, Cambridge, USA

}
\examples{
data(bloodp)

# Permutation test

perm.test(bp ~ group, data=bloodp)
perm.test(bp ~ group, data=bloodp, alternative="greater")
perm.test(bp ~ group, data=bloodp, exact=FALSE)

# Wilcoxon-Mann-Whitney test

wilcox.exact(bp ~ group, data=bloodp, conf.int=TRUE, alternative="l")
wilcox.exact(bp ~ group, data=bloodp, conf.int=TRUE)

# compute the v.d. Waerden test

sc <- cscores(bloodp$bp, type="NormalQuantile")
X <- sum(sc[bloodp$group == "group2"])
round(pperm(X, sc, 11), 4) 
round(pperm(X, sc, 11, simulate=TRUE), 4)
round(pperm(X, sc, 11, alternative="two.sided"), 4) 
round(pperm(X, sc, 11, alternative="two.sided", simulate=TRUE), 4)

# use scores mapped into integers (cf. dperm)

sc <- cscores(bloodp$bp, type="NormalQuantile", int=TRUE)
X <- sum(sc[bloodp$group == "group2"])
round(pperm(X, sc, 11), 4)      
round(pperm(X, sc, 11, alternative="two.sided"), 4)

}
\keyword{datasets}
