\name{ears}
\alias{ears}
\docType{data}
\title{ Survival of Ventilating Tubes }
\description{
  Survival times of ventilating tubes of left and right ears in 78 children
with otitis media.
}
\usage{data(ears)}
\format{
  A data frame with 78 observations on the following 5 variables.
  \describe{
    \item{left}{Survival time in month of tube in left ear.}
    \item{lcens}{Censoring indicator for left ear: \code{0} censored and
                 \code{1} event.}
    \item{right}{Survival time in month of tube in right ear.}
    \item{rcens}{Censoring indicator for right ear: \code{0} censored and
                 \code{1} event.}
    \item{group}{a factor with levels \code{control} and \code{treat}. }}
}
\source{
  Sin-Ho Jung and Jong-Hyeon Jeong (2003). Rank tests for clustered survival
data. \emph{Lifetime Data Analysis}, \bold{9}, 21-33.
}
\references{
  V.M. Howie and R.H. Schwarz (1983). Acute otitis media: One year in
general pediatric practice. \emph{American Journal of Diseases in Children},
\bold{137}, 155-158.

  D.W. Teele, J.O. Klein, B. Rosner et al. (1989). Epidemiology of otitis
media during the first seven years of life in children in greater Boston.
\emph{Journal of Infectious Diseases}, \bold{160}, 89-94.

}
\examples{
data(ears)
if (require(survival, quietly=TRUE)) {
  ls <- cscores(Surv(ears$left, ears$lcens), int=TRUE)
  perm.test(ls ~ group, data=ears)
}

}
\keyword{datasets}
