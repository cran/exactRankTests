\name{irank}
\alias{irank}
\title{ Integer Ranks }
\description{
  Compute the number of elements less or equal the elements in a given vector.
}
\usage{
irank(x, ox = NULL)
}
\arguments{
  \item{x}{a numeric vector. }
  \item{ox}{\code{order(x)}, optionally (for efficiency in case
            \code{order(x)} is already known). }
}
\value{
  A vector of integers. 
}
\examples{
x <- rnorm(10)
irank(x)
rank(x)
x <- c(1,2,3,3,0)
irank(x)
rank(x)
}
\keyword{univar}
