\name{dperm}
\alias{dperm}
\alias{pperm}
\alias{qperm}
\alias{rperm}
\title{ Distribution of One and Two Sample Permutation Tests }
\description{

Density, distribution function and quantile function 
for the distribution of one and two sample 
permutation tests using the Shift-Algorithm by
Streitberg & R\"ohmel.

}
\usage{
dperm(x, scores, m, paired=NULL, tol = 0.01, fact=NULL, density=FALSE,
      simulate=FALSE, B=10000)
pperm(q, scores, m, paired=NULL, tol = 0.01, fact=NULL,
      alternative=c("less", "greater", "two.sided"), pprob=FALSE,
      simulate=FALSE, B=10000)
qperm(p, scores, m, paired=NULL, tol = 0.01, fact=NULL, 
      simulate=FALSE, B=10000)
rperm(n, scores, m)
}
\arguments{
  \item{x, q}{ vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{scores}{arbitrary scores of the observations 
                of the \code{x} (first \code{m} elements) and 
                \code{y} sample. }
  \item{m}{sample size of the \code{x} sample. If \code{m = length(x)}
           scores of paired observations are assumed.}
  \item{paired}{logical. Indicates if paired observations are used. Needed
                to discriminate between a paired problem and the distribution 
                of the total sum of the scores (which has mass 1 at the 
                point \code{sum(scores)}).}.
  \item{tol}{real. Real valued scores are mapped into integers by rounding
             after multiplication with an appropriate factor.
             Make sure that the absolute difference between the 
             each possible test statistic for the original scores and the
             rounded scores is less than \code{tol}. This might not be 
             possible due to memory/time limitations, a warning is given in
             this case.} 
  \item{fact}{real. If \code{fact} is given, real valued scores are mapped 
              into integers using \code{fact} as factor. \code{tol} is 
              ignored in this case.}
  \item{n}{number of random observations to generate.}
  \item{alternative}{character indicating whether the probability 
                     \eqn{P(T \le q)} (\code{less}), 
                     \eqn{P(T \ge q)} (\code{greater}) or a two-sided 
                     p-value (\code{two.sided}) should be computed in
                     \code{pperm}.} 
  \item{pprob}{logical. Indicates if the probability \eqn{P(T = q)} should
               be computed additionally.}
  \item{density}{logical. When  \code{x} is a scalar and \code{density} is
                 \code{TRUE}, \code{dperm} returns the density for all
                 possible statistics less or equal \code{x} as a data frame.}
  \item{simulate}{logical. Use conditional Monte-Carlo to 
                  compute the distribution.}
  \item{B}{number of Monte-Carlo replications to be used.}
}
\details{

  The exact distribution of the sum of the first \code{m} scores is
evaluated using the Shift-Algorithm by Streitberg & R\"ohmel under the
hypothesis of exchangeability (or, equivalent, the hypothesis that all
permutations of the scores are equally likely). The algorithm is able
to deal with tied scores, so the conditional distribution can be
evaluated. 

The algorithm is defined for positive integer valued scores only. 
There are two ways dealing with real valued scores. 
First, one can try to find integer valued scores that lead to statistics 
which differ not more than \code{tol}
from the statistics computed for the original scores. This can be done as
follows.  

Without loss of generality let \eqn{a_i > 0} denote real valued scores in
reverse ordering and
\eqn{f} a positive factor (this is the \code{fact} argument). 
Let \eqn{R_i = f \cdot a_i - round(f \cdot a_i)}.  Then 

\deqn{ \sum_{i=1}^m f \cdot a_i = \sum_{i=1}^m round(f \cdot a_i) - R_i. }

Clearly, the maximum difference between \eqn{1/f \sum_{i=1}^m f \cdot a_i} and
\eqn{1/f \sum_{i=1}^n round(f \cdot a_i)} is given by 
\eqn{|\sum_{i=1}^m R_i|}. Therefore one searches for \eqn{f} with 

\deqn{ |\sum_{i=1}^m R_i| \le \sum_{i=1}^m |R_i| \le tol.}

If \eqn{f} induces more that 100.000 columns in the Shift-Algorithm by 
Streitberg  & R\"ohmel, \eqn{f} is restricted to the largest integer 
that does not. 

The second idea is to map the scores into integers by taking the
integer part of \eqn{a_i N / \max(a_i)} (Hothorn & Lausen, 2002). 
This induces additional ties, but the shape of the
scores is very similar. That means we do not try to approximate something
but use a different test (with integer valued scores), serving for the 
same purpose
(due to a similar shape of the scores). However, this has to be done prior
to calling \code{pperm} (see the examples).

Exact two-sided p-values are computed as suggested in the StatXact-5 
manual, page 225, equation (9.31) and equation (8.18), p. 179 (paired case).
In detail: For the paired case the two-sided p-value is just twice the
one-sided one. For the independent sample case the two sided p-value is
defined as
\deqn{p_2 = P( |T - E(T)| \ge | q - E(T) |)}
where \eqn{q} is the quantile passed to \code{pperm}.

}
\value{

\code{dperm} gives the density, \code{pperm} gives the distribution
function and \code{qperm} gives the quantile function. If \code{pprob} is
true, \code{pperm} returns a list with elements
\item{PVALUE}{the probability specified by \code{alternative}.}
\item{PPROB}{the probability \eqn{P(T = q)}.}
\code{rperm} is a wrapper to \code{\link{sample}}.

}
\references{ 

	Bernd Streitberg & Joachim R\"ohmel (1986),
	Exact distributions for permutations and rank tests: 
	An introduction to some recently published algorithms. 
	\emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

	Bernd Streitberg & Joachim R\"ohmel (1987),
	Exakte Verteilungen f\"ur Rang- und Randomisierungstests 
	im allgemeinen $c$-Stichprobenfall.
	\emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

	Torsten Hothorn (2001),
	On exact rank tests in \R.
	\emph{R News} \bold{1}(1), 11--12.
	
        Cyrus R. Mehta & Nitin R. Patel (2001),
        \emph{StatXact-5 for Windows.}
        Manual, Cytel Software Cooperation, Cambridge, USA

        Torsten Hothorn & Berthold Lausen (2002),
        On the exact distribution of maximally selected rank statistics.
        \emph{Computational Statistics \& Data Analysis}, (in press).
        Preprint available from
        \url{http://www.mathpreprints.com/math/Preprint/hothorn/20020227/2}



}

\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

# exact one-sided p-value of the Wilcoxon test for a tied sample

x <- c(0.5, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9)
y <- c(0.5, 1.0, 1.2, 1.2, 1.4, 1.5, 1.9, 2.0)
r <- cscores(c(x,y), type="Wilcoxon")
pperm(sum(r[seq(along=x)]), r, 7)

# Compare the exact algorithm as implemented in ctest and the
# Shift-Algorithm by Streitberg & Roehmel for untied samples
 
# Wilcoxon:

n <- 10
x <- rnorm(n, 2)
y <- rnorm(n, 3)
r <- cscores(c(x,y), type="Wilcoxon")

# exact distribution using the Shift-Algorithm

dwexac <- dperm((n*(n+1)/2):(n^2 + n*(n+1)/2), r, n)
sum(dwexac)           # should be something near 1 :-)

# exact distribution using dwilcox

dw <- dwilcox(0:(n^2), n, n)

# compare the two distributions:

plot(dw, dwexac, main="Wilcoxon", xlab="dwilcox", ylab="dperm")      
# should give a "perfect" line

# Wilcoxon signed rank test

n <- 10
x <- rnorm(n, 5)
y <- rnorm(n, 5)
r <- cscores(abs(x - y), type="Wilcoxon")
pperm(sum(r[x - y > 0]), r, length(r))
wilcox.test(x,y, paired=TRUE, alternative="less")
psignrank(sum(r[x - y > 0]), length(r))

# Ansari-Bradley

n <- 10
x <- rnorm(n, 2, 1)
y <- rnorm(n, 2, 2)

# exact distribution using the Shift-Algorithm

sc <- cscores(c(x,y), type="Ansari")
dabexac <- dperm(0:(n*(2*n+1)/2), sc, n)
sum(dabexac)

# real scores are allowed (but only result in an approximation)
# e.g. v.d. Waerden test

n <- 10
x <- rnorm(n)
y <- rnorm(n)
scores <- cscores(c(x,y), type="NormalQuantile")
X <- sum(scores[seq(along=x)])  # <- v.d. Waerden normal quantile statistic

# critical value, two-sided test

abs(qperm(0.025, scores, length(x)))

# p-values

p1 <- pperm(X, scores, length(x), alternative="two.sided")

# generate integer valued scores with the same shape as normal quantile
# scores, this no longer v.d.Waerden, but something very similar

scores <- cscores(c(x,y), type="NormalQuantile", int=TRUE)

X <- sum(scores[seq(along=x)])
p2 <- pperm(X, scores, length(x), alternative="two.sided")

# compare p1 and p2

p1 - p2

}
\keyword{distribution}
