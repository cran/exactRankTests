\name{dperm}
\alias{dperm}
\alias{pperm}
\alias{qperm}
\alias{rperm}
\title{Distribution of One and Two Sample Permutation Tests}
\description{
Density, distribution function and quantile function 
for the distribution of one and two sample 
permutation tests using the Shift-Algorithm by
Streitberg and R\"ohmel.
}
\usage{
dperm(x, scores, m, paired=NULL, tol = 0.01, fact=NULL)
pperm(q, scores, m, paired=NULL, tol = 0.01, fact=NULL,
      alternative=c("less", "greater", "two.sided"), pprob=FALSE)
qperm(p, scores, m, paired=NULL, tol = 0.01, fact=NULL)
rperm(n, scores, m)
}
\arguments{
  \item{x, q}{ vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{scores}{ranks, midranks or real valued scores of the observations 
of the \code{x} (first \code{m} elements) and \code{y} sample. }
  \item{m}{ sample size of the \code{x} sample. If \code{m = length(x)}
scores of paired observations are assumed.}
  \item{paired}{ logical. Indicates if paired observations are used. Needed
to discriminate between a paired problem and the distribution of the total sum
of the scores (which has mass 1 at the point \code{sum(scores)})}.
  \item{tol}{ real. real valued scores are mapped into integers by multiplication.
Make sure that the absolute difference between the "true" quantile and the
approximated quantile is less than \code{tol}. This might not be possible
due to memory/time limitations.} 
  \item{fact}{ real. If \code{fact} is given, real valued scores are mapped into
integers using \code{fact} as factor. \code{tol} is ignored.}
  \item{n}{number of observations.}
  \item{alternative}{character indicating whether the probability \eqn{P(T \le q)}
(\code{less}), \eqn{P(T \ge q)} (\code{greater}) or a two sided p-value
(\code{two.sided}) should be computed in
\code{pperm}.} 
  \item{pprob}{logical. Indicates if the probability \eqn{P(T = q)} should
be computed additionally.}
}
\details{
  The exact distribution of the sum of the first \code{m} scores is
evaluated using the Shift-Algorithm by Streitberg and R\"ohmel under the
hypothesis of exchangeability (or, equivalent, the hypothesis that all
permutations of the scores are equally likely). The algorithm is able
to deal with tied scores, so the conditional distribution can be
evaluated. 

The algorithm is defined for positive integer valued scores only. 
There are two ways dealing with real valued scores. 
First, one can try to find integer
valued scores that lead to quantiles which differ not more than \code{tol}
from the quantiles induced by the original scores. This can be done as
follows.  

Without loss of generality let \eqn{a_i > 0} denote real valued scores and
\eqn{f} a positive factor. Let \eqn{R_i = a_i - round(f \cdot a_i)}.  
Then 

\deqn{ \sum_{i=1}^m f \cdot a_i = \sum_{i=1}^m round(f \cdot a_i) - R_i. }

Clearly, the maximum difference between \eqn{\sum_{i=1}^m f \cdot a_i} and
\eqn{\sum_{i=1}^n round(f \cdot a_i)} is given by \eqn{|\sum_{i=1}^m R_{(i)}|} or
\eqn{|\sum_{i=m+1}^N R_{(i)}|}, respectively. Therefore one searches for
\eqn{f} with 

\deqn{ \max(|\sum_{i=1}^m R_{(i)}|, |\sum_{i=m+1}^N R_{(i)}|)/f \le tol.}

If \eqn{f} induces more that 20.000 columns in the Streitberg-R\"ohmel
Shift-Algorithm, \eqn{f} is restricted to the largest integer that does not. 

The second idea is to map the scores itself into
\eqn{\{1, \dots, N\}}. This induces additional ties, but the shape of the
scores is very similar. That means we do not try to approximate something
but use a different test (with integer scores), serving for the same purpose
(due to a similar shape of the scores). However, this has to be done prior
to calling \code{pperm} (see the examples).

Exact two-sided p-values are computed as suggested in the StatXact-4 
manual, page 209, equation (9.32) and equation (8.19), p. 165 (paired case).
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

	Bernd Streitberg and Joachim R\"ohmel (1986). 
	Exact Distributions for Permutations and Rank Tests: 
	An Introduction to Some Recently Published Algorithms. 
	\emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

	Bernd Streitberg and Joachim R\"ohmel (1987).
	Exakte Verteilungen f\"ur Rang- und Randomisierungstests 
	im allgemeinen $c$-Stichprobenfall.
	\emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

	Torsten Hothorn (2001).
	On Exact Rank Tests in \R.
	\emph{R News} \bold{1}(1), 11--12.
	
	Cyrus R. Mehta and Nitin R. Patel (1998).
	\emph{StatXact-4 for Windows.}
	Manual, Cytel Software Cooperation, Cambridge, USA

}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

# exact one-sided p-value of the Wilcoxon test for a tied sample

x <- c(0.5, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9)
y <- c(0.5, 1.0, 1.2, 1.2, 1.4, 1.5, 1.9, 2.0)
r <- rank(c(x,y))
pperm(sum(r[seq(along=x)]), r, 7)

# Compare the exact algorithm as implemented in ctest and the
# Streitberg-Roehmel for untied samples
 
# Wilcoxon:

n <- 10
x <- rnorm(n, 2)
y <- rnorm(n, 3)
r <- rank(c(x,y))

# exact distribution using Streitberg-Roehmel

dwexac <- dperm((n*(n+1)/2):(n^2 + n*(n+1)/2), r, n)
su <- sum(dwexac)           # should be something near 1 :-)
su
if (su != 1) stop("sum(dwexac) not equal 1")

# exact distribution using dwilcox

dw <- dwilcox(0:(n^2), n, n)

# compare the two distributions:

plot(dw, dwexac, main="Wilcoxon", xlab="dwilcox", ylab="dperm")      
# should give a "perfect" line

# Wilcoxon signed rank test

n <- 10
x <- rnorm(n, 5)
y <- rnorm(n, 5)
r <- rank(abs(x - y))
pperm(sum(r[x - y > 0]), r, length(r))
wilcox.test(x,y, paired=T, alternative="less")
psignrank(sum(r[x - y > 0]), length(r))

# Ansari-Bradley

n <- 10
x <- rnorm(n, 2, 1)
y <- rnorm(n, 2, 2)

# exact distribution using Streitberg-Roehmel

r <- rank(c(x,y))
sc <- pmin(r, 2*n - r +1)
dabexac <- dperm(0:(n*(2*n+1)/2), sc, n)
sum(dabexac)
tr <- which(dabexac > 0)

# exact distribution using dansari (wrapper to ansari.c in ctest)

dab <- dansari(0:(n*(2*n+1)/2), n, n)

# compare the two distributions:

plot(dab[tr], dabexac[tr], main="Ansari", xlab="dansari", ylab="dperm")

# real scores are allowed (but only result in an approximation)
# e.g. v.d. Waerden test

n <- 10
x <- rnorm(n)
y <- rnorm(n)
N <- length(x) + length(y)
r <- rank(c(x,y))
scores <- qnorm(r/(N+1))
X <- sum(scores[seq(along=x)])  # <- v.d. Waerden normal quantile statistic

# critical value, two-sided test

abs(qperm(0.025, scores, length(x)))

# p-values

p1 <- pperm(X, scores, length(x), alternative="two.sided")

# generate integer valued scores with the same shape as normal quantile
# scores, this no longer v.d.Waerden, but something very similar

scores <- scores - min(scores)
scores <- round(scores*N/max(scores))

X <- sum(scores[seq(along=x)])
p2 <- pperm(X, scores, length(x), alternative="two.sided")

# compare p1 and p2

p1 - p2

# the blood pressure example from StatXact manual, page 221:

treat <- c(94, 108, 110, 90)
contr <- c(80, 94, 85, 90, 90, 90, 108, 94, 78, 105, 88)

# compute the v.d. Waerden test and compare the results to StatXact-4 for
# Windows:

r <- rank(c(contr, treat))
sc <- qnorm(r/16)
X <- sum(sc[seq(along=contr)])
round(pperm(X, sc, 11), 4) 	# == 0.0462 (StatXact)
round(pperm(X, sc, 11, alternative="two.sided"), 4) 	# == 0.0799 (StatXact)

# the alternative method returns:

sc <- sc - min(sc)
sc <- round(sc*16/max(sc))
X <- sum(sc[seq(along=contr)])

round(pperm(X, sc, 11), 4)      # compare to 0.0462 
round(pperm(X, sc, 11, alternative="two.sided"), 4)     # compare to 0.0799

\testonly{

# paired observations

hansi <- c()
seppl <- c()
for (i in 1:10)
{
	m <- sample(10:50, 1)
	score <- sample(m)
	val <- sample(0:m, 1)
	# cat("m: ", m, "n: ", n, " val: ", val, "\n")
	hansi <- c(hansi,  psignrank(val, m))
	cat("psignrank: ", hansi[length(hansi)])
	seppl <- c(seppl, pperm(val, score, m, alternative="less"))
	cat(" pperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

\testonly{
  stopifnot(max(abs(hansi - seppl)) <= 1e-10)
}

hansi <- c()
seppl <- c()
for (i in 1:10)
{
        m <- sample(10:50, 1)
        score <- sample(m)
        prob <- runif(1)
        # cat("m: ", m, "n: ", n, " prob: ", prob, "\n")
        hansi <- c(hansi,  qsignrank(prob, m))
        cat("qwilcox: ", hansi[length(hansi)])
        seppl <- c(seppl, qperm(prob, score, m))
        cat(" qperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

\testonly{
  stopifnot(max(abs(hansi - seppl)) <= 1e-10)
}

# independent observations

hansi <- c()
seppl <- c()
for (i in 1:10)
{
	m <- sample(10:50, 1)
	if (runif(1) < 0.5)
		n <- sample(10:50, 1)
	else	
		n <- m
	score <- sample(n+m)
	val <- sample(0:(m*n), 1)
	# cat("m: ", m, "n: ", n, " val: ", val, "\n")
	hansi <- c(hansi,  pwilcox(val, m, n))
	cat("pwilcox: ", hansi[length(hansi)])
	seppl <- c(seppl, pperm(val + m*(m+1)/2, score, m,
alternative="less"))
	cat(" pperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

\testonly{
  stopifnot(max(abs(hansi - seppl)) <= 1e-10) 
}

hansi <- c()
seppl <- c()
for (i in 1:10)
{
        m <- sample(10:50, 1)
        if (runif(1) < 0.5)
                n <- sample(10:50, 1)
        else
                n <- m
        score <- sample(n+m)
        prob <- runif(1)
        # cat("m: ", m, "n: ", n, " prob: ", prob, "\n")
        hansi <- c(hansi,  qwilcox(prob, m, n))
        cat("qwilcox: ", hansi[length(hansi)])
        seppl <- c(seppl, qperm(prob, score, m) - m*(m+1)/2)
        cat(" qperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

\testonly{
  stopifnot(max(abs(hansi - seppl)) <= 1e-10) 
}

\testonly{
  # bugged me
  stopifnot(pperm(36, c(16,15,5,11,14), 5, alternative="t") == 0.625)
}

}

}
\keyword{distribution}
