attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", pos.to.env(2), pos = length(search()))
assign("ptime", proc.time(), env = .CheckExEnv)
postscript("exactRankTests-Examples.ps")
assign("par.postscript", par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
library('exactRankTests')
rm(list = ls(all = TRUE)); .Random.seed <- c(0,rep(7654,3))
###--- >>> `dperm' <<<----- Distribution of Permutation Tests

	## alias	 help(dperm)
	## alias	 help(pperm)
	## alias	 help(pperm2)
	## alias	 help(qperm)
	## alias	 help(rperm)

##___ Examples ___:


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

p1 <- pperm2(X, scores, length(x))

# generate integer valued scores with the same shape as normal quantile
# scores, this no longer v.d.Waerden, but something very similar

scores <- scores - min(scores)
scores <- round(scores*N/max(scores))

X <- sum(scores[seq(along=x)])
p2 <- pperm2(X, scores, length(x))

# compare p1 and p2

p1 - p2

# the blood pressure example from StatXact:

treat <- c(94, 108, 110, 90)
contr <- c(80, 94, 85, 90, 90, 90, 108, 94, 78, 105, 88)

# compute the v.d. Waerden test and compare the results to StatXact:

r <- rank(c(contr, treat))
sc <- qnorm(r/16)
X <- sum(sc[seq(along=contr)])
round(pperm(X, sc, 11), 4)      # == 0.0462 (StatXact)
round(pperm2(X, sc, 11), 4)     # == 0.0799 (StatXact)

# the alternative method returns:

sc <- sc - min(sc)
sc <- round(sc*16/max(sc))
X <- sum(sc[seq(along=contr)])

round(pperm(X, sc, 11), 4)      # compare to 0.0462 
round(pperm2(X, sc, 11), 4)     # compare to 0.0799

# check if [qp]wilcox and [pq]signrank are equal with [qp]perm

source(system.file("test1.R", pkg="exactRankTests")[1])
source(system.file("test2.R", pkg="exactRankTests")[1])


## Keywords: 'exact distribution'.


rm(list = ls(all = TRUE)); .Random.seed <- c(0,rep(7654,3))
###--- >>> `wilcox.exact' <<<----- Wilcoxon Rank Sum and Signed Rank Tests

	## alias	 help(wilcox.exact)

##___ Examples ___:

## One-sample test.
## Hollander & Wolfe (1973), 29f.
## Hamilton depression scale factor measurements in 9 patients with
##  mixed anxiety and depression, taken at the first (x) and second
##  (y) visit after initiation of a therapy (administration of a
##  tranquilizer).
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
wilcox.exact(x, y, paired = TRUE, alternative = "greater")
wilcox.exact(y - x, alternative = "less")    # The same.

stopifnot(wilcox.test(y-x)$p.value == wilcox.exact(y -x)$p.value)

## Two-sample test.
## Hollander & Wolfe (1973), 69f.
## Permeability constants of the human chorioamnion (a placental
##  membrane) at term (x) and between 12 to 26 weeks gestational
##  age (y).  The alternative of interest is greater permeability
##  of the human chorioamnion for the term pregnancy.
x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
wilcox.exact(x, y, alternative = "g")        # greater

x <- rnorm(10)
y <- rnorm(10, 2)
wilcox.exact(x, y, conf.int = TRUE)

# Data from the StatXact-4 manual, page 221, diastolic blood pressure

treat <- c(94, 108, 110, 90)
contr <- c(80, 94, 85, 90, 90, 90, 108, 94, 78, 105, 88)

# StatXact 4 for Windows: p.value = 0.0989, point prob = 0.019

wilcox.exact(contr, treat)

# StatXact 4 for Windows: p.value = 0.0542, point prob = 0.019
 
wilcox.exact(contr, treat, alternative="less") 

# works in R-1.2.3 only
#  
# x <- rnorm(10,3)
# y <- rnorm(10)
# stopifnot(all.equal(wilcox.test(x,y,conf.int=T), wilcox.exact(x,y,conf.int=T)))


## Keywords: 'htest'.


cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
dev.off(); quit('no')
