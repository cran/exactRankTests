\name{glioma}
\alias{glioma}
\docType{data}
\title{ Malignant Glioma Pilot Study}
\description{
  A non-randomized pilot study on malignant glioma patients with 
  pretargeted adjuvant radioimmunotherapy using Yttrium-90-biotin.
}
\usage{data(glioma)}
\format{
  A data frame with 37 observations on the following 7 variables.
  \describe{
    \item{No.}{patient number.}
    \item{Age}{patients ages in years.}
    \item{Sex}{a factor with levels \code{F}(emale) and \code{M}(ale). }
    \item{Histology}{a factor with levels \code{GBM} (grade IV) and 
                     \code{Grade3} (grade III)}
    \item{Survival}{survival times in month.}
    \item{Cens}{censoring indicator: \code{0} censored and \code{1} dead.}
    \item{Group}{a factor with levels \code{Control} and \code{RIT}.}
  }
}
\details{
  The primary endpoint of this small pilot study is survival. 
Survival times are tied, the usual asymptotic log-rank test may be
inadequate in this setup.
Therefore, a permutation test (via Monte-Carlo sampling) was conducted in
the original paper. The data are taken from Tables 1 and 2 of Grana et al.
(2002).
}
\source{

    C. Grana, M. Chinol, C. Robertson, C. Mazzetta, M. Bartolomei, C. De
    Cicco, M. Fiorenza, M. Gatti, P. Caliceti & G. Paganelli (2002),
    Pretargeted adjuvant radioimmunotherapy with Yttrium-90-biotin in malignant
    glioma patients: A pilot study. \emph{British Journal of Cancer},
    \bold{86}(2), 207--212.
}
\examples{
data(glioma)

if(require(survival, quietly = TRUE)) {

  par(mfrow=c(1,2))

  # Grade III glioma
  g3 <- glioma[glioma$Histology == "Grade3",]

  # Plot Kaplan-Meier curves
  plot(survfit(Surv(Survival, Cens) ~ Group, data=g3), 
       main="Grade III Glioma", lty=c(2,1), 
       legend.text=c("Control", "Treated"),
       legend.bty=1, ylab="Probability", 
       xlab="Survival Time in Month")

  # log-rank test
  survdiff(Surv(Survival, Cens) ~ Group, data=g3)

  # permutation test with integer valued log-rank scores
  lsc <- cscores(Surv(g3$Survival, g3$Cens), int=TRUE) 
  perm.test(lsc ~ Group, data=g3) 

  # permutation test with real valued log-rank scores
  lsc <- cscores(Surv(g3$Survival, g3$Cens), int=FALSE)
  tr <- (g3$Group == "RIT")
  T <- sum(lsc[tr])
  pperm(T, lsc, sum(tr), alternative="tw")
  pperm(T, lsc, sum(tr), alternative="tw", simulate=TRUE)

  # Grade IV glioma
  gbm <- glioma[glioma$Histology == "GBM",] 

  # Plot Kaplan-Meier curves
  plot(survfit(Surv(Survival, Cens) ~ Group, data=gbm), 
       main="Grade IV Glioma", lty=c(2,1), 
       legend.text=c("Control", "Treated"),
       legend.bty=1, legend.pos=1, ylab="Probability", 
       xlab="Survival Time in Month")
   
  # log-rank test
  survdiff(Surv(Survival, Cens) ~ Group, data=gbm)

  # permutation test with integer valued log-rank scores
  lsc <- cscores(Surv(gbm$Survival, gbm$Cens), int=TRUE)
  perm.test(lsc ~ Group, data=gbm)

  # permutation test with real valued log-rank scores 
  lsc <- cscores(Surv(gbm$Survival, gbm$Cens), int=FALSE) 
  tr <- (gbm$Group == "RIT")
  T <- sum(lsc[tr])
  pperm(T, lsc, sum(tr), alternative="tw")
  pperm(T, lsc, sum(tr), alternative="tw", simulate=TRUE)
}
}
\keyword{datasets}
