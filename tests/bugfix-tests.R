
library(exactRankTests)

# From Ruggero Bellio <ruggero.bellio@dss.uniud.it>
# 20.02.2002

x <- c(0,87, 0, 0, 0, 0, 8, 0,16, 0,64, 0, 0,37, 200, 0,0,32, 5 ,0, 0,14, 0,
       0, 4, 0, 0,50, 0, 0,99, 0, 0,12,36,22, 4,13, 0,10,70, 0, 0, 0, 0,19,
       0, 4, 3, 4,23, 0, 0, 0,0, 0, 0 ,14, 4, 8, 0,11,42, 0, 0, 0, 0, 0,18,
       11,29, 0, 0)

y <- c(0, 0,14, 0, 0,16, 0, 145, 0, 0,16, 0, 0,47, 0, 0,0, 144, 0 , 2, 0,13,
       0, 0, 114, 0, 6, 0, 0, 4, 9, 0, 0,17, 0,0,42, 0 , 0, 0, 0, 0, 0, 0,25,
       0, 1, 178, 0, 0, 0, 6, 0,7,76, 252, 0 , 0,23, 0, 0,68, 0, 0, 0,24, 0,
       59, 0, 0, 0, 0,16)

wilcox.exact(x, y, paired=TRUE, conf.int=TRUE)  # was broken

p1 <- wilcox.exact(x, y, paired=TRUE, alternative="g", mu=-50)$p.value
p2 <- wilcox.exact(x, y, paired=TRUE, alternative="l", mu=-50)$p.value
p <- p1 + p2
names(p) <- NULL

stopifnot(all.equal(p, 1))

p1 <- perm.test(x, y, paired=TRUE, alternative="g")$p.value
p2 <- perm.test(x, y, paired=TRUE, alternative="l")$p.value
p <- p1 + p2
names(p) <- NULL

stopifnot(all.equal(p, 1))
