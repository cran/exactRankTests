
n <- 10
p1 <- c()

for (i in 1:1000) {
x <- rnorm(n)
y <- rnorm(n)   
N <- length(x) + length(y)
r <- rank(c(x,y))
scores <- qnorm(r/(N+1))
scores <- scores - min(scores)
scores <- round(scores*N/max(scores))

X <- sum(scores[seq(along=x)])
p2 <- pperm(X, scores, length(x))
p1 <- c(p1,2*ifelse(p2 > 0.5, 1 - p2, p2))
print(i)

}


