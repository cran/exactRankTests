
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
	seppl <- c(seppl, pperm(val, score, m))
	cat(" pperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

if (max(abs(hansi - seppl)) > 1e-10) stop("difference pperm to big!")

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

if (max(abs(hansi - seppl)) > 1e-10) stop("difference qperm to big!")

