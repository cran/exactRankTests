dansari <- function(x, m ,n)
{
	.C("dansari", as.integer(length(x)), d = as.double(x), as.integer(m),
	as.integer(n), PACKAGE="ctest")$d
}
