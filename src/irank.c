/*

  $Id: irank.c,v 1.3 2002/11/18 10:33:30 hothorn Exp $
  
  irank: integer valued ranks
  modified from R-1.5.1/src/main/sort.c

*/
      
#include <R.h>
#include <Rinternals.h>


SEXP irank(SEXP x, SEXP orderx)
{
    SEXP rk;
    double *tx;
    int *to;
    int i, j, k, n;

    if (!isVector(x))
	error("Argument is not a vector");
    n = LENGTH(x);
    if (!isVector(orderx)) /* || length(orderx) != n) */
        error("orderx is not a vector of the same length as x");
    PROTECT(rk = allocVector(INTSXP, n));
    UNPROTECT(1);
    if (n > 0) {
        tx = REAL(x);
        to = INTEGER(orderx);
	i = 0;
	while (i < n) {
	    j = i;
	    while ((j < n - 1) && (tx[to[j]] == tx[to[j+1]]))
		j++;
	    if (i != j) {
		for (k = i; k <= j; k++)
		    /*
		    return the number of observations 
		    less or equal tx[to[k]]
		    */
		    INTEGER(rk)[to[k]] = j + 1;
	    }
	    else
		INTEGER(rk)[to[i]] = i + 1;
	    i = j + 1;
	}
    }
    return rk;
}
