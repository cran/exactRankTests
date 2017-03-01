
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cpermdist1(SEXP);
extern SEXP cpermdist2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_irank(SEXP, SEXP);
extern SEXP sim2is(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cpermdist1", (DL_FUNC) &cpermdist1, 1},
    {"cpermdist2", (DL_FUNC) &cpermdist2, 5},
    {"C_irank",    (DL_FUNC) &C_irank,    2},
    {"sim2is",     (DL_FUNC) &sim2is,     3},
    {NULL, NULL, 0}
};

void R_init_exactRankTests(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
