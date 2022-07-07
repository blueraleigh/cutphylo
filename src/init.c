#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP C_cutphylo(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_cutphylo, 9),
    {NULL, NULL, 0}
};

void attribute_visible R_init_cutphylo(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
