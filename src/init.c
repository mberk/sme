#include "SME.h"

#include <R_ext/Rdynload.h>
#include <Rinternals.h>

static R_NativePrimitiveArgType SME_t[] =
{
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP
};

static R_NativePrimitiveArgType SMEMultiple_t[] =
{
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP
};

static R_NativePrimitiveArgType SMEOptimization_t[] =
{
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP
};

static R_NativePrimitiveArgType SMEOptimizationMultiple_t[] =
{
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};

static const R_CMethodDef cMethods[] =
{
  {"SME", (DL_FUNC) &SME, 23, SME_t},
  {"SMEMultiple", (DL_FUNC) &SMEMultiple, 26, SMEMultiple_t},
  {"SMEOptimization", (DL_FUNC) &SMEOptimization, 25, SMEOptimization_t},
  {"SMEOptimizationMultiple", (DL_FUNC) &SMEOptimizationMultiple, 27, SMEOptimizationMultiple_t},
  {NULL, NULL, 0, NULL}
};

void R_init_sme(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
