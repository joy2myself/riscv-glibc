/* SLEEF-component composite comparison entry for log10p1f.
   SLEEF provides xlog10f/xlog1pf, but no native whole-routine log10p1f.  */
#define LIBMVEC_TYPE float
#define LIBMVEC_FUNC log10p1f_original
#define LIBMVEC_IMPL "v_f_log10p1f_original_skeleton.c"
#include "v_math_importer.h"
