/* SLEEF-component composite comparison entry for log2p1f.
   SLEEF provides xlog2f/xlog1pf, but no native whole-routine log2p1f.  */
#define LIBMVEC_TYPE float
#define LIBMVEC_FUNC log2p1f_original
#define LIBMVEC_IMPL "v_f_log2p1f_original_skeleton.c"
#include "v_math_importer.h"
