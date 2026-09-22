/* SLEEF-component composite comparison entry for asinpif.
   SLEEF provides xasinf_u1, but no native whole-routine asinpif entry.  */
#define LIBMVEC_TYPE float
#define LIBMVEC_FUNC asinpif_original
#define LIBMVEC_IMPL "v_f_asinpif_original_skeleton.c"
#include "v_math_importer.h"
