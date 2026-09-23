/* Source-original SLEEF comparison entry for fastsinf.  */
#define LIBMVEC_TYPE float
#define LIBMVEC_FUNC fastsinf_original
#define LIBMVEC_IMPL "v_f_fastsinf_original_skeleton.c"
/* Shared SLEEF trigonometric helpers for all LMUL instances.  */
#define RVV_NEED_SLEEF_F32_TRIG_REDUCTION 1
#include "v_math_importer.h"
