/* Utilities for RVV libmvec routines.
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#ifndef _V_MATH_H
#define _V_MATH_H

#include <riscv_vector.h>
#include <math_private.h>
#include <assert.h>

#include <v_math_names.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#define _NEED_UNDEF_GNU_SOURCE_V_MATH
#endif

#include <math.h>
#ifndef NAN
_Static_assert (0, "NaN not available on this architecture");
#endif

#define __PASTE2_BASE(A, B) A##B
#define __PASTE2(A, B) __PASTE2_BASE (A, B)
#define __PASTE3_BASE(A, B, C) A##B##C
#define __PASTE3(A, B, C) __PASTE3_BASE (A, B, C)
#define __PASTE4_BASE(A, B, C, D) A##B##C##D
#define __PASTE4(A, B, C, D) __PASTE4_BASE (A, B, C, D)
#define __PASTE5_BASE(A, B, C, D, E) A##B##C##D##E
#define __PASTE5(A, B, C, D, E) __PASTE5_BASE (A, B, C, D, E)
static const uint64_t class_NaN = 0x300;
static const uint64_t class_posZero = 0x10;
static const uint64_t class_posDenorm = 0x20;
static const uint64_t class_negative = 0x7;
#define RVV_SLEEF_TRIG_PI_F 0x1.921fb6p+1f

#ifdef _NEED_UNDEF_GNU_SOURCE_V_MATH
#undef _GNU_SOURCE
#undef _NEED_UNDEF_GNU_SOURCE_V_MATH
#endif

#endif

#ifndef LMUL
# error "Must assign an LMUL before including v_math.h"
#endif
#ifndef BIT_WIDTH
# error "Must assign BIT_WIDTH before including v_math.h"
#endif

/* Undefine previous LMUL-dependent macros to allow redefinition */
#undef MAKE_VTYPE
#undef MAKE_TYPE
#undef MAKE_FUNC
#undef MAKE_VBOOL
#undef VSET
#undef MAKE_REINTERPRET
#undef FLOAT
#undef VFLOAT
#undef INT
#undef VINT
#undef UINT
#undef VUINT
#undef VBOOL
#undef F_AS_I
#undef F_AS_U
#undef I_AS_F
#undef U_AS_F
#undef I_AS_U
#undef U_AS_I
#undef VMVI_VX
#undef VFMV_VF
#undef EXP_BIAS
#undef MAN_LEN

#define MAKE_VTYPE(A) __PASTE3 (A, BIT_WIDTH, __PASTE3 (m, LMUL, _t))
#define MAKE_TYPE(A) __PASTE3 (A, BIT_WIDTH, _t)
#define MAKE_FUNC(A) __PASTE3 (A, BIT_WIDTH, __PASTE2 (m, LMUL))

#if (BIT_WIDTH == 64)
# define EXP_BIAS 1023
# define MAN_LEN 52
# ifndef _V_MATH_SCALAR_64
#  define _V_MATH_SCALAR_64
typedef double scalar_t;
static const scalar_t fp_sNaN = __builtin_nans ("");
static const scalar_t fp_negZero = -0.;
static const scalar_t fp_posZero = 0.0;
static const scalar_t fp_posOne = 0x1.0p0;
# endif
#elif (BIT_WIDTH == 32)
# define EXP_BIAS 127
# define MAN_LEN 23
# ifndef _V_MATH_SCALAR_32
#  define _V_MATH_SCALAR_32
typedef float scalar_t;
static const scalar_t fp_sNaN = __builtin_nansf ("");
static const scalar_t fp_negZero = -0.f;
static const scalar_t fp_posZero = 0.0f;
static const scalar_t fp_posOne = 0x1.0p0f;
# endif
#else
# error "requested BIT_WIDTH unsupported"
#endif

/* For MAKE_VBOOL, the value is BIT_WIDTH/LMUL for the current element size. */
#if (BIT_WIDTH == 64)
# if (LMUL == 1)
#  define MAKE_VBOOL(A) __PASTE3 (A, 64, _t)
# elif (LMUL == 2)
#  define MAKE_VBOOL(A) __PASTE3 (A, 32, _t)
# elif (LMUL == 4)
#  define MAKE_VBOOL(A) __PASTE3 (A, 16, _t)
# elif (LMUL == 8)
#  define MAKE_VBOOL(A) __PASTE3 (A, 8, _t)
# else
#  error "requested LMUL unsupported"
# endif
#elif (BIT_WIDTH == 32)
# if (LMUL == 1)
#  define MAKE_VBOOL(A) __PASTE3 (A, 32, _t)
# elif (LMUL == 2)
#  define MAKE_VBOOL(A) __PASTE3 (A, 16, _t)
# elif (LMUL == 4)
#  define MAKE_VBOOL(A) __PASTE3 (A, 8, _t)
# elif (LMUL == 8)
#  define MAKE_VBOOL(A) __PASTE3 (A, 4, _t)
# else
#  error "requested LMUL unsupported"
# endif
#else
# error "requested BIT_WIDTH unsupported"
#endif

#define VSET __PASTE2 (__riscv_vsetvl_e, __PASTE3 (BIT_WIDTH, m, LMUL))
#define MAKE_REINTERPRET(A, B)                                                \
  __PASTE5 (__riscv_vreinterpret_v_, A, __PASTE4 (BIT_WIDTH, m, LMUL, _), B,  \
            __PASTE3 (BIT_WIDTH, m, LMUL))

#define FLOAT MAKE_TYPE (float)
#define VFLOAT MAKE_VTYPE (vfloat)
#define INT MAKE_TYPE (int)
#define VINT MAKE_VTYPE (vint)
#define UINT MAKE_TYPE (uint)
#define VUINT MAKE_VTYPE (vuint)
#define VBOOL MAKE_VBOOL (vbool)

#define F_AS_I MAKE_REINTERPRET (f, i)
#define F_AS_U MAKE_REINTERPRET (f, u)
#define I_AS_F MAKE_REINTERPRET (i, f)
#define U_AS_F MAKE_REINTERPRET (u, f)
#define I_AS_U MAKE_REINTERPRET (i, u)
#define U_AS_I MAKE_REINTERPRET (u, i)

#define VMVI_VX MAKE_FUNC (__riscv_vmv_v_x_i)
#define VFMV_VF MAKE_FUNC (__riscv_vfmv_v_f_f)

#define IDENTIFY_F IDENTIFY

/* The source-shaped acosf helper is embedded by the acospif comparator.
   Select its public body from the importer function name so the helper
   skeleton has no routine-local enable/disable macro.  */
#ifndef RVV_MATH_IS_ACOSF_HELPERS_ONLY_acospif_original
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY_acospif_original 1
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_ACOSF_HELPERS_ONLY_##name
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY(name) \
  RVV_MATH_IS_ACOSF_HELPERS_ONLY_INNER (name)
#endif

#ifndef RVV_LIBMVEC_TRIG_DOUBLE_HELPERS
#define RVV_LIBMVEC_TRIG_DOUBLE_HELPERS
#define PIBY2_HI 0x1.921fb54442d18p+0
#define PIBY2_MID 0x1.1a62633145c07p-54
#define RVV_PI_Q60 0x3243f6a8885a308dLL
#define RVV_PIBY2_Q60 0x1921fb54442d1847LL
#define RVV_LN2_HI_D 0x1.62e42fefa39efp-1
#define RVV_LN2_LO_D 0x1.abc9e3b39803fp-56
#define RVV_SLEEF_LOG_U1_LN2_HI_F 0.69314718246459960938f
#define RVV_SLEEF_LOG_U1_LN2_LO_F -1.904654323148236017e-09f
#define RVV_SLEEF_LOG_U1_POLY0_F 0.3027294874e+0f
#define RVV_SLEEF_LOG_U1_POLY1_F 0.3996108174e+0f
#define RVV_SLEEF_LOG_U1_POLY2_F 0.6666694880e+0f

#endif /* RVV_LIBMVEC_TRIG_DOUBLE_HELPERS */

/* Shared source-shaped acosf names and constants.  The clean finite kernel
   and the independent SLEEF comparator use distinct prefixes, but both
   source operation graphs are instantiated by the same importer LMUL.  */
#ifndef RVV_ACOSF_COMMON_DEFINITIONS
# define RVV_ACOSF_COMMON_DEFINITIONS
# define RVV_ACOSF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ACOSF_PIBY2_HI_F 0x1.921fb6p+0f
# define RVV_ACOSF_PIBY2_LO_F -0x1.777a5cp-25f
# define RVV_ACOSF_PI_HI_F RVV_SLEEF_TRIG_PI_F
# define RVV_ACOSF_PI_LO_F -0x1.777a5cp-24f
# define RVV_ACOSF_C0 0.1666677296e+0f
# define RVV_ACOSF_C1 0.7495029271e-1f
# define RVV_ACOSF_C2 0.4547423869e-1f
# define RVV_ACOSF_C3 0.2424046025e-1f
# define RVV_ACOSF_C4 0.4197454825e-1f

/* The pinned SLEEF source passes the five polynomial literals as double
   expressions to vcast_vf_f(float).  Keep them unsuffixed here so the
   conversion remains observable under a non-RNE caller FRM.  */
# define RVV_SLEEF_ORIG_ACOSF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_SLEEF_ORIG_ACOSF_C0 0.1666677296e+0
# define RVV_SLEEF_ORIG_ACOSF_C1 0.7495029271e-1
# define RVV_SLEEF_ORIG_ACOSF_C2 0.4547423869e-1
# define RVV_SLEEF_ORIG_ACOSF_C3 0.2424046025e-1
# define RVV_SLEEF_ORIG_ACOSF_C4 0.4197454825e-1
# define RVV_SLEEF_ORIG_ACOSF_PI2_HI 3.1415927410125732422f / 2
# define RVV_SLEEF_ORIG_ACOSF_PI2_LO -8.7422776573475857731e-08f / 2
# define RVV_SLEEF_ORIG_ACOSF_PI_HI 3.1415927410125732422f
# define RVV_SLEEF_ORIG_ACOSF_PI_LO -8.7422776573475857731e-08f
#endif

/* Target-owned source-shaped SLEEF DF/FMA helpers.  These macros are
   kept in the common math header so skeletons do not carry a private
   helper header.  The operation graphs are unchanged from the former
   former private SLEEF helper translation unit.  */
#ifndef RVV_SLEEF_ORIGINAL_DF_HELPERS
#define RVV_SLEEF_ORIGINAL_DF_HELPERS

#define RVV_SLEEF_ORIG_MANGLE_INNER(base, lmul) base##lmul
#define RVV_SLEEF_ORIG_MANGLE(base, lmul) \
  RVV_SLEEF_ORIG_MANGLE_INNER (base, lmul)

#define RVV_SLEEF_ORIG_DFADD_VF_VF(x, y, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x, _y, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_x, _s, (vlen)), _y, (vlen));              \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)),                   \
          VFSUB_VV (_y_hi, _v, (vlen)), (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)      \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x_hi, _y, (vlen));                               \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x_hi, _s, (vlen)), _y, (vlen));        \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _x_lo, (vlen));                                  \
    }                                                                         \
  while (0)

/* Source dfsub_vf2_vf2_vf2: pair subtraction with the source's sequential
   high-part subtraction and low-part correction order.  */
#define RVV_SLEEF_ORIG_DFSUB_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFSUB_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _t = VFSUB_VV (_x_hi, _s, (vlen));                               \
      _t = VFSUB_VV (_t, _y_hi, (vlen));                                      \
      _t = VFADD_VV (_t, _x_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFSUB_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFSCALE_VF2_VF(x_hi, x_lo, scale, r_hi, r_lo, vlen)   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _scale = (scale);                                                \
      (r_hi) = VFMUL_VV (_x_hi, _scale, (vlen));                              \
      (r_lo) = VFMUL_VV (_x_lo, _scale, (vlen));                              \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFMUL_VF_VF(x, y, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFMUL_VV (_x, _y, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMSUB_VV (_x, _y, _s, (vlen));                                \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFMUL_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _residual = VFMSUB_VV (_x_hi, _y_hi, _s, (vlen));                \
      VFLOAT _cross_hi = VFMADD_VV (_x_lo, _y_hi, _residual, (vlen));         \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_x_hi, _y_lo, _cross_hi, (vlen));                   \
    }                                                                         \
  while (0)

#endif /* RVV_SLEEF_ORIGINAL_DF_HELPERS */
