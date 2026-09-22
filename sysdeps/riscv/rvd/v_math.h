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
static const uint64_t class_Inf = 0x81;
static const uint64_t class_Zero = 0x18;
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
#undef VMVU_VX
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
#define VMVU_VX MAKE_FUNC (__riscv_vmv_v_x_u)
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

#define PIBY2_INV 0x1.45f306dc9c883p-1
#define PIBY2_HI 0x1.921fb54442d18p+0
#define PIBY2_MID 0x1.1a62633145c07p-54
#define PIBY2_LO -0x1.f1976b7ed8fbcp-110
#define RVV_PI_Q60 0x3243f6a8885a308dLL
#define RVV_PIBY2_Q60 0x1921fb54442d1847LL
#define RVV_ONE_OVER_PI_HI_D 0x1.45f306dc9c883p-2
#define RVV_ONE_OVER_PI_Q63 0x28be60db9391054aLL
#define RVV_LN2_HI_D 0x1.62e42fefa39efp-1
#define RVV_LN2_LO_D 0x1.abc9e3b39803fp-56
#define RVV_LOG2_HI_D 0x1.62e42fefa4000p-1
#define RVV_LOG2_LO_D -0x1.8432a1b0e2634p-43
#define RVV_COSPI_POLY_RIGHT_SEED_A -0x1.8f5dd75850673p-37
#define RVV_COSPI_POLY_RIGHT_SEED_B 0x1.5d8b5ae12066ap-33
#define RVV_COSPI_POLY_RIGHT_OUTER_A -0x1.27e4f72551e3dp-22
#define RVV_COSPI_POLY_RIGHT_OUTER_B 0x1.71de35553ddb6p-19
#define RVV_COSPI_POLY_RIGHT_INNER_A 0x1.1ee950032f74cp-29
#define RVV_COSPI_POLY_RIGHT_INNER_B -0x1.ae5e4b94836f8p-26
#define RVV_COSPI_POLY_LEFT_SEED_A 0x1.a01a019b77545p-16
#define RVV_COSPI_POLY_LEFT_SEED_B -0x1.a01a019be932ap-13
#define RVV_COSPI_POLY_LEFT_OUTER_A 0x1.5555555555546p-5
#define RVV_COSPI_POLY_LEFT_OUTER_B -0x1.5555555555548p-3
#define RVV_COSPI_POLY_LEFT_INNER_A -0x1.6c16c16c1450cp-10
#define RVV_COSPI_POLY_LEFT_INNER_B 0x1.111111110f730p-7
#define RVV_SLEEF_LOG_U1_LN2_HI_F 0.69314718246459960938f
#define RVV_SLEEF_LOG_U1_LN2_LO_F -1.904654323148236017e-09f
#define RVV_SLEEF_LOG_U1_POLY0_F 0.3027294874e+0f
#define RVV_SLEEF_LOG_U1_POLY1_F 0.3996108174e+0f
#define RVV_SLEEF_LOG_U1_POLY2_F 0.6666694880e+0f

#define EXCEPTION_HANDLING_TRIG(vx, expo_x, special_args, vy_special, vlen)  \
  do                                                                          \
    {                                                                         \
      VUINT vclass = __riscv_vfclass ((vx), (vlen));                          \
      VBOOL nan_inf;                                                          \
      IDENTIFY (vclass, class_NaN | class_Inf, nan_inf, (vlen));             \
      VBOOL small_x                                                           \
        = __riscv_vmsleu ((expo_x), EXP_BIAS - MAN_LEN - 5, (vlen));         \
      (special_args) = __riscv_vmor (nan_inf, small_x, (vlen));               \
      if (__riscv_vcpop ((special_args), (vlen)) > 0)                         \
        {                                                                     \
          VBOOL id_mask;                                                      \
          IDENTIFY (vclass, class_Inf, id_mask, (vlen));                      \
          (vy_special)                                                        \
            = __riscv_vfmadd (small_x, (vx), 0x1.0p-60, (vx), (vlen));        \
          VFLOAT vy_nan_inf                                                   \
            = __riscv_vfmerge ((vx), fp_sNaN, id_mask, (vlen));               \
          vy_nan_inf                                                          \
            = __riscv_vfadd (nan_inf, vy_nan_inf, vy_nan_inf, (vlen));        \
          (vy_special)                                                        \
            = __riscv_vmerge ((vy_special), vy_nan_inf, nan_inf, (vlen));     \
          (vx) = __riscv_vfmerge ((vx), fp_posZero, (special_args), (vlen));  \
        }                                                                     \
    }                                                                         \
  while (0)

/* The 2^500 * (2/pi) table is owned by v_d_cos_data.c.  */

/* Trigonometric argument reduction for |x| >= 2^24.  */
#define LARGE_ARGUMENT_REDUCTION_Piby2(vx, vlen, x_large, n_xlarge, r_xlarge, \
                                       r_delta_xlarge)                        \
  do                                                                          \
    {                                                                         \
      VBOOL x_not_large = __riscv_vmnot ((x_large), vlen);                    \
      VFLOAT vx_work = __riscv_vfmerge ((vx), 0x1.0p30, x_not_large, vlen);   \
      VINT lsb_x = U_AS_I (__riscv_vsrl (F_AS_U (vx_work), 52, vlen));        \
      lsb_x = __riscv_vand (lsb_x, 0x7ff, vlen);                              \
      lsb_x = __riscv_vsub (lsb_x, 1075, vlen);                               \
      VUINT expo_mask = VMVU_VX (0x7ff, vlen);                                \
      expo_mask = __riscv_vsll (expo_mask, 52, vlen);                         \
      VFLOAT scale = U_AS_F (__riscv_vand (F_AS_U (vx_work), expo_mask,       \
                                           vlen));                            \
      scale = __riscv_vfmul (scale, 0x1.0p-500, vlen);                        \
      expo_mask = __riscv_vnot (expo_mask, vlen);                             \
      VUINT expo_1 = VMVU_VX (0x3ff, vlen);                                   \
      expo_1 = __riscv_vsll (expo_1, 52, vlen);                               \
      vx_work = U_AS_F (__riscv_vand (F_AS_U (vx_work), expo_mask, vlen));    \
      vx_work = U_AS_F (__riscv_vor (F_AS_U (vx_work), expo_1, vlen));        \
      VUINT j_start = I_AS_U (__riscv_vmax (lsb_x, 2, vlen));                 \
      j_start = __riscv_vsub (j_start, 2, vlen);                              \
      j_start = __riscv_vmul (j_start, 20165, vlen);                          \
      j_start = __riscv_vsrl (j_start, 20, vlen);                             \
      VUINT ind = __riscv_vsll (j_start, 3, vlen);                            \
      VFLOAT two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind,    \
                                           vlen);                             \
      VUINT peg_expo = VMVU_VX (1524, vlen);                                  \
      peg_expo = __riscv_vnmsac (peg_expo, 52, j_start, vlen);                \
      VFLOAT peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      VFLOAT s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);              \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      VFLOAT ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);               \
      VFLOAT prod_0 = s;                                                      \
      VFLOAT prod_1 = ss;                                                     \
      prod_0 = __riscv_vfmul (prod_0, scale, vlen);                           \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_1 = __riscv_vfadd (prod_1, s, vlen);                               \
      VFLOAT prod_2 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_1 = __riscv_vfmul (prod_1, scale, vlen);                           \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_2 = __riscv_vfadd (prod_2, s, vlen);                               \
      VFLOAT prod_3 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_2 = __riscv_vfmul (prod_2, scale, vlen);                           \
      VBOOL ignore_prod_0 = __riscv_vmfge (scale, 0x1.0p-446, vlen);          \
      prod_0 = __riscv_vfmerge (prod_0, 0.0, ignore_prod_0, vlen);            \
      VFLOAT flt_n = __riscv_vfmul (prod_0, 0x1.0p-12, vlen);                 \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_0 = __riscv_vfnmsac (prod_0, 0x1.0p12, flt_n, vlen);               \
      flt_n = __riscv_vfmul (prod_1, 0x1.0p-12, vlen);                        \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_1 = __riscv_vfnmsac (prod_1, 0x1.0p12, flt_n, vlen);               \
      flt_n = __riscv_vfadd (prod_1, prod_2, vlen);                           \
      flt_n = __riscv_vfadd (flt_n, prod_0, vlen);                            \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_0 = __riscv_vfsub (prod_0, flt_n, vlen);                           \
      VFLOAT r_hi_local = __riscv_vfadd (prod_0, prod_1, vlen);               \
      VFLOAT r_lo_local = __riscv_vfsub (prod_0, r_hi_local, vlen);           \
      r_lo_local = __riscv_vfadd (r_lo_local, prod_1, vlen);                  \
      VFLOAT tmp_1 = __riscv_vfadd (r_hi_local, prod_2, vlen);                \
      VFLOAT tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                 \
      tmp_2 = __riscv_vfadd (tmp_2, prod_2, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_3 = __riscv_vfadd (prod_3, s, vlen);                               \
      VFLOAT prod_4 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_3 = __riscv_vfmul (prod_3, scale, vlen);                           \
      tmp_1 = __riscv_vfadd (r_hi_local, prod_3, vlen);                       \
      tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                        \
      tmp_2 = __riscv_vfadd (tmp_2, prod_3, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      prod_4 = __riscv_vfadd (prod_4, s, vlen);                               \
      prod_4 = __riscv_vfmul (prod_4, scale, vlen);                           \
      tmp_1 = __riscv_vfadd (r_hi_local, prod_4, vlen);                       \
      tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                        \
      tmp_2 = __riscv_vfadd (tmp_2, prod_4, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      (r_xlarge) = __riscv_vfmul (r_hi_local, PIBY2_HI, vlen);                \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmsub (r_hi_local, PIBY2_HI, (r_xlarge), vlen);            \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmacc ((r_delta_xlarge), PIBY2_MID, r_hi_local, vlen);     \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmacc ((r_delta_xlarge), PIBY2_HI, r_lo_local, vlen);      \
    }                                                                         \
  while (0)

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
