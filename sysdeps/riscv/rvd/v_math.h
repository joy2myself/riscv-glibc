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

static const uint64_t class_sNaN = 0x100;
static const uint64_t class_NaN = 0x300;
static const uint64_t class_negInf = 0x1;
static const uint64_t class_Inf = 0x81;
static const uint64_t class_negZero = 0x8;
static const uint64_t class_Zero = 0x18;
static const uint64_t class_posZero = 0x10;
static const uint64_t class_posDenorm = 0x20;
static const uint64_t class_Subnormal = 0x24;
static const uint64_t class_negative = 0x7;
static const uint64_t class_finite_neg = 0x6;
static const uint64_t class_finite_pos = 0x60;

/* Compile-time constants shared by the source-shaped single-precision
   expm1-family graphs.  These values are scalar constants, rather than
   indexed data, so keep them in the math header instead of a data TU.  */
#define RVV_SLEEF_EXPM1F_R_LN2 0x1.715476p+0f
#define RVV_SLEEF_EXPM1F_L2U 0x1.62e400p-1f
#define RVV_SLEEF_EXPM1F_L2L 0x1.7f7d1cp-20f
#define RVV_SLEEF_EXPM1F_UNDERFLOW -104.0f
#define RVV_SLEEF_EXPM1F_POLY0 0.1980960224e-3f
#define RVV_SLEEF_EXPM1F_POLY1 0.1394256484e-2f
#define RVV_SLEEF_EXPM1F_POLY2 0.8333456703e-2f
#define RVV_SLEEF_EXPM1F_POLY3 0.4166637361e-1f
#define RVV_SLEEF_EXPM1F_TAYLOR_C \
  0.166666659414234244790680580464e+0f
#define RVV_SLEEF_EXPM1F_MAX_FINITE 0x1.fffffep+127f
#define RVV_SLEEF_EXPM1F_HALF 0.5f
#define RVV_SLEEF_EXPM1F_TWO 2.0f

/* The expk2f large-input coshf path uses the source split multiplied by two.
   Keep the already-rounded binary32 products as constants.  */
#define RVV_SLEEF_COSHF_TWO_L2U 0x1.62e400p+0f
#define RVV_SLEEF_COSHF_TWO_L2L 0x1.7f7d1cp-19f

#define RVV_LIBMVEC_COSHF_LARGE_INPUT 88.0f
#define RVV_LIBMVEC_COSHF_TINY_INPUT 0x1p-48f
#define RVV_LIBMVEC_COSHF_OVERFLOW 89.415985f
#define RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_HI RVV_SLEEF_COSHF_TWO_L2U
#define RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_LO RVV_SLEEF_COSHF_TWO_L2L
#define RVV_SLEEF_COSHF_SOURCE_THRESHOLD 89.0f

/* Scalar constants shared by the SLEEF single-precision trigonometric
   reduction paths.  Keep the actual rempif lookup table in
   v_f_sleef_trig_data.c; these values are compile-time constants, not
   indexed data.  */
#define RVV_SLEEF_TRIG_ONE_OVER_PI 0x1.45f306p-2f
#define RVV_SLEEF_TRIG_PI_F 0x1.921fb6p+1f
#define RVV_SLEEF_TRIG_RANGE_MAX2 125.0f
#define RVV_SLEEF_TRIG_REMPI_TWO23 8388608.0f
#define RVV_SLEEF_TRIG_REMPI_FOUR 4.0f
#define RVV_SLEEF_TRIG_REMPI_NEG_QUARTER -0.25f
#define RVV_SLEEF_TRIG_REMPI_NEG_FOUR -4.0f
#define RVV_SLEEF_TRIG_REMPI_TWO_PI_HI \
  (3.1415927410125732422f * 2.0f)
#define RVV_SLEEF_TRIG_REMPI_TWO_PI_LO \
  (-8.7422776573475857731e-08f * 2.0f)
#define RVV_SLEEF_TRIG_REMPI_SMALL 0.7f

#define RVV_SLEEF_COSF_NEG_PI_A2_HALF \
  (-3.1414794921875f * 0.5f)
#define RVV_SLEEF_COSF_NEG_PI_B2_HALF \
  (-0.00011315941810607910156f * 0.5f)
#define RVV_SLEEF_COSF_NEG_PI_C2_HALF \
  (-1.9841872589410058936e-09f * 0.5f)
#define RVV_SLEEF_COSF_REMPI_HALF_PI_HI \
  (-3.1415927410125732422f * 0.5f)
#define RVV_SLEEF_COSF_REMPI_HALF_PI_LO \
  (-8.7422776573475857731e-08f * -0.5f)
#define RVV_SLEEF_COSF_REDUCE_NEG_HALF -0.5f
#define RVV_SLEEF_COSF_REDUCE_TWO 2.0f
#define RVV_SLEEF_COSF_POLY_0 2.6083159809786593541503e-06f
#define RVV_SLEEF_COSF_POLY_1 -0.0001981069071916863322258f
#define RVV_SLEEF_COSF_POLY_2 0.00833307858556509017944336f
#define RVV_SLEEF_COSF_POLY_3 -0.166666597127914428710938f

/* Scalar constants for the SLEEF single-precision cbrtf graph.  */
#define RVV_SLEEF_CBRTF_VILOGBK_THRESH 5.421010862427522e-20f
#define RVV_SLEEF_CBRTF_VILOGBK_SCALE 1.8446744073709552e19f
#define RVV_SLEEF_CBRTF_VILOGBK_MASK 0xffu
#define RVV_SLEEF_CBRTF_VILOGBK_SMALL_BIAS 0x7f
#define RVV_SLEEF_CBRTF_VILOGBK_LARGE_BIAS (64 + 0x7f)
#define RVV_SLEEF_CBRTF_Q2_1_HI 1.2599210739135742188f
#define RVV_SLEEF_CBRTF_Q2_1_LO -2.4018701694217270415e-08f
#define RVV_SLEEF_CBRTF_Q2_2_HI 1.5874010324478149414f
#define RVV_SLEEF_CBRTF_Q2_2_LO 1.9520385308169352356e-08f
#define RVV_SLEEF_CBRTF_ONE_THIRD (1.0f / 3.0f)
#define RVV_SLEEF_CBRTF_NEG_ONE_THIRD (-1.0f / 3.0f)
#define RVV_SLEEF_CBRTF_NEG_TWO_THIRDS (-2.0f / 3.0f)
#define RVV_SLEEF_CBRTF_QUOTIENT_OFFSET_F 6144.0f
#define RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR_F 3.0f
#define RVV_SLEEF_CBRTF_EXPONENT_OFFSET 1
#define RVV_SLEEF_CBRTF_REMAINDER_ONE 1
#define RVV_SLEEF_CBRTF_REMAINDER_TWO 2
#define RVV_SLEEF_CBRTF_QUOTIENT_OFFSET 6144
#define RVV_SLEEF_CBRTF_QUOTIENT_MUL 2731
#define RVV_SLEEF_CBRTF_QUOTIENT_SHIFT 13
#define RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR 3
#define RVV_SLEEF_CBRTF_EXP_ADJUST 2048
#define RVV_SLEEF_CBRTF_POLY0_F -0.601564466953277587890625f
#define RVV_SLEEF_CBRTF_POLY1_F 2.8208892345428466796875f
#define RVV_SLEEF_CBRTF_POLY2_F -5.532182216644287109375f
#define RVV_SLEEF_CBRTF_POLY3_F 5.898262500762939453125f
#define RVV_SLEEF_CBRTF_POLY4_F -3.8095417022705078125f
#define RVV_SLEEF_CBRTF_POLY5_F 2.2241256237030029296875f

/* SLEEF's vcast_vf_f converts unsuffixed double expressions to float at the
   source operation point.  Keep that conversion tied to the surrounding
   source-shaped graph: GCC otherwise may hoist a data-backed conversion
   across the explicit-RNE vector-conversion window.  */
static inline float
RVV_SLEEF_ORIG_F64_TO_F32 (double value)
{
  float result;
  __asm__ __volatile__ ("fcvt.s.d %0, %1"
                        : "=f" (result)
                        : "f" (value)
                        : "memory");
  return result;
}

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
#undef F_AS_I_F
#undef F_AS_U
#undef F_AS_U_F
#undef I_AS_F
#undef I_AS_F_F
#undef U_AS_F
#undef U_AS_F_F
#undef I_AS_U
#undef I_AS_U_F
#undef U_AS_I
#undef U_AS_I_F
#undef VMVI_VX
#undef VMVI_VX_F
#undef VMVU_VX
#undef VMVU_VX_F
#undef VFMV_VF
#undef VFMV_VF_F
#undef EXP_BIAS
#undef MAN_LEN
#undef EXP_BIAS_F
#undef MAN_LEN_F
#undef BIT_WIDTH_F
#undef fp_negOne_f
#undef fp_sNaN_f

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
static const scalar_t fp_negInf = -__builtin_inf ();
static const scalar_t fp_negZero = -0.;
static const scalar_t fp_posZero = 0.0;
static const scalar_t fp_posOne = 0x1.0p0;
static const scalar_t fp_negOne = -0x1.0p0;
static const scalar_t fp_posHalf = 0x1.0p-1;
static const scalar_t fp_posInf = __builtin_inf ();
# endif
#elif (BIT_WIDTH == 32)
# define EXP_BIAS 127
# define MAN_LEN 23
# ifndef _V_MATH_SCALAR_32
#  define _V_MATH_SCALAR_32
typedef float scalar_t;
static const scalar_t fp_sNaN = __builtin_nansf ("");
static const scalar_t fp_negInf = -__builtin_inff ();
static const scalar_t fp_negZero = -0.f;
static const scalar_t fp_posZero = 0.0f;
static const scalar_t fp_posOne = 0x1.0p0f;
static const scalar_t fp_negOne = -0x1.0p0f;
static const scalar_t fp_posHalf = 0x1.0p-1f;
static const scalar_t fp_posInf = __builtin_inff ();
# endif
#else
# error "requested BIT_WIDTH unsupported"
#endif

#define EXP_BIAS_F EXP_BIAS
#define MAN_LEN_F MAN_LEN
#define BIT_WIDTH_F BIT_WIDTH
#define fp_posOne_f fp_posOne
#define fp_posZero_f fp_posZero
#define fp_posInf_f fp_posInf
#define fp_negOne_f fp_negOne
#define fp_sNaN_f fp_sNaN

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

/* Float-qualified spellings used by source-shaped SLEEF adapters.  The
   importer already selects the element width, so these are aliases rather
   than a second intrinsic family.  */
#define F_AS_I_F F_AS_I
#define F_AS_U_F F_AS_U
#define I_AS_F_F I_AS_F
#define U_AS_F_F U_AS_F
#define I_AS_U_F I_AS_U
#define U_AS_I_F U_AS_I
#define VMVI_VX_F VMVI_VX
#define VMVU_VX_F VMVU_VX
#define VFMV_VF_F VFMV_VF

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
#define RVV_PIBY2_Q61 0x3243f6a8885a308dLL
#define RVV_ONE_OVER_PI_HI_D 0x1.45f306dc9c883p-2
#define RVV_ONE_OVER_PI_LO_D -0x1.6b01ec5417056p-56
#define RVV_ONE_OVER_PI_Q63 0x28be60db9391054aLL
#define RVV_ATAN2_PI_075_HI_D 0x1.2d97c7f3321d2p+1
#define RVV_ATAN2_PIBY4_HI_D 0x1.921fb54442d18p-1
#define RVV_ATAN2_PI_HI_D 0x1.921fb54442d18p+1
#define RVV_ATAN2_PI_LO_D 0x1.1a62633145c07p-53
#define RVV_ATAN2_PI_Q61 0x6487ed5110b4611aLL
#define RVV_ATAN2_HALF_Q61 0x1000000000000000LL
#define RVV_ATAN2_ONE_Q61 0x2000000000000000LL
#define RVV_ATAN2_ONE_OVER_4_D 0x1.0p-2
#define RVV_ATAN2_THREE_OVER_4_D 0x1.8p-1
#define RVV_LN2_INV_D 0x1.71547652b82fep+0
#define RVV_LN2_HI_D 0x1.62e42fefa39efp-1
#define RVV_LN2_LO_D 0x1.abc9e3b39803fp-56
#define RVV_LOG2_HI_D 0x1.62e42fefa4000p-1
#define RVV_LOG2_LO_D -0x1.8432a1b0e2634p-43
#define RVV_LOG2_BY2_HI_D 0x1.62e42fefa4000p-2
#define RVV_LOG2_BY2_LO_D -0x1.8432a1b0e2634p-44
#define RVV_ONE_Q60 0x1000000000000000LL
#define RVV_ONE_OVER_3_D 0x1.5555555555555p-2
#define RVV_THIRD_Q62 0x1555555555555555LL
#define RVV_ONE_Q62 0x4000000000000000LL
#define RVV_CBRT_2_Q62 0x50a28be635ca2b89LL
#define RVV_CBRT_4_Q62 0x6597fa94f5b8f20bLL
#define RVV_COSPI_R_PRIME_SCALE -0x1.0p-1
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
#define RVV_COSPI_SMALL_EXP_OFFSET 5
#define RVV_COSPI_LARGE_EXP_OFFSET 53
#define RVV_COSHF_X_MAX_D 0x1.65p9
#define RVV_COSHF_SMALL_EXP_OFFSET 5
#define RVV_COSHF_SCALE_LIMIT 60
#define RVV_COSHF_LARGE_N 50
#define RVV_COSH_P_EVEN0_D 0x1.555555555555ap-3
#define RVV_COSH_P_EVEN1_D 0x1.111111110ef6ap-7
#define RVV_COSH_P_EVEN2_D 0x1.a01a01b32b633p-13
#define RVV_COSH_P_EVEN3_D 0x1.71ddef82f4beep-19
#define RVV_COSH_P_EVEN4_D 0x1.af6eacd796f0bp-26
#define RVV_COSH_P_ODD0_D 0x1.5555555553aefp-5
#define RVV_COSH_P_ODD1_D 0x1.6c16c17a09506p-10
#define RVV_COSH_P_ODD2_D 0x1.a019b37a2b3dfp-16
#define RVV_COSH_P_ODD3_D 0x1.289788d8bdadfp-22
#define RVV_EXP2_OVERFLOW_BOUND_D 0x1.0p10
#define RVV_EXP2_UNDERFLOW_BOUND_D -0x1.f4p10
#define RVV_EXP2_POLY_HI0_D 0.7073164598085707425e-8
#define RVV_EXP2_POLY_HI1_D 0.4434359082926529454e-9
#define RVV_EXP2_POLY_MID_HI_A0_D 0.1321543872511327615e-5
#define RVV_EXP2_POLY_MID_HI_A1_D 0.1017819260921760451e-6
#define RVV_EXP2_POLY_MID_HI_B0_D 0.1540353045101147808e-3
#define RVV_EXP2_POLY_MID_HI_B1_D 0.1525273353517584730e-4
#define RVV_EXP2_POLY_MID_LO_A0_D 0.9618129107597600536e-2
#define RVV_EXP2_POLY_MID_LO_A1_D 0.1333355814670499073e-2
#define RVV_EXP2_POLY_MID_LO_B0_D 0.2402265069591012214e+0
#define RVV_EXP2_POLY_MID_LO_B1_D 0.5550410866482046596e-1
#define RVV_EXP2F_CLEAN_MAX_INPUT_F 128.0f
#define RVV_EXP2F_CLEAN_MIN_INPUT_F -150.0f
#define RVV_EXP2F_CLEAN_MAX_FINITE_F 0x1.fffffep127f
#define RVV_EXP2F_CLEAN_MIN_NORMAL_F 0x1.0p-126f
#define RVV_EXP2F_SOURCE_MAX_INPUT_F 128.0f
#define RVV_EXP2F_SOURCE_MIN_INPUT_F -150.0f
#define RVV_EXP2F_CLEAN_POLY0_F 0.1535920892e-3f
#define RVV_EXP2F_CLEAN_POLY1_F 0.1339262701e-2f
#define RVV_EXP2F_CLEAN_POLY2_F 0.9618384764e-2f
#define RVV_EXP2F_CLEAN_POLY3_F 0.5550347269e-1f
#define RVV_EXP2F_CLEAN_POLY4_F 0.2402264476e+0f
#define RVV_EXP2F_CLEAN_POLY5_F 0.6931471825e+0f
#define RVV_EXP2F_SOURCE_POLY0_D +0.1535920892e-3
#define RVV_EXP2F_SOURCE_POLY1_D +0.1339262701e-2
#define RVV_EXP2F_SOURCE_POLY2_D +0.9618384764e-2
#define RVV_EXP2F_SOURCE_POLY3_D +0.5550347269e-1
#define RVV_EXP2F_SOURCE_POLY4_D +0.2402264476e+0
#define RVV_EXP2F_SOURCE_POLY5_D +0.6931471825e+0
#define RVV_SLEEF_LOG_U1_LN2_HI_F 0.69314718246459960938f
#define RVV_SLEEF_LOG_U1_LN2_LO_F -1.904654323148236017e-09f
#define RVV_SLEEF_LOG_U1_POLY0_F 0.3027294874e+0f
#define RVV_SLEEF_LOG_U1_POLY1_F 0.3996108174e+0f
#define RVV_SLEEF_LOG_U1_POLY2_F 0.6666694880e+0f
#define RVV_LOG1PF_LN2_HI_F RVV_SLEEF_LOG_U1_LN2_HI_F
#define RVV_LOG1PF_LN2_LO_F RVV_SLEEF_LOG_U1_LN2_LO_F
#define RVV_COSPIF_TINY_EXPONENT 99u
#define RVV_COSPIF_LARGE_EXPONENT 151u
#define RVV_COSPIF_TRIG_RANGE_MAX4_F 8.0e6f
#define RVV_COSPIF_POLY_COS0_F 0.3093842054e-6f
#define RVV_COSPIF_POLY_SIN0_F -0.2430611801e-7f
#define RVV_COSPIF_POLY_COS1_F -0.3657307388e-4f
#define RVV_COSPIF_POLY_SIN1_F 0.3590577080e-5f
#define RVV_COSPIF_POLY_COS2_F 0.2490393585e-2f
#define RVV_COSPIF_POLY_SIN2_F -0.3259917721e-3f
#define RVV_COSPIF_DF_FIRST_COS_HI_F -0.080745510756969451904f
#define RVV_COSPIF_DF_FIRST_COS_LO_F -1.3373665339076936258e-09f
#define RVV_COSPIF_DF_FIRST_SIN_HI_F 0.015854343771934509277f
#define RVV_COSPIF_DF_FIRST_SIN_LO_F 4.4940051354032242811e-10f
#define RVV_COSPIF_DF_SECOND_COS_HI_F 0.78539818525314331055f
#define RVV_COSPIF_DF_SECOND_COS_LO_F -2.1857338617566484855e-08f
#define RVV_COSPIF_DF_SECOND_SIN_HI_F -0.30842512845993041992f
#define RVV_COSPIF_DF_SECOND_SIN_LO_F -9.0728339030733922277e-09f

/* Fixed coefficients for the single-precision gamma helpers.  Preserve
   the source float/double types and the explicit conversion points.  */
/* Shared SLEEF single-precision logk2f polynomial.  */
#define RVV_SLEEF_LOGK2F_POLY0_F 0.2392828464508056640625f
#define RVV_SLEEF_LOGK2F_POLY1_F 0.28518211841583251953125f
#define RVV_SLEEF_LOGK2F_POLY2_F 0.400005877017974853515625f
#define RVV_SLEEF_LOGK2F_POLY3_F 0.666666686534881591796875f

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

/* The source-shaped asinf helper is embedded by the asinpif comparator.
   Select its public body from the importer function name so the helper
   skeleton has no routine-local enable/disable macro.  */
#ifndef RVV_MATH_IS_ASINF_HELPERS_ONLY_asinpif_original
# define RVV_MATH_IS_ASINF_HELPERS_ONLY_asinpif_original 1
# define RVV_MATH_IS_ASINF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_ASINF_HELPERS_ONLY_##name
# define RVV_MATH_IS_ASINF_HELPERS_ONLY(name) \
  RVV_MATH_IS_ASINF_HELPERS_ONLY_INNER (name)
#endif

/* The source-shaped cosf helper is embedded by the fastcosf comparator and
   by the clean fastcosf implementation.  Select its public body from the
   importer function name so nested skeletons need no control defines.  */
#ifndef RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf
# define RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf 1
# define RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf_original 1
# define RVV_MATH_IS_COSF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_COSF_HELPERS_ONLY_##name
# define RVV_MATH_IS_COSF_HELPERS_ONLY(name) \
  RVV_MATH_IS_COSF_HELPERS_ONLY_INNER (name)
#endif

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

/* Shared source-shaped asinf names and constants.  The clean importer and
   the independent SLEEF comparator use separate prefixes for their helper
   symbols, while both mangle through the current importer LMUL.  The
   `_original` polynomial literals intentionally remain unsuffixed: pinned
   SLEEF passes them as double expressions to vcast_vf_f(float), so their
   conversion remains source-shaped under a non-RNE caller FRM.  */
#ifndef RVV_ASINF_COMMON_DEFINITIONS
# define RVV_ASINF_COMMON_DEFINITIONS
# define RVV_ASINF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ASINF_PIBY4_HI 3.1415927410125732422f / 4.0f
# define RVV_ASINF_PIBY4_LO -8.7422776573475857731e-08f / 4.0f
# define RVV_ASINF_C0 0.1666677296e+0f
# define RVV_ASINF_C1 0.7495029271e-1f
# define RVV_ASINF_C2 0.4547423869e-1f
# define RVV_ASINF_C3 0.2424046025e-1f
# define RVV_ASINF_C4 0.4197454825e-1f
# define RVV_SLEEF_ORIG_ASINF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_SLEEF_ORIG_ASINF_C0 0.1666677296e+0
# define RVV_SLEEF_ORIG_ASINF_C1 0.7495029271e-1
# define RVV_SLEEF_ORIG_ASINF_C2 0.4547423869e-1
# define RVV_SLEEF_ORIG_ASINF_C3 0.2424046025e-1
# define RVV_SLEEF_ORIG_ASINF_C4 0.4197454825e-1
# define RVV_SLEEF_ORIG_ASINF_PI4_HI 3.1415927410125732422f / 4.0f
# define RVV_SLEEF_ORIG_ASINF_PI4_LO \
  -8.7422776573475857731e-08f / 4.0f
#endif

/* Source-shaped helper names used by the asinpif composite.  Keep this
   importer-dependent mangle in the common math header; the routine skeleton
   itself should contain only the finite operation graph and data references.  */
#ifndef RVV_ASINPIF_COMMON_DEFINITIONS
# define RVV_ASINPIF_COMMON_DEFINITIONS
# define RVV_ASINPIF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ASINPIF_ONE_OVER_PI_HI_F RVV_SLEEF_TRIG_ONE_OVER_PI
# define RVV_ASINPIF_ONE_OVER_PI_LO_F 0x1.b93910p-27f
# define RVV_ASINPIF_PI4_HI_F RVV_ASINF_PIBY4_HI
# define RVV_ASINPIF_PI4_LO_F RVV_ASINF_PIBY4_LO
# define RVV_ASINPIF_HALF_F 0.5f
# define RVV_ASINPIF_TWO_F 2.0f
# define RVV_ASINPIF_TINY_SQUARE_LIMIT_F 0x1p-63f
# define RVV_ASINPIF_POLY0_F 0.1666677296e+0f
# define RVV_ASINPIF_POLY1_F 0.7495029271e-1f
# define RVV_ASINPIF_POLY2_F 0.4547423869e-1f
# define RVV_ASINPIF_POLY3_F 0.2424046025e-1f
# define RVV_ASINPIF_POLY4_F 0.4197454825e-1f
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

#define RVV_SLEEF_ORIG_DFNORMALIZE_VF2(a_hi, a_lo, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _a_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _a_lo, (vlen));       \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFNEG_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      (r_hi) = VFNEG_V (_a_hi, (vlen));                                      \
      (r_lo) = VFNEG_V (_a_lo, (vlen));                                      \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF_VF(x, y, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x, _y, (vlen));                                  \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                      \
          VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)),                   \
          VFSUB_VV (_y, _v, (vlen)), (vlen));                                 \
    }                                                                         \
  while (0)

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

#define RVV_SLEEF_ORIG_DFADD2_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x_hi, _y, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_y, _v, (vlen)), (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _x_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF2_VF_S(x_hi, x_lo, y, r_hi, r_lo, vlen)       \
  RVV_SLEEF_ORIG_DFADD2_VF2_VF ((x_hi), (x_lo), VFMV_VF ((y), (vlen)),       \
                                (r_hi), (r_lo), (vlen))

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

#define RVV_SLEEF_ORIG_DFADD_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x, _s, (vlen)), _y_hi, (vlen));        \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x_hi, _s, (vlen)), _y_hi, (vlen));     \
      _t = VFADD_VV (_t, _x_lo, (vlen));                                      \
      _t = VFADD_VV (_t, _y_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

/* Source dfadd2_vf2_vf2_vf2: compensated pair-plus-pair addition.  */
#define RVV_SLEEF_ORIG_DFADD2_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,     \
                                      vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _v = VFSUB_VV (_s, _x_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_y_hi, _v, (vlen)), (vlen));                              \
      VFLOAT _lo = VFADD_VV (_x_lo, _y_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _lo, (vlen));                                    \
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

#define RVV_SLEEF_ORIG_DFMUL_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFMUL_VV (_x_hi, _y, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_x_lo, _y, VFMSUB_VV (_x_hi, _y, _s, (vlen)),       \
                          (vlen));                                           \
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

/* Source dfmul_vf_vf2_vf2: return the collapsed scalar product of two
   double-float pairs with the source's nested FMA association.  */
#define RVV_SLEEF_ORIG_DFMUL_VF2_VF2_SCALAR(x_hi, x_lo, y_hi, y_lo, result,   \
                                            vlen)                              \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _low = VFMADD_VV (_x_lo, _y_hi,                                 \
                               VFMUL_VV (_x_hi, _y_lo, (vlen)), (vlen));      \
      (result) = VFMADD_VV (_x_hi, _y_hi, _low, (vlen));                      \
    }                                                                         \
  while (0)

/* Source dfrec_vf2_vf in the FMA branch: reciprocal correction for a scalar
   high part retains the source's nested negative-product sequence.  */
#define RVV_SLEEF_ORIG_DFREC_VF2_VF(x, r_hi, r_lo, vlen)                     \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _x, (vlen));              \
      VFLOAT _u = VFNMSUB_VV (_x, _t, VFMV_VF (1.0f, (vlen)), (vlen));         \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_t, _u, (vlen));                                     \
    }                                                                         \
  while (0)

/* Source dfrec_vf2_vf2 in the FMA branch: reciprocal correction retains the
   source's nested negative-product sequence.  */
#define RVV_SLEEF_ORIG_DFREC_VF2_VF2(x_hi, x_lo, r_hi, r_lo, vlen)            \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _x_hi, (vlen));            \
      VFLOAT _u = VFNMSUB_VV (_x_hi, _t, VFMV_VF (1.0f, (vlen)), (vlen));      \
      _u = VFNMSAC_VV (_u, _x_lo, _t, (vlen));                                \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_t, _u, (vlen));                                     \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFDIV_VF2_VF2(n_hi, n_lo, d_hi, d_lo, r_hi, r_lo,     \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _n_hi = (n_hi);                                                  \
      VFLOAT _n_lo = (n_lo);                                                  \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _s = VFMUL_VV (_n_hi, _t, (vlen));                               \
      VFLOAT _u = VFMSUB_VV (_t, _n_hi, _s, (vlen));                          \
      VFLOAT _v = VFNMSUB_VV (_d_hi, _t, VFMV_VF (1.0f, (vlen)), (vlen));     \
      _v = VFNMSUB_VV (_d_lo, _t, _v, (vlen));                                \
      VFLOAT _nlo_t = VFMADD_VV (_n_lo, _t, _u, (vlen));                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_s, _v, _nlo_t, (vlen));                            \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFSQU_VF2(x_hi, x_lo, r_hi, r_lo, vlen)                \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_x_hi, _x_hi, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (VFADD_VV (_x_hi, _x_hi, (vlen)), _x_lo,             \
                          VFMSUB_VV (_x_hi, _x_hi, _s, (vlen)), (vlen));      \
    }                                                                         \
  while (0)

/* Source dfsqrt_vf2_vf2 in common/df.h.  The pinned RVV source does not
   define ENABLE_RECSQRT_SP, so this is the source's sqrt/reciprocal branch:
   sqrt(high+low), square it as a source DF pair, add it to the input pair,
   multiply by the source DF reciprocal, then scale by one half.  */
#define RVV_SLEEF_ORIG_DFSQRT_VF2_VF2(d_hi, d_lo, r_hi, r_lo, vlen)           \
  do                                                                          \
    {                                                                         \
      VFLOAT _d_hi = (d_hi);                                                   \
      VFLOAT _d_lo = (d_lo);                                                   \
      VFLOAT _t = __riscv_vfsqrt (VFADD_VV (_d_hi, _d_lo, (vlen)), (vlen));    \
      VFLOAT _tsq_hi, _tsq_lo;                                                \
      RVV_SLEEF_ORIG_DFMUL_VF_VF (_t, _t, _tsq_hi, _tsq_lo, (vlen));           \
      VFLOAT _sum_hi, _sum_lo;                                                \
      RVV_SLEEF_ORIG_DFADD2_VF2_VF2                                          \
        (_d_hi, _d_lo, _tsq_hi, _tsq_lo, _sum_hi, _sum_lo, (vlen));            \
      VFLOAT _rec_hi, _rec_lo;                                                \
      RVV_SLEEF_ORIG_DFREC_VF2_VF (_t, _rec_hi, _rec_lo, (vlen));              \
      VFLOAT _root_hi, _root_lo;                                              \
      RVV_SLEEF_ORIG_DFMUL_VF2_VF2                                           \
        (_sum_hi, _sum_lo, _rec_hi, _rec_lo, _root_hi, _root_lo, (vlen));      \
      RVV_SLEEF_ORIG_DFSCALE_VF2_VF                                           \
        (_root_hi, _root_lo, VFMV_VF (0.5f, (vlen)),                          \
         (r_hi), (r_lo), (vlen));                                             \
    }                                                                         \
  while (0)

#endif /* RVV_SLEEF_ORIGINAL_DF_HELPERS */

/* Shared asinhf constants and importer mangles.  Keep the clean and
   source-original operation graphs independently named while placing their
   public definitions in the common math header.  */
#ifndef RVV_ASINHF_COMMON_DEFINITIONS
# define RVV_ASINHF_COMMON_DEFINITIONS
# define RVV_ASINHF_SQRT_FLT_MAX 1.8446743e+19f
# define RVV_ASINHF_L2U 0.69314718246459960938f
# define RVV_ASINHF_L2L -1.904654323148236017e-09f
# define RVV_ASINHF_LOG_P0 RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_ASINHF_LOG_P1 RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_ASINHF_LOG_P2 RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_ASINHF_LOG_P3 RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_ASINHF_CLEAN_LOGK2F \
  __PASTE2 (rvv_libmvec_asinhf_logk2f_m, LMUL)
# define RVV_ASINHF_CLEAN_LOGF_LARGE \
  __PASTE2 (rvv_libmvec_asinhf_logf_large_m, LMUL)
# define RVV_SLEEF_ORIG_ASINHF_SQRT_FLT_MAX 1.8446743e+19f
# define RVV_SLEEF_ORIG_ASINHF_L2U 0.69314718246459960938f
# define RVV_SLEEF_ORIG_ASINHF_L2L -1.904654323148236017e-09f
# define RVV_SLEEF_ORIG_ASINHF_LOG_P0 RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P1 RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P2 RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P3 RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_SLEEF_ORIG_ASINHF_MANGLE(base, lmul) \
  RVV_SLEEF_ORIG_MANGLE (base, lmul)
#endif

/* Shared importer mangles for the single-precision atanh source-shaped
   helpers.  The operation bodies remain in the routine skeletons, while
   these LMUL-dependent names live in the common math header so the skeletons
   contain no routine-local helper macro block.  */
#ifndef RVV_ATANHF_COMMON_DEFINITIONS
# define RVV_ATANHF_COMMON_DEFINITIONS
# define RVV_ATANHF_LOG_POLY0_F RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_ATANHF_LOG_POLY1_F RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_ATANHF_LOG_POLY2_F RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_ATANHF_LOG_POLY3_F RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_LIBMVEC_ATANHF_VILOGBK \
  __PASTE2 (rvv_libmvec_atanhf_vilogbk_m, LMUL)
# define RVV_LIBMVEC_ATANHF_VPOW2I \
  __PASTE2 (rvv_libmvec_atanhf_vpow2i_m, LMUL)
# define RVV_LIBMVEC_ATANHF_LOGK2F \
  __PASTE2 (rvv_libmvec_atanhf_logk2f_m, LMUL)
#endif

/* The double-precision atan2 and atan2pi skeletons share one operation graph.
   Keep their mode-specific constants and fixed-point post-selection in this
   common header; the skeleton itself should not carry a compile-control
   define merely because the importer has a different public function name.  */
#ifndef RVV_MATH_IS_ATAN2_FAMILY_atan2
# define RVV_MATH_IS_ATAN2_FAMILY_atan2 1
# define RVV_MATH_IS_ATAN2_FAMILY_atan2pi 1
# define RVV_MATH_IS_ATAN2_FAMILY_INNER(name) \
  RVV_MATH_IS_ATAN2_FAMILY_##name
# define RVV_MATH_IS_ATAN2_FAMILY(name) \
  RVV_MATH_IS_ATAN2_FAMILY_INNER (name)
#endif

#ifndef RVV_MATH_IS_ATAN2PI_atan2pi
# define RVV_MATH_IS_ATAN2PI_atan2pi 1
# define RVV_MATH_IS_ATAN2PI_INNER(name) RVV_MATH_IS_ATAN2PI_##name
# define RVV_MATH_IS_ATAN2PI(name) RVV_MATH_IS_ATAN2PI_INNER (name)
#endif

#undef RVV_ATAN2_Y_NONZERO_X_ZERO
#undef RVV_ATAN2_Y_ZERO_X_NEGATIVE
#undef RVV_ATAN2_POSINF_POSINF
#undef RVV_ATAN2_NEGINF_INF
#undef RVV_ATAN2_NEGINF_FINITE
#undef RVV_ATAN2_FINITE_NEGATIVE_ZERO
#undef RVV_ATAN2_Y_INF_FINITE
#undef RVV_ATAN2_PI_HI
#undef RVV_ATAN2_PI_LO
#undef RVV_ATAN2_PIBY2_HI
#undef RVV_ATAN2_FIXED_HALF
#undef RVV_ATAN2_FIXED_ONE
#undef RVV_ATAN2_POST_SCALE

#if (BIT_WIDTH == 64) && RVV_MATH_IS_ATAN2_FAMILY (V_MATH_FUNC)
# if RVV_MATH_IS_ATAN2PI (V_MATH_FUNC)
#  define RVV_ATAN2_Y_NONZERO_X_ZERO (fp_posHalf)
#  define RVV_ATAN2_Y_ZERO_X_NEGATIVE (fp_posOne)
#  define RVV_ATAN2_POSINF_POSINF RVV_ATAN2_ONE_OVER_4_D
#  define RVV_ATAN2_NEGINF_INF RVV_ATAN2_THREE_OVER_4_D
#  define RVV_ATAN2_NEGINF_FINITE (fp_posOne)
#  define RVV_ATAN2_FINITE_NEGATIVE_ZERO (fp_posOne)
#  define RVV_ATAN2_Y_INF_FINITE (fp_posHalf)
#  define RVV_ATAN2_FIXED_HALF RVV_ATAN2_HALF_Q61
#  define RVV_ATAN2_FIXED_ONE RVV_ATAN2_ONE_Q61
#  define RVV_ATAN2_POST_SCALE(value, vlen) \
  __riscv_vsmul ((value), RVV_ONE_OVER_PI_Q63, 1, (vlen))
# else
#  define RVV_ATAN2_Y_NONZERO_X_ZERO PIBY2_HI
#  define RVV_ATAN2_Y_ZERO_X_NEGATIVE RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_POSINF_POSINF RVV_ATAN2_PIBY4_HI_D
#  define RVV_ATAN2_NEGINF_INF RVV_ATAN2_PI_075_HI_D
#  define RVV_ATAN2_NEGINF_FINITE RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_FINITE_NEGATIVE_ZERO RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_Y_INF_FINITE PIBY2_HI
#  define RVV_ATAN2_FIXED_HALF RVV_PIBY2_Q61
#  define RVV_ATAN2_FIXED_ONE RVV_ATAN2_PI_Q61
#  define RVV_ATAN2_POST_SCALE(value, vlen) (value)
# endif
# define RVV_ATAN2_PI_HI RVV_ATAN2_PI_HI_D
# define RVV_ATAN2_PI_LO RVV_ATAN2_PI_LO_D
# define RVV_ATAN2_PIBY2_HI PIBY2_HI
#endif

/* Clean asinhf is the first importer that needs this subset of the generic
   float-float operation family.  Keep only the operations used by asinhf
   here; later routines can extend the common family at their first use.  */
#ifndef RVV_ASINHF_CLEAN_DF_HELPERS
# define RVV_ASINHF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF_S(a_hi, a_lo, b, r_hi, r_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, VFMV_VF ((b), (vlen)), (vlen));            \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),      \
                     VFRSUB_VF (_v, (b), (vlen)), (vlen)),                    \
          _a_lo, (vlen));                                                     \
    }                                                                         \
  while (0)

# define DFADD_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _b_hi, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _b_hi, (vlen)),              \
          VFADD_VV (_a_lo, _b_lo, (vlen)), (vlen));                            \
    }                                                                         \
  while (0)

# define DFSQU_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      (r_hi) = VFMUL_VV (_a_hi, _a_hi, (vlen));                               \
      (r_lo) = VFMSUB_VV (_a_hi, _a_hi, (r_hi), (vlen));                      \
      (r_lo) = VFMACC_VV ((r_lo), VFADD_VV (_a_hi, _a_hi, (vlen)), _a_lo,    \
                          (vlen));                                           \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b, _s, (vlen));                          \
      _t = VFMACC_VV (_t, _a_lo, _b, (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_a_hi, _b_hi, (vlen));                            \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b_hi, _s, (vlen));                       \
      _t = VFMACC_VV (_t, _a_hi, _b_lo, (vlen));                              \
      _t = VFMACC_VV (_t, _a_lo, _b_hi, (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFDIV_VF2_VF2_VF2(n_hi, n_lo, d_hi, d_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _n_hi = (n_hi);                                                  \
      VFLOAT _n_lo = (n_lo);                                                  \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _s = VFMUL_VV (_n_hi, _t, (vlen));                               \
      VFLOAT _u = VFNMSUB_VV (_s, _d_hi, _n_hi, (vlen));                      \
      _u = VFNMSAC_VV (_u, _s, _d_lo, (vlen));                                \
      _u = VFADD_VV (_u, _n_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMUL_VV (_u, _t, (vlen));                                     \
    }                                                                         \
  while (0)

# define DFSCALE_VF2_VF2_S(a_hi, a_lo, s, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      (r_hi) = VFMUL_VF ((a_hi), (s), (vlen));                                \
      (r_lo) = VFMUL_VF ((a_lo), (s), (vlen));                                \
    }                                                                         \
  while (0)

# define DFSCALE_VF2_VF2_VF(a_hi, a_lo, s, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _s = (s);                                                        \
      (r_hi) = VFMUL_VV ((a_hi), _s, (vlen));                                 \
      (r_lo) = VFMUL_VV ((a_lo), _s, (vlen));                                 \
    }                                                                         \
  while (0)

#endif /* RVV_ASINHF_CLEAN_DF_HELPERS */

/* The atanf family uses these two generic float-float operations from the
   common helper family.  Keep them in v_math.h so the routine skeletons do
   not need a private helper header; the remaining generic operations can be
   added here by their first real user.  */
#ifndef RVV_VFLOAT2_COMMON_HELPERS
# define RVV_VFLOAT2_COMMON_HELPERS

# define VSEL_VF2_VO(o, a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)             \
  do                                                                          \
    {                                                                         \
      (r_hi) = __riscv_vmerge ((b_hi), (a_hi), (o), (vlen));                   \
      (r_lo) = __riscv_vmerge ((b_lo), (a_lo), (o), (vlen));                   \
    }                                                                         \
  while (0)

# define DFNEG_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                           \
  do                                                                          \
    {                                                                         \
      (r_hi) = VFNEG_V ((a_hi), (vlen));                                       \
      (r_lo) = VFNEG_V ((a_lo), (vlen));                                       \
    }                                                                         \
  while (0)

#endif /* RVV_VFLOAT2_COMMON_HELPERS */

/* Float-float operations used by the clean cosf kernel.  Keep these
   importer-dependent operations in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSF_CLEAN_DF_HELPERS
# define RVV_COSF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF_VF(a, b, r_hi, r_lo, vlen)                            \
  do                                                                          \
    {                                                                         \
      VFLOAT _a = (a);                                                        \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a, _b, (vlen));                                 \
      VFLOAT _v = VFSUB_VV (_s, _a, (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                      \
          VFSUB_VV (_a, VFSUB_VV (_s, _v, (vlen)), (vlen)),                  \
          VFSUB_VV (_b, _v, (vlen)), (vlen));                                \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a_hi, _b, (vlen));                              \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                              \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_b, _v, (vlen)), (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _a_lo, (vlen));                                 \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)       \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _b_hi, (vlen));                           \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                              \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_b_hi, _v, (vlen)), (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, VFADD_VV (_a_lo, _b_lo, (vlen)), (vlen));       \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF_VF(a, b, r_hi, r_lo, vlen)                              \
  do                                                                          \
    {                                                                         \
      VFLOAT _a = (a);                                                        \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a, _b, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMSUB_VV (_a, _b, _s, (vlen));                                \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b, _s, (vlen));                          \
      _t = VFMACC_VV (_t, _a_lo, _b, (vlen));                                \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFMUL_VF_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _y_hi = (y_hi);                                                 \
      VFLOAT _y_lo = (y_lo);                                                 \
      VFLOAT _t = VFMUL_VV (_x_hi, _y_lo, (vlen));                           \
      _t = VFMACC_VV (_t, _x_lo, _y_hi, (vlen));                             \
      (r) = VFMACC_VV (_t, _x_hi, _y_hi, (vlen));                            \
    }                                                                         \
  while (0)

# define DFNORMALIZE_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _a_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _a_lo, (vlen));       \
    }                                                                         \
  while (0)

#endif /* RVV_COSF_CLEAN_DF_HELPERS */

/* Float-float reciprocal used by the clean coshf kernel.  Keep this
   importer-dependent operation in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSHF_CLEAN_DF_HELPERS
# define RVV_COSHF_CLEAN_DF_HELPERS

# define DFREC_VF2_VF2(d_hi, d_lo, r_hi, r_lo, vlen)                         \
  do                                                                          \
    {                                                                         \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _u = VFNMSUB_VV (_t, _d_hi, VFMV_VF (1.0f, (vlen)), (vlen));      \
      _u = VFNMSAC_VV (_u, _t, _d_lo, (vlen));                                \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_u, _t, (vlen));                                     \
    }                                                                         \
  while (0)

#endif /* RVV_COSHF_CLEAN_DF_HELPERS */

/* Float-float additions used by the clean cospif kernel.  Keep these
   importer-dependent operations in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSPIF_CLEAN_DF_HELPERS
# define RVV_COSPIF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                      \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),      \
                    VFSUB_VV (_b, _v, (vlen)), (vlen)),                       \
          _a_lo, (vlen));                                                     \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)), \
                            VFSUB_VV (_y_hi, _v, (vlen)), (vlen));             \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#endif /* RVV_COSPIF_CLEAN_DF_HELPERS */

/* Fixed scalar coefficients and endpoint constants for the single-precision
   atan family.  */
#ifndef RVV_ATANF_CONSTANTS
# define RVV_ATANF_CONSTANTS
/* Polynomial shared by the single-precision atan kernels.  */
# define RVV_ATANF_POLY0_F -0.00176397908944636583328247f
# define RVV_ATANF_POLY1_F 0.0107900900766253471374512f
# define RVV_ATANF_POLY2_F -0.0309564601629972457885742f
# define RVV_ATANF_POLY3_F 0.0577365085482597351074219f
# define RVV_ATANF_POLY4_F -0.0838950723409652709960938f
# define RVV_ATANF_POLY5_F 0.109463557600975036621094f
# define RVV_ATANF_POLY6_F -0.142626821994781494140625f
# define RVV_ATANF_POLY7_F 0.199983194470405578613281f
# define RVV_ATANF_POLY8_F -0.333332866430282592773438f

/* Input bit patterns for atanpif endpoint repairs.  */
# define RVV_ATANPIF_FIX_NEG_0_U 0xbf7fb455u
# define RVV_ATANPIF_FIX_NEG_1_U 0xbf7e6ad9u
# define RVV_ATANPIF_FIX_NEG_2_U 0xbf78335fu
# define RVV_ATANPIF_FIX_NEG_3_U 0xbf7670b7u
# define RVV_ATANPIF_FIX_NEG_4_U 0xbf762a13u
# define RVV_ATANPIF_FIX_NEG_5_U 0xbf71578au
# define RVV_ATANPIF_FIX_NEG_6_U 0xbf70b248u
# define RVV_ATANPIF_FIX_NEG_7_U 0xbf6eda9du
# define RVV_ATANPIF_FIX_NEG_8_U 0xbf6e9333u
# define RVV_ATANPIF_FIX_NEG_9_U 0xbf6941b1u
# define RVV_ATANPIF_FIX_POS_0_U 0x3f77419du
# define RVV_ATANPIF_FIX_POS_1_U 0x3f70f56fu
# define RVV_ATANPIF_FIX_POS_2_U 0x3f72d2e7u
# define RVV_ATANPIF_FIX_POS_3_U 0x3f7349e2u
# define RVV_ATANPIF_FIX_POS_4_U 0x3f77c250u
# define RVV_ATANPIF_FIX_POS_5_U 0x3f7bf501u
# define RVV_ATANPIF_FIX_POS_6_U 0x3f7eead8u

# define RVV_ATANF_PI_F RVV_SLEEF_TRIG_PI_F
# define RVV_ATANF_PI2_F 0x1.921fb6p+0f
# define RVV_ATANF_PI4_F 0x1.921fb6p-1f
# define RVV_ATANF_PI2_HI_F 1.5707963705062866211f
# define RVV_ATANF_PI2_LO_F -4.3711388286737928865e-08f
# define RVV_ATANF_INV_PI_HI_F RVV_SLEEF_TRIG_ONE_OVER_PI
# define RVV_ATANF_INV_PI_LO_F 0x1.b93910p-27f
# define RVV_ATANF_TINY_F 2.9387372783541830947e-39f
# define RVV_ATANF_SCALE_F 16777216.0f
# define RVV_ATANF_NORM_PI_F 1.0f
# define RVV_ATANF_NORM_PI2_F 0.5f
# define RVV_ATANF_NORM_PI4_F 0.25f
# define RVV_ATANPIF_BIAS_F 0x1.0p-25f
#endif

/* LMUL-qualified helpers used by the clean single-precision atan family.
   Importer wrappers may explicitly activate the block, while direct
   skeleton users select it from V_MATH_FUNC.  Per-LMUL guards are required
   because v_math.h is re-included for every importer instantiation.  */
#ifndef RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atanf
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atanf 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atan2f 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atan2pif 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_INNER(name)                  \
  RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_##name
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG(name)                        \
  RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_INNER (name)
#endif

/* Shared coshf importer names.  The clean source-shaped helper and the
   source-original comparator use distinct function prefixes, but both need
   the same per-LMUL naming discipline.  Keep these names here so neither
   coshf skeleton carries a routine-local macro block.  */
#undef RVV_COSHF_SOURCE_MANGLE
#define RVV_COSHF_SOURCE_MANGLE(base) __PASTE3 (base, _m, LMUL)
#undef RVV_COSHF_EXPK2F_VF2
#define RVV_COSHF_EXPK2F_VF2 \
  __PASTE4 (coshf_expk2f_vf2_, LMUL, _, SIMDLEN)
#undef RVV_SLEEF_ORIG_COSHF_MANGLE
#define RVV_SLEEF_ORIG_COSHF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
#undef RVV_SLEEF_ORIG_COSHF_VPOW2I
#define RVV_SLEEF_ORIG_COSHF_VPOW2I \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_vpow2i_m)
#undef RVV_SLEEF_ORIG_COSHF_VLDEXP2
#define RVV_SLEEF_ORIG_COSHF_VLDEXP2 \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_vldexp2_m)
#undef RVV_SLEEF_ORIG_COSHF_EXPK2F
#define RVV_SLEEF_ORIG_COSHF_EXPK2F \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_expk2f_m)

/* Shared source-association helper names for the single-precision cospif
   finite graph.  The helper body remains in the clean skeleton, while its
   importer-qualified name is owned here.  */
#undef RVV_COSPIF_SOURCE_MANGLE
#define RVV_COSPIF_SOURCE_MANGLE(base) __PASTE3 (base, _m, LMUL)
#undef RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR
#define RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR \
  RVV_COSPIF_SOURCE_MANGLE (cospif_source_dfmul_pair_pair)

#if defined(RVV_LIBMVEC_ATAN2F_HELPERS)                                        \
  || RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M1
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M1
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M2
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M2
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M4
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M4
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M8
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M8
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
static inline VINT
RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
  (VFLOAT f0, VFLOAT f1, VINT x, VINT y, size_t vlen)
{
  VBOOL m = __riscv_vmflt (f0, f1, vlen);
  return __riscv_vmerge (y, x, m, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL)
  (VFLOAT d, VFLOAT m, size_t vlen)
{
  VUINT abs_d = __riscv_vand (F_AS_U (d), 0x7fffffff, vlen);
  VBOOL is_inf = __riscv_vmseq (abs_d, 0x7f800000, vlen);
  VUINT sign = __riscv_vand (F_AS_U (d), 0x80000000, vlen);
  VUINT out_bits = __riscv_vor (F_AS_U (m), sign, vlen);
  VFLOAT out = U_AS_F (out_bits);
  return __riscv_vmerge (VFMV_VF (0.0f, vlen), out, is_inf, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (atan2kf_u1, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT x_hi, VFLOAT x_lo, VFLOAT *r_hi,
   VFLOAT *r_lo, size_t vlen)
{
  VINT q = RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
    (x_hi, VFMV_VF (0.0f, vlen), VMVI_VX (-2, vlen), VMVI_VX (0, vlen),
     vlen);
  VBOOL p = __riscv_vmflt (x_hi, 0.0f, vlen);
  VUINT sign = VMVU_VX (0x80000000, vlen);
  VUINT flip = __riscv_vmerge (VMVU_VX (0, vlen), sign, p, vlen);
  x_hi = U_AS_F (__riscv_vxor (F_AS_U (x_hi), flip, vlen));
  x_lo = U_AS_F (__riscv_vxor (F_AS_U (x_lo), flip, vlen));

  q = RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
    (x_hi, y_hi, __riscv_vadd (q, 1, vlen), q, vlen);
  p = __riscv_vmflt (x_hi, y_hi, vlen);

  VFLOAT neg_x_hi, neg_x_lo;
  DFNEG_VF2_VF2 (x_hi, x_lo, neg_x_hi, neg_x_lo, vlen);

  VFLOAT s_hi, s_lo, t_hi, t_lo;
  VSEL_VF2_VO (p, neg_x_hi, neg_x_lo, y_hi, y_lo, s_hi, s_lo, vlen);
  VSEL_VF2_VO (p, y_hi, y_lo, x_hi, x_lo, t_hi, t_lo, vlen);

  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (t_hi, t_lo, t_hi, t_lo, vlen);

  VFLOAT u = VFMV_VF (RVV_ATANF_POLY0_F, vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY1_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY2_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY3_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY4_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY5_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY6_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY7_F, vlen), vlen);

  VFLOAT u_t = __riscv_vfmul (u, t_hi, vlen);
  VFLOAT poly_hi, poly_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_ATANF_POLY8_F, vlen), u_t, poly_hi, poly_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (t_hi, t_lo, poly_hi, poly_lo, t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, sum_hi, sum_lo, t_hi, t_lo, vlen);

  VFLOAT qf = __riscv_vfcvt_f (q, vlen);
  VFLOAT add_hi, add_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ATANF_PI2_HI_F, vlen),
     VFMV_VF (RVV_ATANF_PI2_LO_F, vlen), qf,
     add_hi, add_lo, vlen);

  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (add_hi, add_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}
#  undef RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
# endif
#endif

/* Source-shaped SLEEF atan-family helpers.  They are kept separately
   mangled from the clean helper block so the comparison objects retain the
   source operation graph while sharing the public math header.  */
#ifndef RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanf_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atan2f_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atan2pif_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_INNER(name)                  \
  RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_##name
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG(name)                        \
  RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_INNER (name)
#endif

#if defined(RVV_SLEEF_ATAN_HELPERS)                                           \
  || RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M1
#   define RVV_SLEEF_ATAN_HELPERS_M1
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M2
#   define RVV_SLEEF_ATAN_HELPERS_M2
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M4
#   define RVV_SLEEF_ATAN_HELPERS_M4
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M8
#   define RVV_SLEEF_ATAN_HELPERS_M8
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
  (VFLOAT f0, VFLOAT f1, VINT x, VINT y, size_t vlen)
{
  VBOOL m = __riscv_vmflt (f0, f1, vlen);
  return __riscv_vmerge (y, x, m, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_visinf2_m, LMUL)
  (VFLOAT d, VFLOAT m, size_t vlen)
{
  VUINT abs_d = __riscv_vand (F_AS_U (d), 0x7fffffff, vlen);
  VBOOL is_inf = __riscv_vmseq (abs_d, 0x7f800000, vlen);
  VUINT sign = __riscv_vand (F_AS_U (d), 0x80000000, vlen);
  VUINT out_bits = __riscv_vor (F_AS_U (m), sign, vlen);
  VUINT selected = __riscv_vmerge (VMVU_VX (0, vlen), out_bits, is_inf,
                                   vlen);
  return U_AS_F (selected);
}

static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan2kf_u1_m, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT x_hi, VFLOAT x_lo,
   VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VINT q = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
    (x_hi, VFMV_VF (0.0f, vlen), VMVI_VX (-2, vlen), VMVI_VX (0, vlen),
     vlen);
  VBOOL p = __riscv_vmflt (x_hi, 0.0f, vlen);
  VUINT sign = VMVU_VX (0x80000000, vlen);
  VUINT flip = __riscv_vmerge (VMVU_VX (0, vlen), sign, p, vlen);
  x_hi = U_AS_F (__riscv_vxor (F_AS_U (x_hi), flip, vlen));
  x_lo = U_AS_F (__riscv_vxor (F_AS_U (x_lo), flip, vlen));

  q = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
    (x_hi, y_hi, __riscv_vadd (q, 1, vlen), q, vlen);
  p = __riscv_vmflt (x_hi, y_hi, vlen);

  VFLOAT neg_x_hi, neg_x_lo;
  RVV_SLEEF_ORIG_DFNEG_VF2 (x_hi, x_lo, neg_x_hi, neg_x_lo, vlen);

  VFLOAT s_hi = __riscv_vmerge (y_hi, neg_x_hi, p, vlen);
  VFLOAT s_lo = __riscv_vmerge (y_lo, neg_x_lo, p, vlen);
  VFLOAT t_hi = __riscv_vmerge (x_hi, y_hi, p, vlen);
  VFLOAT t_lo = __riscv_vmerge (x_lo, y_lo, p, vlen);

  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (s_hi, s_lo, t_hi, t_lo,
                                s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (t_hi, t_lo, t_hi, t_lo, vlen);

  VFLOAT u = VFMV_VF (RVV_ATANF_POLY0_F, vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY1_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY2_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY3_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY4_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY5_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY6_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY7_F, vlen), vlen);

  VFLOAT u_t = VFMUL_VV (u, t_hi, vlen);
  VFLOAT poly_hi, poly_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_ATANF_POLY8_F, vlen), u_t, poly_hi, poly_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (t_hi, t_lo, poly_hi, poly_lo,
                                t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo,
                               sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (s_hi, s_lo, sum_hi, sum_lo,
                                t_hi, t_lo, vlen);

  VFLOAT qf = __riscv_vfcvt_f (q, vlen);
  VFLOAT add_hi, add_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ATANF_PI2_HI_F, vlen),
     VFMV_VF (RVV_ATANF_PI2_LO_F, vlen), qf,
     add_hi, add_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (add_hi, add_lo, t_hi, t_lo,
                                *r_hi, *r_lo, vlen);
}
#  undef RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
# endif
#endif
/* The double-precision cbrt importer keeps its source EH envelope in this
   shared helper so the skeleton contains only the operation graph and its
   data reference.  */
#ifndef RVV_MATH_IS_CBRT_cbrt
# define RVV_MATH_IS_CBRT_cbrt 1
# define RVV_MATH_IS_CBRT_INNER(name) RVV_MATH_IS_CBRT_##name
# define RVV_MATH_IS_CBRT(name) RVV_MATH_IS_CBRT_INNER (name)
#endif

#if RVV_MATH_IS_CBRT (V_MATH_FUNC)
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_exception_handling, LMUL)
  (VFLOAT *vx, VBOOL *special_args, VFLOAT *vy_special, VINT *n_adjust,
   size_t vlen)
{
  VUINT vclass = __riscv_vfclass (*vx, vlen);
  IDENTIFY (vclass, 0x399, *special_args, vlen);
  VBOOL denorm;
  IDENTIFY (vclass, 0x24, denorm, vlen);
  VBOOL special_n_denorm = __riscv_vmor (*special_args, denorm, vlen);
  *n_adjust = __riscv_vxor (*n_adjust, *n_adjust, vlen);
  if (__riscv_vcpop (special_n_denorm, vlen) > 0)
    {
      VFLOAT vx_normal = __riscv_vfmul
        (denorm, *vx, VFMV_VF (0x1.0p60, vlen), vlen);
      *vx = __riscv_vmerge (*vx, vx_normal, denorm, vlen);
      *n_adjust = __riscv_vmerge (*n_adjust, -20, denorm, vlen);
      *vy_special = __riscv_vfadd (*special_args, *vx, *vx, vlen);
      *vx = __riscv_vfmerge (*vx, fp_posOne, *special_args, vlen);
    }
}

/* Keep the cbrt polynomial operation graph unchanged for active finite lanes
   while preventing sanitized exceptional lanes from executing the FMA.  The
   initial scalar polynomial constants are materialized by the caller so this
   helper uses the same vector-vector fused operation as the common PSTEP
   spelling.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL)
  (VBOOL active, double coeff, VFLOAT x, VFLOAT poly, size_t vlen)
{
  return VFMADD_VV_M (active, poly, x, VFMV_VF (coeff, vlen), vlen);
}
#endif
