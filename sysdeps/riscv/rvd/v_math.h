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
static const uint64_t class_posZero = 0x10;
static const uint64_t class_posDenorm = 0x20;
static const uint64_t class_negative = 0x7;

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
