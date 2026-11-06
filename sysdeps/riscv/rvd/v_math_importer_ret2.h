/* RVV libmvec 1-input/2-output implementation importer.
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

#include <v_math_variants.h>

#ifndef LIBMVEC_TYPE
# error "LIBMVEC_TYPE must be defined for RVV libmvec import"
#endif
#ifndef LIBMVEC_FUNC
# error "LIBMVEC_FUNC must be defined for RVV libmvec import"
#endif
#ifndef LIBMVEC_IMPL
# error "LIBMVEC_IMPL must be defined for RVV libmvec import"
#endif

#define RVV_LIBMVEC_IS_VALID_TYPE_double 1
#define RVV_LIBMVEC_IS_VALID_TYPE_float 1
#define RVV_LIBMVEC_IS_VALID_TYPE_INNER(t) RVV_LIBMVEC_IS_VALID_TYPE_##t
#define RVV_LIBMVEC_IS_VALID_TYPE(t) RVV_LIBMVEC_IS_VALID_TYPE_INNER (t)

#if !RVV_LIBMVEC_IS_VALID_TYPE (LIBMVEC_TYPE)
# error "Unsupported LIBMVEC_TYPE for RVV libmvec import (expected float or double)"
#endif

#define RVV_LIBMVEC_IS_DOUBLE_double 1
#define RVV_LIBMVEC_IS_DOUBLE_float 0
#define RVV_LIBMVEC_IS_DOUBLE_INNER(t) RVV_LIBMVEC_IS_DOUBLE_##t
#define RVV_LIBMVEC_IS_DOUBLE(t) RVV_LIBMVEC_IS_DOUBLE_INNER (t)

#define RVV_LIBMVEC_BIT_WIDTH_double 64
#define RVV_LIBMVEC_BIT_WIDTH_float 32
#define RVV_LIBMVEC_BIT_WIDTH_INNER(t) RVV_LIBMVEC_BIT_WIDTH_##t
#define RVV_LIBMVEC_BIT_WIDTH(t) RVV_LIBMVEC_BIT_WIDTH_INNER (t)

#define BIT_WIDTH RVV_LIBMVEC_BIT_WIDTH (LIBMVEC_TYPE)

#if RVV_LIBMVEC_IS_DOUBLE (LIBMVEC_TYPE)

# define RVV_LIBMVEC_IMPL_NAME_INNER(fun, lmul)                             \
  __rvv_##fun##_w64_m##lmul##_impl
# define RVV_LIBMVEC_IMPL_NAME(fun, lmul)                                   \
  RVV_LIBMVEC_IMPL_NAME_INNER (fun, lmul)

# define RVV_LIBMVEC_EXPORT_NAME(lmul, simdlen, fun)                        \
  V_NAME_D1_D1 (lmul, simdlen, fun)
# define RVV_LIBMVEC_VFLOAT(lmul) vfloat64m##lmul##_t
# define RVV_LIBMVEC_VSETVL(lmul) __riscv_vsetvl_e64m##lmul

# define RVV_LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F64 (X)

#else

# define RVV_LIBMVEC_IMPL_NAME_INNER(fun, lmul)                             \
  __rvv_##fun##_w32_m##lmul##_impl
# define RVV_LIBMVEC_IMPL_NAME(fun, lmul)                                   \
  RVV_LIBMVEC_IMPL_NAME_INNER (fun, lmul)

# define RVV_LIBMVEC_EXPORT_NAME(lmul, simdlen, fun)                        \
  V_NAME_F1_F1 (lmul, simdlen, fun)
# define RVV_LIBMVEC_VFLOAT(lmul) vfloat32m##lmul##_t
# define RVV_LIBMVEC_VSETVL(lmul) __riscv_vsetvl_e32m##lmul

# define RVV_LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F32 (X)

#endif

#define V_MATH_FUNC LIBMVEC_FUNC
#define V_MATH_IMPL LIBMVEC_IMPL

#undef LMUL
#define LMUL 1
#include V_MATH_IMPL

#undef LMUL
#define LMUL 2
#include V_MATH_IMPL

#undef LMUL
#define LMUL 4
#include V_MATH_IMPL

#undef LMUL
#define LMUL 8
#include V_MATH_IMPL

#define RVV_LIBMVEC_DEFINE_WRAPPER(lmul, simdlen)                             \
  void                                                                        \
  RVV_LIBMVEC_EXPORT_NAME (lmul, simdlen, V_MATH_FUNC)                        \
    (RVV_LIBMVEC_VFLOAT (lmul) x, RVV_LIBMVEC_VFLOAT (lmul) *out0,            \
     RVV_LIBMVEC_VFLOAT (lmul) *out1)                                         \
  {                                                                           \
    size_t vlen = RVV_LIBMVEC_VSETVL (lmul) (simdlen);                        \
    RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, lmul) (x, out0, out1, vlen);          \
  }

#define X(lmul, simdlen) RVV_LIBMVEC_DEFINE_WRAPPER (lmul, simdlen)
RVV_LIBMVEC_VARIANTS (X)
#undef X

#undef RVV_LIBMVEC_DEFINE_WRAPPER
#undef RVV_LIBMVEC_VARIANTS
#undef RVV_LIBMVEC_VSETVL
#undef RVV_LIBMVEC_VFLOAT
#undef RVV_LIBMVEC_EXPORT_NAME
#undef RVV_LIBMVEC_IMPL_NAME
#undef RVV_LIBMVEC_IMPL_NAME_INNER

#undef V_MATH_IMPL
#undef V_MATH_FUNC

#undef BIT_WIDTH
#undef RVV_LIBMVEC_BIT_WIDTH
#undef RVV_LIBMVEC_BIT_WIDTH_INNER
#undef RVV_LIBMVEC_BIT_WIDTH_double
#undef RVV_LIBMVEC_BIT_WIDTH_float

#undef RVV_LIBMVEC_IS_DOUBLE
#undef RVV_LIBMVEC_IS_DOUBLE_INNER
#undef RVV_LIBMVEC_IS_DOUBLE_double
#undef RVV_LIBMVEC_IS_DOUBLE_float

#undef RVV_LIBMVEC_IS_VALID_TYPE
#undef RVV_LIBMVEC_IS_VALID_TYPE_INNER
#undef RVV_LIBMVEC_IS_VALID_TYPE_double
#undef RVV_LIBMVEC_IS_VALID_TYPE_float
