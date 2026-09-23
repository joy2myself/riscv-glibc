/* Single-precision source-original SLEEF sinpif skeleton.

   This is a target-owned translation of SLEEF sinpifk/xsinpif_u05.  The
   signed quotient correction, polynomial constants, double-float
   reconstruction, and source late selections are kept in source order.
   Only the RVV representation, LMUL/VL plumbing, and glibc importer ABI are
   adapted; no SLEEF source file is included or linked.

   SLEEF's source files are distributed under their own license.  The
   corresponding attribution is retained in the source-reference tree.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT d, size_t vlen)
{
  return RVV_LIBMVEC_IMPL_NAME (sinpif_original_kernel, LMUL) (d, vlen);
}

#undef LMUL
#undef SIMDLEN
