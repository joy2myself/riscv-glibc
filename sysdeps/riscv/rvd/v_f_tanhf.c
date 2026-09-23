/* Single-precision vector tanhf function.
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

#define LIBMVEC_TYPE float
#define LIBMVEC_FUNC tanhf
#define LIBMVEC_IMPL "v_f_tanhf_skeleton.c"
/* Shared components for this translation unit and all LMUL instances.  */
#define RVV_NEED_CLEAN_F32_HYPERBOLIC_EXP 1
#include "v_math_importer.h"
