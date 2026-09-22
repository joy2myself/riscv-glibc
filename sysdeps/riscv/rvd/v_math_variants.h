/* RVV libmvec LMUL/SIMDLEN variants.
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

#ifndef V_MATH_VARIANTS_H
#define V_MATH_VARIANTS_H

#define RVV_LIBMVEC_VARIANTS_F64(X)                                          \
  X (1, 2)                                                                   \
  X (1, 4)                                                                   \
  X (2, 2)                                                                   \
  X (2, 4)                                                                   \
  X (2, 8)                                                                   \
  X (4, 4)                                                                   \
  X (4, 8)                                                                   \
  X (4, 16)                                                                  \
  X (8, 8)                                                                   \
  X (8, 16)                                                                  \
  X (8, 32)

#define RVV_LIBMVEC_VARIANTS_F32(X)                                          \
  X (1, 4)                                                                   \
  X (1, 8)                                                                   \
  X (2, 4)                                                                   \
  X (2, 8)                                                                   \
  X (2, 16)                                                                  \
  X (4, 8)                                                                   \
  X (4, 16)                                                                  \
  X (4, 32)                                                                  \
  X (8, 16)                                                                  \
  X (8, 32)                                                                  \
  X (8, 64)

#endif
