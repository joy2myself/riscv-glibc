/* Vector ABI test for RISC-V libmvec sinpi.
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

#define LIBMVEC_TYPE double
#define LIBMVEC_FUNC sinpi
/* Finite ABI fixtures generated independently with 256-bit MPFR.
   This is an entry/lane smoke check, not the full accuracy suite.  */
static const LIBMVEC_TYPE libmvec_test_input[] =
  {
    0x1.0000000000000p-5, 0x1.0000000000000p-4,
    0x1.8000000000000p-4, 0x1.0000000000000p-3,
    0x1.4000000000000p-3, 0x1.8000000000000p-3,
    0x1.c000000000000p-3, 0x1.0000000000000p-2,
  };

static const LIBMVEC_TYPE libmvec_test_expected[] =
  {
    0x1.917a6bc29b42cp-4, 0x1.8f8b83c69a60bp-3,
    0x1.294062ed59f06p-2, 0x1.87de2a6aea963p-2,
    0x1.e2b5d3806f63bp-2, 0x1.1c73b39ae68c8p-1,
    0x1.44cf325091dd6p-1, 0x1.6a09e667f3bcdp-1,
  };

#define LIBMVEC_ABI_TOLERANCE 0x1.6849b86a12b9bp-46
#include "test-vector-abi-arg1.h"
