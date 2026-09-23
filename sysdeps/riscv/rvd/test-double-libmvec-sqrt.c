/* Vector ABI test for RISC-V libmvec sqrt.
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
#define LIBMVEC_FUNC sqrt
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
    0x1.6a09e667f3bcdp-3, 0x1.0000000000000p-2,
    0x1.3988e1409212ep-2, 0x1.6a09e667f3bcdp-2,
    0x1.94c583ada5b53p-2, 0x1.bb67ae8584caap-2,
    0x1.deeea11683f49p-2, 0x1.0000000000000p-1,
  };

#define LIBMVEC_ABI_TOLERANCE 0x1.6849b86a12b9bp-46
#define LIBMVEC_ORIGINAL sqrt_original
#include "test-vector-abi-arg1.h"
