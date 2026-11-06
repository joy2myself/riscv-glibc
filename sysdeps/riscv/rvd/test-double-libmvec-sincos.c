/* Vector ABI test for RISC-V libmvec sincos.
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
#define LIBMVEC_FUNC sincos
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
    0x1.ffeaaaeeee86fp-6, 0x1.ffaaaeeed4edbp-5,
    0x1.7f701032550e4p-4, 0x1.feaaeee86ee36p-4,
    0x1.3eb312c5d66cbp-3, 0x1.7dc102fbaf2b5p-3,
    0x1.bc6f84edc6199p-3, 0x1.faaeed4f31577p-3,
  };

static const LIBMVEC_TYPE libmvec_test_expected2[] =
  {
    0x1.ffc00155527d3p-1, 0x1.ff0015549f4d3p-1,
    0x1.fdc06bf7e6b9bp-1, 0x1.fc015527d5bd3p-1,
    0x1.f9c340a7cc428p-1, 0x1.f706bdf9ece1cp-1,
    0x1.f3cc7c3b3d16ep-1, 0x1.f01549f7deea1p-1,
  };

#define LIBMVEC_ABI_TOLERANCE 0x1.6849b86a12b9bp-46
#include "test-vector-abi-ret2.h"
