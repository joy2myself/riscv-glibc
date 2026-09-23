/* Vector ABI test for RISC-V libmvec sinh.
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
#define LIBMVEC_FUNC sinh
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
    0x1.000aaacccd00dp-5, 0x1.002aacccd9cddp-4,
    0x1.8090103411660p-4, 0x1.00aaccd00d2f1p-3,
    0x1.414dbd8f81999p-3, 0x1.8241036ac51ddp-3,
    0x1.c394db89e8f7fp-3, 0x1.02accd9d08102p-2,
  };

#define LIBMVEC_ABI_TOLERANCE 0x1.6849b86a12b9bp-46
#include "test-vector-abi-arg1.h"
