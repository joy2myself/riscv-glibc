/* Vector ABI test for RISC-V libmvec tanf.
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
#define LIBMVEC_FUNC tanf
/* Finite ABI fixtures generated independently with 256-bit MPFR.
   This is an entry/lane smoke check, not the full accuracy suite.  */
static const LIBMVEC_TYPE libmvec_test_input[] =
  {
    0x1.0000000000000p-5f, 0x1.0000000000000p-4f,
    0x1.8000000000000p-4f, 0x1.0000000000000p-3f,
    0x1.4000000000000p-3f, 0x1.8000000000000p-3f,
    0x1.c000000000000p-3f, 0x1.0000000000000p-2f,
  };

static const LIBMVEC_TYPE libmvec_test_expected[] =
  {
    0x1.0015580000000p-5f, 0x1.0055780000000p-4f,
    0x1.8121040000000p-4f, 0x1.01577a0000000p-3f,
    0x1.42a13e0000000p-3f, 0x1.8490700000000p-3f,
    0x1.c7490a0000000p-3f, 0x1.05785a0000000p-2f,
  };

#define LIBMVEC_ABI_TOLERANCE 0x1.92a737110e454p-19
#define LIBMVEC_ORIGINAL tanf_original
#include "test-vector-abi-arg1.h"
