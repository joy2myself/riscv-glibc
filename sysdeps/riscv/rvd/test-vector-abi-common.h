/* Shared support for the RVV libmvec ABI smoke tests.
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

#ifndef TEST_VECTOR_ABI_COMMON_H
#define TEST_VECTOR_ABI_COMMON_H

#include <math.h>
#include <riscv_vector.h>
#include <stddef.h>
#include <stdio.h>
#include <support/check.h>
#include <support/test-driver.h>
#include <v_math_names.h>
#include <v_math_variants.h>

#define LIBMVEC_PASTE4(a, b, c, d) LIBMVEC_PASTE4_INNER (a, b, c, d)
#define LIBMVEC_PASTE4_INNER(a, b, c, d) a##b##c##d
#define LIBMVEC_PASTE5(a, b, c, d, e) LIBMVEC_PASTE5_INNER (a, b, c, d, e)
#define LIBMVEC_PASTE5_INNER(a, b, c, d, e) a##b##c##d##e
#define LIBMVEC_PASTE6(a, b, c, d, e, f) \
  LIBMVEC_PASTE6_INNER (a, b, c, d, e, f)
#define LIBMVEC_PASTE6_INNER(a, b, c, d, e, f) a##b##c##d##e##f
#define LIBMVEC_STRINGIFY(x) LIBMVEC_STRINGIFY_INNER (x)
#define LIBMVEC_STRINGIFY_INNER(x) #x

#define RVV_LIBMVEC_ELEM_WIDTH_double 64
#define RVV_LIBMVEC_ELEM_WIDTH_float 32
#define RVV_LIBMVEC_ELEM_WIDTH_INNER(t) RVV_LIBMVEC_ELEM_WIDTH_##t
#define RVV_LIBMVEC_ELEM_WIDTH(t) RVV_LIBMVEC_ELEM_WIDTH_INNER (t)
#define LIBMVEC_ELEM_WIDTH RVV_LIBMVEC_ELEM_WIDTH (LIBMVEC_TYPE)

#define LIBMVEC_VTYPE_WIDTH(width, lmul) \
  LIBMVEC_PASTE5 (vfloat, width, m, lmul, _t)
#define LIBMVEC_VSETVL_WIDTH(width, lmul) \
  LIBMVEC_PASTE4 (__riscv_vsetvl_e, width, m, lmul)
#define LIBMVEC_VLE_WIDTH(width, lmul) \
  LIBMVEC_PASTE6 (__riscv_vle, width, _v_f, width, m, lmul)
#define LIBMVEC_VSE_WIDTH(width, lmul) \
  LIBMVEC_PASTE6 (__riscv_vse, width, _v_f, width, m, lmul)
#define LIBMVEC_VTYPE(lmul) LIBMVEC_VTYPE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VSETVL(lmul) LIBMVEC_VSETVL_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VLE(lmul) LIBMVEC_VLE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VSE(lmul) LIBMVEC_VSE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)

#if LIBMVEC_ELEM_WIDTH == 64
# define LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F64 (X)
# define LIBMVEC_RET2_NAME(lmul, simdlen, fun) \
  V_NAME_D1_D1 (lmul, simdlen, fun)
#elif LIBMVEC_ELEM_WIDTH == 32
# define LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F32 (X)
# define LIBMVEC_RET2_NAME(lmul, simdlen, fun) \
  V_NAME_F1_F1 (lmul, simdlen, fun)
#else
# error "Unsupported LIBMVEC_TYPE for RVV libmvec ABI test"
#endif

#define LIBMVEC_MAX_LANES 64
#define LIBMVEC_NSAMPLES \
  (sizeof (libmvec_test_input) / sizeof (libmvec_test_input[0]))

_Static_assert (LIBMVEC_NSAMPLES >= 2, "distinct lane samples required");
_Static_assert (sizeof (libmvec_test_input)
                == sizeof (libmvec_test_expected), "fixture sizes differ");

static unsigned int libmvec_abi_calls;
extern int test_vector_abi (void);

/* OUTPUT has a guard before the active lanes and an untouched tail after
   them.  These check the test's scalar stores, not the pointer-output
   implementation's memory bounds.  Numerical bounds are deliberately loose
   ABI smoke bounds; the full accuracy suite remains a separate test.  */
static void
libmvec_check_output (const LIBMVEC_TYPE *output,
                      const LIBMVEC_TYPE *expected, size_t first,
                      size_t lanes, const char *symbol)
{
  TEST_VERIFY (output[0] == 12345);
  for (size_t i = lanes + 1; i < LIBMVEC_MAX_LANES + 2; i++)
    TEST_VERIFY (output[i] == 12345);

  for (size_t i = 0; i < lanes; i++)
    {
      size_t sample = (first + i) % LIBMVEC_NSAMPLES;
      double got = output[i + 1];
      double ref = expected[sample];
      double bound = LIBMVEC_ABI_TOLERANCE * fmax (1.0, fabs (ref));
      if (!isfinite (got) || fabs (got - ref) > bound)
        FAIL_EXIT1 ("%s lane %zu sample %zu: got %a, expected %a",
                    symbol, i, sample, got, ref);
    }
}

static int
libmvec_abi_result (void)
{
  return libmvec_abi_calls == 0 ? EXIT_UNSUPPORTED : 0;
}

#endif
