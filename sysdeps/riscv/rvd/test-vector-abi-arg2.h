/* Test for vector ABI with 2 arguments.
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

#include <riscv_vector.h>
#include <stddef.h>
#include <support/test-driver.h>

#include <v_math_variants.h>
#include <v_math_names.h>

/* Since libsupport_nonshared.a is placed before test-libmvec*.o, which
   defines do_test, reference support_test_main here to include it to
   avoid undefined reference to support_test_main.  The libmvec ABI test
   doesn't need other symbols in libsupport_nonshared.a.  */
__typeof (support_test_main) *support_test_main_p = support_test_main;

#define N 1000
LIBMVEC_TYPE x[N], a[N], b[N];

#define LIBMVEC_PASTE4(a, b, c, d) LIBMVEC_PASTE4_INNER (a, b, c, d)
#define LIBMVEC_PASTE4_INNER(a, b, c, d) a##b##c##d
#define LIBMVEC_PASTE5(a, b, c, d, e) LIBMVEC_PASTE5_INNER (a, b, c, d, e)
#define LIBMVEC_PASTE5_INNER(a, b, c, d, e) a##b##c##d##e
#define LIBMVEC_PASTE6(a, b, c, d, e, f)                                     \
  LIBMVEC_PASTE6_INNER (a, b, c, d, e, f)
#define LIBMVEC_PASTE6_INNER(a, b, c, d, e, f) a##b##c##d##e##f

#define RVV_LIBMVEC_IS_VALID_TYPE_double 1
#define RVV_LIBMVEC_IS_VALID_TYPE_float 1
#define RVV_LIBMVEC_IS_VALID_TYPE_INNER(t) RVV_LIBMVEC_IS_VALID_TYPE_##t
#define RVV_LIBMVEC_IS_VALID_TYPE(t) RVV_LIBMVEC_IS_VALID_TYPE_INNER (t)

#if !RVV_LIBMVEC_IS_VALID_TYPE (LIBMVEC_TYPE)
# error "Unsupported LIBMVEC_TYPE for RVV libmvec ABI test (expected float or double)"
#endif

#define RVV_LIBMVEC_ELEM_WIDTH_double 64
#define RVV_LIBMVEC_ELEM_WIDTH_float 32
#define RVV_LIBMVEC_ELEM_WIDTH_INNER(t) RVV_LIBMVEC_ELEM_WIDTH_##t
#define RVV_LIBMVEC_ELEM_WIDTH(t) RVV_LIBMVEC_ELEM_WIDTH_INNER (t)

#define LIBMVEC_ELEM_WIDTH RVV_LIBMVEC_ELEM_WIDTH (LIBMVEC_TYPE)

#define LIBMVEC_VTYPE_WIDTH(width, lmul)                                     \
  LIBMVEC_PASTE5 (vfloat, width, m, lmul, _t)
#define LIBMVEC_VSETVL_WIDTH(width, lmul)                                    \
  LIBMVEC_PASTE4 (__riscv_vsetvl_e, width, m, lmul)
#define LIBMVEC_VLE_WIDTH(width, lmul)                                       \
  LIBMVEC_PASTE6 (__riscv_vle, width, _v_f, width, m, lmul)
#define LIBMVEC_VSE_WIDTH(width, lmul)                                       \
  LIBMVEC_PASTE6 (__riscv_vse, width, _v_f, width, m, lmul)

#define LIBMVEC_VTYPE(lmul) LIBMVEC_VTYPE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VSETVL(lmul) LIBMVEC_VSETVL_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VLE(lmul) LIBMVEC_VLE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)
#define LIBMVEC_VSE(lmul) LIBMVEC_VSE_WIDTH (LIBMVEC_ELEM_WIDTH, lmul)

#define LIBMVEC_DECLARE_VARIANT(lmul, simdlen)                               \
  extern LIBMVEC_VTYPE (lmul) V_NAME_D2 (lmul, simdlen, LIBMVEC_FUNC)        \
    (LIBMVEC_VTYPE (lmul), LIBMVEC_VTYPE (lmul))

#define LIBMVEC_CALL_VARIANT(lmul, simdlen)                                  \
  do                                                                         \
  {                                                                          \
      size_t vl = LIBMVEC_VSETVL (lmul) (simdlen);                           \
      LIBMVEC_VTYPE (lmul) va = LIBMVEC_VLE (lmul) (a, vl);                  \
      LIBMVEC_VTYPE (lmul) vb = LIBMVEC_VLE (lmul) (b, vl);                  \
      LIBMVEC_VTYPE (lmul) vr = V_NAME_D2 (lmul, simdlen, LIBMVEC_FUNC)      \
        (va, vb);                                                            \
      LIBMVEC_VSE (lmul) (x, vr, vl);                                        \
  }                                                                          \
  while (0)

#if LIBMVEC_ELEM_WIDTH == 64
# define LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F64 (X)
#elif LIBMVEC_ELEM_WIDTH == 32
# define LIBMVEC_VARIANTS(X) RVV_LIBMVEC_VARIANTS_F32 (X)
#else
# error "Unsupported LIBMVEC_ELEM_WIDTH for RVV libmvec ABI test"
#endif

#define X(lmul, simdlen) LIBMVEC_DECLARE_VARIANT (lmul, simdlen);
LIBMVEC_VARIANTS (X)
#undef X

int
test_vector_abi (void)
{
  int i;
  for (i = 0; i < N; i++)
    {
      a[i] = i / 3;
      b[i] = i / 5 + 1;
    }

#define X(lmul, simdlen) LIBMVEC_CALL_VARIANT (lmul, simdlen);
  LIBMVEC_VARIANTS (X)
#undef X

  return 0;
}

#include <test-vector-abi.h>
