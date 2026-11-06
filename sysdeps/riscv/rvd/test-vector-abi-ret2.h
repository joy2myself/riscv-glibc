/* Test for vector ABI with 1 input and 2 outputs.
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

#include <test-vector-abi-common.h>

_Static_assert (sizeof (libmvec_test_input)
                == sizeof (libmvec_test_expected2), "fixture sizes differ");

#define LIBMVEC_DECLARE_VARIANT(lmul, simdlen, fun) \
  _Static_assert (simdlen <= LIBMVEC_MAX_LANES, "lane buffer too small"); \
  extern void LIBMVEC_RET2_NAME (lmul, simdlen, fun) \
    (LIBMVEC_VTYPE (lmul), LIBMVEC_VTYPE (lmul) *, LIBMVEC_VTYPE (lmul) *)

#define X(lmul, simdlen) \
  LIBMVEC_DECLARE_VARIANT (lmul, simdlen, LIBMVEC_FUNC);
LIBMVEC_VARIANTS (X)
#undef X
#ifdef LIBMVEC_ORIGINAL
# define X(lmul, simdlen) \
  LIBMVEC_DECLARE_VARIANT (lmul, simdlen, LIBMVEC_ORIGINAL);
LIBMVEC_VARIANTS (X)
# undef X
#endif

#define LIBMVEC_CALL_VARIANT(lmul, simdlen, fun) \
  do \
    { \
      const char *symbol = \
        LIBMVEC_STRINGIFY (LIBMVEC_RET2_NAME (lmul, simdlen, fun)); \
      size_t vl = LIBMVEC_VSETVL (lmul) (simdlen); \
      if (vl != simdlen) \
        { \
          printf ("SKIP %s: requested VL %u, available %zu\n", \
                  symbol, (unsigned int) simdlen, vl); \
          break; \
        } \
      libmvec_abi_calls++; \
      for (size_t first = 0; first < LIBMVEC_NSAMPLES; first += simdlen) \
        { \
          LIBMVEC_TYPE input[LIBMVEC_MAX_LANES]; \
          LIBMVEC_TYPE output[LIBMVEC_MAX_LANES + 2]; \
          LIBMVEC_TYPE output2[LIBMVEC_MAX_LANES + 2]; \
          for (size_t i = 0; i < simdlen; i++) \
            input[i] = libmvec_test_input[(first + i) % LIBMVEC_NSAMPLES]; \
          for (size_t i = 0; i < LIBMVEC_MAX_LANES + 2; i++) \
            output[i] = output2[i] = 12345; \
          vl = LIBMVEC_VSETVL (lmul) (simdlen); \
          LIBMVEC_VTYPE (lmul) vx = LIBMVEC_VLE (lmul) (input, vl); \
          LIBMVEC_VTYPE (lmul) vr = LIBMVEC_VLE (lmul) (output + 1, vl); \
          LIBMVEC_VTYPE (lmul) vr2 = LIBMVEC_VLE (lmul) (output2 + 1, vl); \
          LIBMVEC_RET2_NAME (lmul, simdlen, fun) (vx, &vr, &vr2); \
          LIBMVEC_VSE (lmul) (output + 1, vr, vl); \
          LIBMVEC_VSE (lmul) (output2 + 1, vr2, vl); \
          libmvec_check_output (output, libmvec_test_expected, first, \
                                simdlen, symbol); \
          libmvec_check_output (output2, libmvec_test_expected2, first, \
                                simdlen, symbol); \
        } \
    } \
  while (0)

int
test_vector_abi (void)
{
#define X(lmul, simdlen) LIBMVEC_CALL_VARIANT (lmul, simdlen, LIBMVEC_FUNC);
  LIBMVEC_VARIANTS (X)
#undef X
#ifdef LIBMVEC_ORIGINAL
# define X(lmul, simdlen) \
  LIBMVEC_CALL_VARIANT (lmul, simdlen, LIBMVEC_ORIGINAL);
  LIBMVEC_VARIANTS (X)
# undef X
#endif
  return libmvec_abi_result ();
}

#include <test-vector-abi.h>
