/* RVV libmvec name helpers.
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

#ifndef V_MATH_NAMES_H
#define V_MATH_NAMES_H

#ifndef V_NAME_F1
# define V_NAME_F1_INNER(lmul, simdlen, fun) _ZGVr##lmul##N##simdlen##v_##fun
# define V_NAME_F1(lmul, simdlen, fun) V_NAME_F1_INNER (lmul, simdlen, fun)
#endif

#ifndef V_NAME_F2
# define V_NAME_F2_INNER(lmul, simdlen, fun) _ZGVr##lmul##N##simdlen##vv_##fun
# define V_NAME_F2(lmul, simdlen, fun) V_NAME_F2_INNER (lmul, simdlen, fun)
#endif

#ifndef V_NAME_D1
# define V_NAME_D1_INNER(lmul, simdlen, fun) _ZGVr##lmul##N##simdlen##v_##fun
# define V_NAME_D1(lmul, simdlen, fun) V_NAME_D1_INNER (lmul, simdlen, fun)
#endif

#ifndef V_NAME_D2
# define V_NAME_D2_INNER(lmul, simdlen, fun) _ZGVr##lmul##N##simdlen##vv_##fun
# define V_NAME_D2(lmul, simdlen, fun) V_NAME_D2_INNER (lmul, simdlen, fun)
#endif

#ifndef V_NAME_D1_D1
# define V_NAME_D1_D1_INNER(lmul, simdlen, fun)                              \
  _ZGVr##lmul##N##simdlen##vl8l8_##fun
# define V_NAME_D1_D1(lmul, simdlen, fun)                                    \
  V_NAME_D1_D1_INNER (lmul, simdlen, fun)
#endif

#endif
