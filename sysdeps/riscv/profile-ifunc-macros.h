/* Common ifunc selection utils
   All versions must be listed in ifunc-impl-list.c.
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

#ifndef _PROFILE_IFUNC_MACROS_H
#define _PROFILE_IFUNC_MACROS_H

#include <profile-level.h>

/* Only list implementations at the level of the minimum build profile or
   higher. I.e if built with profile lower than RVA20, then include all
   implementations. On the other hand if built with profile=RVA22/RVA23 only
   include RVA22/RVA23 implementations. If there is no implementation at or
   above the minimum build profile level, then include the highest profile
   level implementation.  */
#if MINIMUM_RISCV_PROFILE_LEVEL <= 23
# define RISCV_IFUNC_IMPL_ADD_RVA23(...) IFUNC_IMPL_ADD (__VA_ARGS__)
#else
# define RISCV_IFUNC_IMPL_ADD_RVA23(...)
#endif

#if MINIMUM_RISCV_PROFILE_LEVEL <= 22
# define RISCV_IFUNC_IMPL_ADD_RVA22(...) IFUNC_IMPL_ADD (__VA_ARGS__)
#else
# define RISCV_IFUNC_IMPL_ADD_RVA22(...)
#endif

#if MINIMUM_RISCV_PROFILE_LEVEL <= 20
# define RISCV_IFUNC_IMPL_ADD_RVA20(...) IFUNC_IMPL_ADD (__VA_ARGS__)
#else
# define RISCV_IFUNC_IMPL_ADD_RVA20(...)
#endif

#if MINIMUM_RISCV_PROFILE_LEVEL <= 0
# define RISCV_IFUNC_IMPL_ADD_INIT(...) IFUNC_IMPL_ADD (__VA_ARGS__)
#else
# define RISCV_IFUNC_IMPL_ADD_INIT(...)
#endif

#define RISCV_PROFILE_COND(cond, name)                                    \
  (((name##_RISCV_PROFILE_LEVEL) <= MINIMUM_RISCV_PROFILE_LEVEL)           \
   || (cond))

#endif
