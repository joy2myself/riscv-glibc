/* Multiple versions of memrchr.
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

#if IS_IN (libc)
/* Redefine memrchr so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# define memrchr __redirect_memrchr
# include <stdint.h>
# include <string.h>
# undef memrchr
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_memrchr) __memrchr;
extern __typeof (__redirect_memrchr) __memrchr_generic attribute_hidden;
extern __typeof (__redirect_memrchr) __memrchr_vector attribute_hidden;

static inline __typeof (__redirect_memrchr) *
select_memrchr_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __memrchr_vector;

  return __memrchr_generic;
}

riscv_libc_ifunc (__memrchr, select_memrchr_ifunc);
libc_hidden_def (__memrchr)
weak_alias (__memrchr, memrchr)
#else
# include <string/memrchr.c>
#endif
