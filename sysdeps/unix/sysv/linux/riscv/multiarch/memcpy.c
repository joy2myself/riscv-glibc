/* Multiple versions of memcpy.
   All versions must be listed in ifunc-impl-list.c.
   Copyright (C) 2017-2026 Free Software Foundation, Inc.
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
/* Redefine memcpy so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef memcpy
# define memcpy __redirect_memcpy
# include <string.h>
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_memcpy) __libc_memcpy;

extern __typeof (__redirect_memcpy) __memcpy_generic attribute_hidden;
extern __typeof (__redirect_memcpy) __memcpy_noalignment attribute_hidden;
extern __typeof (__redirect_memcpy) __memcpy_vector attribute_hidden;

static inline __typeof (__redirect_memcpy) *
select_memcpy_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __memcpy_vector;

  if (RISCV_PROFILE_COND (fast_unaligned, RVA20))
    return __memcpy_noalignment;

  return __memcpy_generic;
}

riscv_libc_ifunc (__libc_memcpy, select_memcpy_ifunc);

# undef memcpy
strong_alias (__libc_memcpy, memcpy);
# ifdef SHARED
__hidden_ver1 (memcpy, __GI_memcpy, __redirect_memcpy)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (memcpy);
# endif
#else
# include <string/memcpy.c>
#endif
