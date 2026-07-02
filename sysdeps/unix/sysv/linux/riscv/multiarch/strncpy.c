/* Multiple versions of strncpy.
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
/* Redefine strncpy so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef strncpy
# define strncpy __redirect_strncpy
# include <stdint.h>
# include <string.h>
# include <ifunc-init.h>
# include <riscv-ifunc.h>
# include <sys/hwprobe.h>

extern __typeof (__redirect_strncpy) __libc_strncpy;

extern __typeof (__redirect_strncpy) __strncpy_generic attribute_hidden;
extern __typeof (__redirect_strncpy) __strncpy_vector attribute_hidden;

static inline __typeof (__redirect_strncpy) *
select_strncpy_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  unsigned long long int v;
  if (__riscv_hwprobe_one (hwprobe_func, RISCV_HWPROBE_KEY_IMA_EXT_0, &v) == 0
      && (v & RISCV_HWPROBE_IMA_V) == RISCV_HWPROBE_IMA_V)
    return __strncpy_vector;
  return __strncpy_generic;
}

riscv_libc_ifunc (__libc_strncpy, select_strncpy_ifunc);

# undef strncpy
strong_alias (__libc_strncpy, strncpy);
# ifdef SHARED
__hidden_ver1 (strncpy, __GI_strncpy, __redirect_strncpy)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (strncpy);
# endif
#else
# include <string/strncpy.c>
#endif
