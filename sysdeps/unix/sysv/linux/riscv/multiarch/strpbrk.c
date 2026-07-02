/* Multiple versions of strpbrk.
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
/* Redefine strpbrk so that the compiler won't complain about the type
   mismatch with the IFUNC selector below.  */
# define strpbrk __redirect_strpbrk
# include <stdint.h>
# include <string.h>
# undef strpbrk
# include <ifunc-init.h>
# include <riscv-ifunc.h>
# include <sys/hwprobe.h>

extern __typeof (__redirect_strpbrk) strpbrk;
extern __typeof (__redirect_strpbrk) __strpbrk_generic attribute_hidden;
extern __typeof (__redirect_strpbrk) __strpbrk_vector attribute_hidden;

static inline __typeof (__redirect_strpbrk) *
select_strpbrk_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  if (riscv_hwprobe_has_vector (hwprobe_func))
    return __strpbrk_vector;
  return __strpbrk_generic;
}

riscv_libc_ifunc (strpbrk, select_strpbrk_ifunc);
# ifdef SHARED
__hidden_ver1 (strpbrk, __GI_strpbrk, __redirect_strpbrk)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (strpbrk);
# endif
#else
# include <string/strpbrk.c>
#endif
