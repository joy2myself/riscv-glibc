/* Multiple versions of mempcpy.
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
/* Redefine mempcpy so that the compiler won't complain about the type
   mismatch with the IFUNC selector in weak_alias, below.  */
# define mempcpy __redirect_mempcpy
# define __mempcpy __redirect___mempcpy
# define NO_MEMPCPY_STPCPY_REDIRECT
# define __NO_STRING_INLINES
# include <stdint.h>
# include <string.h>
# undef mempcpy
# undef __mempcpy
# include <ifunc-init.h>
# include <riscv-ifunc.h>
# include <sys/hwprobe.h>

extern __typeof (__redirect_mempcpy) __mempcpy;
extern __typeof (__redirect_mempcpy) __mempcpy_generic attribute_hidden;
extern __typeof (__redirect_mempcpy) __mempcpy_vector attribute_hidden;

static inline __typeof (__redirect_mempcpy) *
select_mempcpy_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  if (riscv_hwprobe_has_vector (hwprobe_func))
    return __mempcpy_vector;
  return __mempcpy_generic;
}

riscv_libc_ifunc (__mempcpy, select_mempcpy_ifunc);
weak_alias (__mempcpy, mempcpy)
# ifdef SHARED
__hidden_ver1 (__mempcpy, __GI___mempcpy, __redirect___mempcpy)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (__mempcpy);
__hidden_ver1 (__mempcpy, __GI_mempcpy, __redirect_mempcpy)
  __attribute__ ((weak, visibility ("hidden"))) __attribute_copy__ (mempcpy);
# endif
#else
# include <string/mempcpy.c>
#endif
