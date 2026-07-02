/* Multiple versions of stpncpy.
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
/* Redefine stpncpy so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef stpncpy
# define stpncpy __redirect_stpncpy
# define __stpncpy __redirect___stpncpy
# include <stdint.h>
# include <string.h>
# undef stpncpy
# undef __stpncpy
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_stpncpy) __stpncpy;
extern __typeof (__redirect_stpncpy) __stpncpy_generic attribute_hidden;
extern __typeof (__redirect_stpncpy) __stpncpy_vector attribute_hidden;

static inline __typeof (__redirect_stpncpy) *
select_stpncpy_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __stpncpy_vector;

  return __stpncpy_generic;
}

riscv_libc_ifunc (__stpncpy, select_stpncpy_ifunc);

# undef stpncpy
weak_alias (__stpncpy, stpncpy);
# ifdef SHARED
__hidden_ver1 (__stpncpy, __GI___stpncpy, __redirect___stpncpy)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (__stpncpy);
# endif
#else
# include <string/stpncpy.c>
#endif
