/* Multiple versions of strcpy.
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
/* Redefine strcpy so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef strcpy
# define strcpy __redirect_strcpy
# include <string.h>
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_strcpy) __libc_strcpy;

extern __typeof (__redirect_strcpy) __strcpy_generic attribute_hidden;
extern __typeof (__redirect_strcpy) __strcpy_vector attribute_hidden;
extern __typeof (__redirect_strcpy) __strcpy_spacemit_x60 attribute_hidden;

static inline __typeof (__redirect_strcpy) *
select_strcpy_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    {
      if (tune == RISCV_CPU_TUNE_SPACEMIT_X60)
	return __strcpy_spacemit_x60;
      return __strcpy_vector;
    }

  return __strcpy_generic;
}

riscv_libc_ifunc (__libc_strcpy, select_strcpy_ifunc);

# undef strcpy
strong_alias (__libc_strcpy, strcpy);
# ifdef SHARED
__hidden_ver1 (strcpy, __GI_strcpy, __redirect_strcpy)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (strcpy);
# endif
#else
# include <string/strcpy.c>
#endif
