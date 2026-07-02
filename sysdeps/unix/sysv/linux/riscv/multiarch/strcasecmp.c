/* Multiple versions of strcasecmp.
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
/* Redefine strcasecmp so that the compiler won't complain about the type
   mismatch with the IFUNC selector in weak_alias, below.  */
# define strcasecmp __redirect_strcasecmp
# define __strcasecmp __redirect___strcasecmp
# include <stdint.h>
# include <string.h>
# undef strcasecmp
# undef __strcasecmp
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_strcasecmp) __strcasecmp;
extern __typeof (__redirect_strcasecmp) __strcasecmp_generic attribute_hidden;
extern __typeof (__redirect_strcasecmp) __strcasecmp_vector attribute_hidden;

static inline __typeof (__redirect_strcasecmp) *
select_strcasecmp_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    return __strcasecmp_vector;

  return __strcasecmp_generic;
}

riscv_libc_ifunc (__strcasecmp, select_strcasecmp_ifunc);
weak_alias (__strcasecmp, strcasecmp)
# ifdef SHARED
__hidden_ver1 (__strcasecmp, __GI___strcasecmp, __redirect___strcasecmp)
  __attribute__ ((visibility ("hidden")));
# endif
#else
# include <string/strcasecmp.c>
#endif
