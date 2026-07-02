/* Multiple versions of strspn.
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
/* Redefine strspn so that the compiler won't complain about the type
   mismatch with the IFUNC selector below.  */
# define strspn __redirect_strspn
# include <stdint.h>
# include <string.h>
# undef strspn
# include <ifunc-init.h>
# include <riscv-ifunc.h>
# include <sys/hwprobe.h>

extern __typeof (__redirect_strspn) strspn;
extern __typeof (__redirect_strspn) __strspn_generic attribute_hidden;
extern __typeof (__redirect_strspn) __strspn_vector attribute_hidden;

static inline __typeof (__redirect_strspn) *
select_strspn_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  if (riscv_hwprobe_has_vector (hwprobe_func))
    return __strspn_vector;
  return __strspn_generic;
}

riscv_libc_ifunc (strspn, select_strspn_ifunc);
# ifdef SHARED
__hidden_ver1 (strspn, __GI_strspn, __redirect_strspn)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (strspn);
# endif
#else
# include <string/strspn.c>
#endif
