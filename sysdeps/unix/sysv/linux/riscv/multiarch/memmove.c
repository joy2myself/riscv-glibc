/* Multiple versions of memmove.
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
/* Redefine memmove so that the compiler won't complain about the type
   mismatch with the IFUNC selector in strong_alias, below.  */
# undef memmove
# define memmove __redirect_memmove
# include <stdint.h>
# include <string.h>
# include <profile-ifunc-macros.h>
# include <riscv-ifunc.h>

extern __typeof (__redirect_memmove) __libc_memmove;

extern __typeof (__redirect_memmove) __memmove_generic attribute_hidden;
extern __typeof (__redirect_memmove) __memmove_c908 attribute_hidden;
extern __typeof (__redirect_memmove) __memmove_spacemit_x60 attribute_hidden;
extern __typeof (__redirect_memmove) __memmove_vector attribute_hidden;

static inline __typeof (__redirect_memmove) *
select_memmove_ifunc (uint64_t dl_hwcap, __riscv_hwprobe_t hwprobe_func)
{
  (void) dl_hwcap;
  (void) hwprobe_func;
  INIT_ARCH ();

  if (RISCV_PROFILE_COND (HAS_VECTOR (), V))
    {
      if (tune == RISCV_CPU_TUNE_SPACEMIT_X60)
	return __memmove_spacemit_x60;
      if (tune == RISCV_CPU_TUNE_C908)
	return __memmove_c908;
      return __memmove_vector;
    }

  return __memmove_generic;
}

riscv_libc_ifunc (__libc_memmove, select_memmove_ifunc);

# undef memmove
strong_alias (__libc_memmove, memmove);
# ifdef SHARED
__hidden_ver1 (memmove, __GI_memmove, __redirect_memmove)
  __attribute__ ((visibility ("hidden"))) __attribute_copy__ (memmove);
# endif
#else
# include <string/memmove.c>
#endif
