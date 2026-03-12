/* Enumerate available IFUNC implementations of a function.  RISCV version.
   Copyright (C) 2024-2026 Free Software Foundation, Inc.
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

#include <ifunc-impl-list.h>
#include <riscv-ifunc.h>
#include <string.h>
#include <cpu-features.h>
#include <ldsodefs.h>
#include <profile-detection.h>

size_t
__libc_ifunc_impl_list (const char *name, struct libc_ifunc_impl *array,
			size_t max)
{
  size_t i = max;

  const struct cpu_features *cpu_features = &GLRO (dl_riscv_cpu_features);
  bool __attribute__ ((unused)) fast_unaligned =
    riscv_hwprobe_fast_unaligned (cpu_features->cpuperf);
  bool __attribute__ ((unused)) rvv_enabled =
    (cpu_features->extensions & RISCV_HWPROBE_IMA_V) != 0;
  bool __attribute__ ((unused)) spacemit_x60_tuned =
    (cpu_features->tune_flags & RISCV_CPU_TUNE_SPACEMIT_X60) != 0;

  IFUNC_IMPL (i, name, memcpy,
	      IFUNC_IMPL_ADD (array, i, memcpy, rvv_enabled,
			      __memcpy_vector)
	      IFUNC_IMPL_ADD (array, i, memcpy, fast_unaligned,
			      __memcpy_noalignment)
	      IFUNC_IMPL_ADD (array, i, memcpy, 1, __memcpy_generic))

  IFUNC_IMPL (i, name, mempcpy,
	      IFUNC_IMPL_ADD (array, i, mempcpy, rvv_enabled,
			      __mempcpy_vector)
	      IFUNC_IMPL_ADD (array, i, mempcpy, 1, __mempcpy_generic))

  IFUNC_IMPL (i, name, memset,
	      IFUNC_IMPL_ADD (array, i, memset, rvv_enabled,
			      __memset_vector)
	      IFUNC_IMPL_ADD (array, i, memset, 1, __memset_generic))

  IFUNC_IMPL (i, name, strcat,
	      IFUNC_IMPL_ADD (array, i, strcat, rvv_enabled,
			      __strcat_vector)
	      IFUNC_IMPL_ADD (array, i, strcat, 1, __strcat_generic))

  IFUNC_IMPL (i, name, strcpy,
	      IFUNC_IMPL_ADD (array, i, strcpy, rvv_enabled && spacemit_x60_tuned,
			      __strcpy_spacemit_x60)
	      IFUNC_IMPL_ADD (array, i, strcpy, rvv_enabled,
			      __strcpy_vector)
	      IFUNC_IMPL_ADD (array, i, strcpy, 1, __strcpy_generic))

  IFUNC_IMPL (i, name, strlen,
	      IFUNC_IMPL_ADD (array, i, strlen, rvv_enabled,
			      __strlen_vector)
	      IFUNC_IMPL_ADD (array, i, strlen, 1, __strlen_generic))

  IFUNC_IMPL (i, name, strcmp,
	      IFUNC_IMPL_ADD (array, i, strcmp, rvv_enabled,
			      __strcmp_vector)
	      IFUNC_IMPL_ADD (array, i, strcmp, 1, __strcmp_generic))

  IFUNC_IMPL (i, name, strncmp,
	      IFUNC_IMPL_ADD (array, i, strncmp, rvv_enabled,
			      __strncmp_vector)
	      IFUNC_IMPL_ADD (array, i, strncmp, 1, __strncmp_generic))

  IFUNC_IMPL (i, name, strncasecmp,
	      IFUNC_IMPL_ADD (array, i, strncasecmp, rvv_enabled,
			      __strncasecmp_vector)
	      IFUNC_IMPL_ADD (array, i, strncasecmp, 1,
			      __strncasecmp_generic))

  IFUNC_IMPL (i, name, strncasecmp_l,
	      IFUNC_IMPL_ADD (array, i, strncasecmp_l, rvv_enabled,
			      __strncasecmp_l_vector)
	      IFUNC_IMPL_ADD (array, i, strncasecmp_l, 1,
			      __strncasecmp_l_generic))

  IFUNC_IMPL (i, name, memccpy,
	      IFUNC_IMPL_ADD (array, i, memccpy, rvv_enabled,
			      __memccpy_vector)
	      IFUNC_IMPL_ADD (array, i, memccpy, 1, __memccpy_generic))

  IFUNC_IMPL (i, name, memcmp,
	      IFUNC_IMPL_ADD (array, i, memcmp, rvv_enabled,
			      __memcmp_vector)
	      IFUNC_IMPL_ADD (array, i, memcmp, 1, __memcmp_generic))

  IFUNC_IMPL (i, name, __memcmpeq,
	      IFUNC_IMPL_ADD (array, i, __memcmpeq, rvv_enabled,
			      ____memcmpeq_vector)
	      IFUNC_IMPL_ADD (array, i, __memcmpeq, 1, ____memcmpeq_generic))

  IFUNC_IMPL (i, name, rawmemchr,
	      IFUNC_IMPL_ADD (array, i, rawmemchr, rvv_enabled,
			      __rawmemchr_vector)
	      IFUNC_IMPL_ADD (array, i, rawmemchr, 1,
			      __rawmemchr_generic))
  IFUNC_IMPL (i, name, memrchr,
	      IFUNC_IMPL_ADD (array, i, memrchr, rvv_enabled,
			      __memrchr_vector)
	      IFUNC_IMPL_ADD (array, i, memrchr, 1, __memrchr_generic))
  IFUNC_IMPL (i, name, stpncpy,
	      IFUNC_IMPL_ADD (array, i, stpncpy, rvv_enabled,
			      __stpncpy_vector)
	      IFUNC_IMPL_ADD (array, i, stpncpy, 1, __stpncpy_generic))

  IFUNC_IMPL (i, name, strcasecmp,
	      IFUNC_IMPL_ADD (array, i, strcasecmp, rvv_enabled,
			      __strcasecmp_vector)
	      IFUNC_IMPL_ADD (array, i, strcasecmp, 1,
			      __strcasecmp_generic))
  IFUNC_IMPL (i, name, strcasecmp_l,
	      IFUNC_IMPL_ADD (array, i, strcasecmp_l, rvv_enabled,
			      __strcasecmp_l_vector)
	      IFUNC_IMPL_ADD (array, i, strcasecmp_l, 1,
			      __strcasecmp_l_generic))
  IFUNC_IMPL (i, name, memchr,
	      IFUNC_IMPL_ADD (array, i, memchr, rvv_enabled,
			      __memchr_vector)
	      IFUNC_IMPL_ADD (array, i, memchr, 1, __memchr_generic))

  IFUNC_IMPL (i, name, strchr,
	      IFUNC_IMPL_ADD (array, i, strchr, rvv_enabled,
			      __strchr_vector)
	      IFUNC_IMPL_ADD (array, i, strchr, 1, __strchr_generic))

  IFUNC_IMPL (i, name, strchrnul,
	      IFUNC_IMPL_ADD (array, i, strchrnul, rvv_enabled,
			      __strchrnul_vector)
	      IFUNC_IMPL_ADD (array, i, strchrnul, 1, __strchrnul_generic))

  IFUNC_IMPL (i, name, strcspn,
	      IFUNC_IMPL_ADD (array, i, strcspn, rvv_enabled,
			      __strcspn_vector)
	      IFUNC_IMPL_ADD (array, i, strcspn, 1, __strcspn_generic))

  IFUNC_IMPL (i, name, strstr,
	      IFUNC_IMPL_ADD (array, i, strstr, rvv_enabled,
			      __strstr_vector)
	      IFUNC_IMPL_ADD (array, i, strstr, 1, __strstr_generic))

  IFUNC_IMPL (i, name, strrchr,
	      IFUNC_IMPL_ADD (array, i, strrchr, rvv_enabled,
			      __strrchr_vector)
	      IFUNC_IMPL_ADD (array, i, strrchr, 1, __strrchr_generic))

  IFUNC_IMPL (i, name, memmove,
	      IFUNC_IMPL_ADD (array, i, memmove, rvv_enabled,
			      __memmove_vector)
	      IFUNC_IMPL_ADD (array, i, memmove, 1, __memmove_generic))

  IFUNC_IMPL (i, name, strncat,
	      IFUNC_IMPL_ADD (array, i, strncat, rvv_enabled,
			      __strncat_vector)
	      IFUNC_IMPL_ADD (array, i, strncat, 1, __strncat_generic))

  IFUNC_IMPL (i, name, strncpy,
	      IFUNC_IMPL_ADD (array, i, strncpy, rvv_enabled,
			      __strncpy_vector)
	      IFUNC_IMPL_ADD (array, i, strncpy, 1, __strncpy_generic))

  IFUNC_IMPL (i, name, strnlen,
	      IFUNC_IMPL_ADD (array, i, strnlen, rvv_enabled,
			      __strnlen_vector)
	      IFUNC_IMPL_ADD (array, i, strnlen, 1, __strnlen_generic))

  IFUNC_IMPL (i, name, strpbrk,
	      IFUNC_IMPL_ADD (array, i, strpbrk, rvv_enabled,
			      __strpbrk_vector)
	      IFUNC_IMPL_ADD (array, i, strpbrk, 1, __strpbrk_generic))
  IFUNC_IMPL (i, name, strspn,
	      IFUNC_IMPL_ADD (array, i, strspn, rvv_enabled,
			      __strspn_vector)
	      IFUNC_IMPL_ADD (array, i, strspn, 1, __strspn_generic))
  return 0;
}
