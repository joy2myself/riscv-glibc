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
#include <profile-ifunc-macros.h>

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
    cpu_features->tune == RISCV_CPU_TUNE_SPACEMIT_X60;

  IFUNC_IMPL (i, name, memcpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memcpy,
					  rvv_enabled && spacemit_x60_tuned,
					  __memcpy_spacemit_x60)
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memcpy, rvv_enabled,
					  __memcpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memcpy, fast_unaligned,
					  __memcpy_noalignment)
	      RISCV_IFUNC_IMPL_ADD_INIT (array, i, memcpy, 1,
					 __memcpy_generic))

  IFUNC_IMPL (i, name, mempcpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, mempcpy, rvv_enabled,
					  __mempcpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, mempcpy, 1,
					  __mempcpy_generic))

  IFUNC_IMPL (i, name, memset,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memset,
					  rvv_enabled && spacemit_x60_tuned,
					  __memset_spacemit_x60)
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memset, rvv_enabled,
					  __memset_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memset, 1,
					  __memset_generic))

  IFUNC_IMPL (i, name, strcat,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strcat, rvv_enabled,
					  __strcat_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strcat, 1,
					  __strcat_generic))

  IFUNC_IMPL (i, name, strcpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strcpy,
					  rvv_enabled && spacemit_x60_tuned,
					  __strcpy_spacemit_x60)
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strcpy, rvv_enabled,
					  __strcpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strcpy, 1,
					  __strcpy_generic))

  IFUNC_IMPL (i, name, strlen,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strlen, rvv_enabled,
					  __strlen_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strlen, 1,
					  __strlen_generic))

  IFUNC_IMPL (i, name, strcmp,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strcmp, rvv_enabled,
					  __strcmp_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strcmp, 1,
					  __strcmp_generic))

  IFUNC_IMPL (i, name, strncmp,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strncmp, rvv_enabled,
					  __strncmp_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strncmp, 1,
					  __strncmp_generic))

  IFUNC_IMPL (i, name, strncasecmp,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strncasecmp, rvv_enabled,
					  __strncasecmp_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strncasecmp, 1,
					  __strncasecmp_generic))

  IFUNC_IMPL (i, name, memccpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memccpy, rvv_enabled,
					  __memccpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memccpy, 1,
					  __memccpy_generic))

  IFUNC_IMPL (i, name, memcmp,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memcmp, rvv_enabled,
					  __memcmp_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memcmp, 1,
					  __memcmp_generic))

  IFUNC_IMPL (i, name, __memcmpeq,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, __memcmpeq, rvv_enabled,
					  ____memcmpeq_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, __memcmpeq, 1,
					  ____memcmpeq_generic))

  IFUNC_IMPL (i, name, rawmemchr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, rawmemchr, rvv_enabled,
					  __rawmemchr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, rawmemchr, 1,
					  __rawmemchr_generic))
  IFUNC_IMPL (i, name, memrchr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memrchr, rvv_enabled,
					  __memrchr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memrchr, 1,
					  __memrchr_generic))
  IFUNC_IMPL (i, name, stpncpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, stpncpy, rvv_enabled,
					  __stpncpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, stpncpy, 1,
					  __stpncpy_generic))

  IFUNC_IMPL (i, name, strcasecmp,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strcasecmp, rvv_enabled,
					  __strcasecmp_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strcasecmp, 1,
					  __strcasecmp_generic))
  IFUNC_IMPL (i, name, memchr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memchr, rvv_enabled,
					  __memchr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memchr, 1,
					  __memchr_generic))

  IFUNC_IMPL (i, name, strchr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strchr, rvv_enabled,
					  __strchr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strchr, 1,
					  __strchr_generic))

  IFUNC_IMPL (i, name, strstr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strstr, rvv_enabled,
					  __strstr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strstr, 1,
					  __strstr_generic))

  IFUNC_IMPL (i, name, strrchr,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strrchr, rvv_enabled,
					  __strrchr_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strrchr, 1,
					  __strrchr_generic))

  IFUNC_IMPL (i, name, memmove,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memmove,
					  rvv_enabled && spacemit_x60_tuned,
					  __memmove_spacemit_x60)
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memmove, rvv_enabled,
					  __memmove_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memmove, 1,
					  __memmove_generic))

  IFUNC_IMPL (i, name, strncat,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strncat, rvv_enabled,
					  __strncat_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strncat, 1,
					  __strncat_generic))

  IFUNC_IMPL (i, name, strncpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strncpy, rvv_enabled,
					  __strncpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strncpy, 1,
					  __strncpy_generic))

  IFUNC_IMPL (i, name, strnlen,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strnlen, rvv_enabled,
					  __strnlen_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strnlen, 1,
					  __strnlen_generic))

  IFUNC_IMPL (i, name, strpbrk,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strpbrk, rvv_enabled,
					  __strpbrk_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strpbrk, 1,
					  __strpbrk_generic))
  IFUNC_IMPL (i, name, strspn,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, strspn, rvv_enabled,
					  __strspn_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, strspn, 1,
					  __strspn_generic))
  return 0;
}
