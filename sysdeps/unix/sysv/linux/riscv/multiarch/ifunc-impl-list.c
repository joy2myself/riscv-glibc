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
    (cpu_features->tune_flags & RISCV_CPU_TUNE_SPACEMIT_X60) != 0;

  IFUNC_IMPL (i, name, memcpy,
	      RISCV_IFUNC_IMPL_ADD_RVA23 (array, i, memcpy, rvv_enabled,
					  __memcpy_vector)
	      RISCV_IFUNC_IMPL_ADD_RVA22 (array, i, memcpy, fast_unaligned,
					  __memcpy_noalignment)
	      RISCV_IFUNC_IMPL_ADD_INIT (array, i, memcpy, 1,
					 __memcpy_generic))

  IFUNC_IMPL (i, name, memset,
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

  return 0;
}
