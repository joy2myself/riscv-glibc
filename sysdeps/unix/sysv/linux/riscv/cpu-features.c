/* Initialize CPU feature data.
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

#include <cpu-features.h>
#include <sys/hwprobe.h>
#include <elf/dl-hwcaps.h>
#include <elf/dl-tunables.h>
#include <ldsodefs.h>
#include <profile-detection.h>

extern void TUNABLE_CALLBACK (set_hwcaps) (tunable_val_t *) attribute_hidden;

static inline uint32_t
riscv_get_tune_flags (const struct cpu_features *cpu_features)
{
  uint32_t tune_flags = RISCV_CPU_TUNE_NONE;

#ifndef RISCV_DISABLE_SPACEMIT_X60_TUNING
  /* SpaceMiT X60 tuning currently applies to all known revisions of this
     marchid, so mimpid is intentionally not part of the match.  */
  if (cpu_features->marchid == RISCV_MARCHID_SPACEMIT_X60)
    tune_flags |= RISCV_CPU_TUNE_SPACEMIT_X60;
#endif

  return tune_flags;
}

static inline void
init_cpu_features (struct cpu_features *cpu_features)
{
  struct riscv_hwprobe pairs[5];

  cpu_features->zicboz_block_size = 0;
  cpu_features->zicbop_block_size = 0;
  cpu_features->profile_level = 0;
  cpu_features->tune_flags = RISCV_CPU_TUNE_NONE;

  pairs[0].key = RISCV_HWPROBE_KEY_MVENDORID;
  pairs[1].key = RISCV_HWPROBE_KEY_MARCHID;
  pairs[2].key = RISCV_HWPROBE_KEY_MIMPID;
  pairs[3].key = RISCV_HWPROBE_KEY_IMA_EXT_0;
  pairs[4].key = RISCV_HWPROBE_KEY_CPUPERF_0;

  if (__riscv_hwprobe (pairs, 5, 0, NULL, 0) == 0)
    {
      cpu_features->mvendorid = pairs[0].value;
      cpu_features->marchid = pairs[1].value;
      cpu_features->mimpid = pairs[2].value;
      cpu_features->extensions = pairs[3].value;
      cpu_features->cpuperf = pairs[4].value;
    }
  else
    {
      cpu_features->mvendorid = 0;
      cpu_features->marchid = 0;
      cpu_features->mimpid = 0;
      cpu_features->extensions = 0;
      cpu_features->cpuperf = 0;
    }

  cpu_features->tune_flags = riscv_get_tune_flags (cpu_features);

  TUNABLE_GET (glibc, cpu, hwcaps, tunable_val_t *,
	       TUNABLE_CALLBACK (set_hwcaps));

  if (cpu_features->extensions & RISCV_HWPROBE_EXT_ZICBOP)
    cpu_features->zicbop_block_size =
      riscv_hwprobe_get_block_size (RISCV_HWPROBE_KEY_ZICBOP_BLOCK_SIZE);

  if (cpu_features->extensions & RISCV_HWPROBE_EXT_ZICBOZ)
    cpu_features->zicboz_block_size =
      riscv_hwprobe_get_block_size (RISCV_HWPROBE_KEY_ZICBOZ_BLOCK_SIZE);

  cpu_features->profile_level =
    riscv_get_profile_level_from_hwprobe (cpu_features->extensions,
					  cpu_features->cpuperf,
					  cpu_features->zicbop_block_size,
					  cpu_features->zicboz_block_size);
}
