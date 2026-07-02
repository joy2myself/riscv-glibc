/* RISC-V profile detection helpers.
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

#ifndef _RISCV_PROFILE_DETECTION_H
#define _RISCV_PROFILE_DETECTION_H

#include <stdbool.h>
#include <stdint.h>
#include <sys/hwprobe.h>

#define RISCV_PROFILE_ZIC64B_MIN_BLOCK_SIZE 64

#define RISCV_PROFILE_MASK_RVA20U64 \
  (RISCV_HWPROBE_IMA_FD | RISCV_HWPROBE_IMA_C | RISCV_HWPROBE_EXT_ZICNTR)

#define RISCV_PROFILE_MASK_RVA22U64 \
  (RISCV_PROFILE_MASK_RVA20U64 | RISCV_HWPROBE_EXT_ZIHINTPAUSE \
   | RISCV_HWPROBE_EXT_ZBA | RISCV_HWPROBE_EXT_ZBB | RISCV_HWPROBE_EXT_ZBS \
   | RISCV_HWPROBE_EXT_ZICBOM | RISCV_HWPROBE_EXT_ZICBOP \
   | RISCV_HWPROBE_EXT_ZICBOZ | RISCV_HWPROBE_EXT_ZFHMIN \
   | RISCV_HWPROBE_EXT_ZKT | RISCV_HWPROBE_EXT_ZIHPM)

#define RISCV_PROFILE_MASK_RVA23U64 \
  (RISCV_PROFILE_MASK_RVA22U64 | RISCV_HWPROBE_IMA_V \
   | RISCV_HWPROBE_EXT_ZVFHMIN | RISCV_HWPROBE_EXT_ZVBB \
   | RISCV_HWPROBE_EXT_ZVKT | RISCV_HWPROBE_EXT_ZIHINTNTL \
   | RISCV_HWPROBE_EXT_ZICOND | RISCV_HWPROBE_EXT_ZIMOP \
   | RISCV_HWPROBE_EXT_ZCMOP | RISCV_HWPROBE_EXT_ZCB \
   | RISCV_HWPROBE_EXT_ZFA | RISCV_HWPROBE_EXT_ZAWRS)

/* Temporary debug masks for downward-trimmed RVA23 validation.  These match
   the subset currently supported by the SpaceMiT X60 lab board, after
   removing extensions that are missing in either the hardware ISA string or
   the local toolchain.  */
#ifdef RISCV_DEBUG_TRIMMED_RVA23
# define RISCV_PROFILE_MASK_TRIMMED_RVA22U64 \
  (RISCV_PROFILE_MASK_RVA20U64 | RISCV_HWPROBE_EXT_ZIHINTPAUSE \
   | RISCV_HWPROBE_EXT_ZBA | RISCV_HWPROBE_EXT_ZBB | RISCV_HWPROBE_EXT_ZBS \
   | RISCV_HWPROBE_EXT_ZICBOM | RISCV_HWPROBE_EXT_ZICBOZ \
   | RISCV_HWPROBE_EXT_ZFHMIN | RISCV_HWPROBE_EXT_ZKT \
   | RISCV_HWPROBE_EXT_ZIHPM)

# define RISCV_PROFILE_MASK_TRIMMED_RVA23U64 \
  (RISCV_PROFILE_MASK_TRIMMED_RVA22U64 | RISCV_HWPROBE_IMA_V \
   | RISCV_HWPROBE_EXT_ZVFHMIN | RISCV_HWPROBE_EXT_ZVKT \
   | RISCV_HWPROBE_EXT_ZICOND)
#endif

static inline unsigned int
riscv_hwprobe_get_block_size (long long int key)
{
  struct riscv_hwprobe pair = { .key = key };

  if (__riscv_hwprobe (&pair, 1, 0, NULL, 0) == 0 && pair.key == key)
    return pair.value;

  return 0;
}

static inline bool
riscv_hwprobe_have_zicclsm (uint64_t cpuperf)
{
  return ((cpuperf & RISCV_HWPROBE_MISALIGNED_MASK)
	  != RISCV_HWPROBE_MISALIGNED_UNSUPPORTED);
}

static inline bool
riscv_hwprobe_fast_unaligned (uint64_t cpuperf)
{
  return ((cpuperf & RISCV_HWPROBE_MISALIGNED_MASK)
	  == RISCV_HWPROBE_MISALIGNED_FAST);
}

static inline unsigned int
riscv_get_profile_level_from_hwprobe (uint64_t extensions, uint64_t cpuperf,
				      unsigned int zicbop_block_size,
				      unsigned int zicboz_block_size)
{
  if (!riscv_hwprobe_have_zicclsm (cpuperf))
    return 0;

  /* Extensions without explicit detection (yet):
     Ziccif, Ziccrse, Ziccamoa, Za128rs.  */
  if ((extensions & RISCV_PROFILE_MASK_RVA20U64)
      != RISCV_PROFILE_MASK_RVA20U64)
    return 0;

#ifdef RISCV_DEBUG_TRIMMED_RVA23
  if ((extensions & RISCV_PROFILE_MASK_TRIMMED_RVA22U64)
      != RISCV_PROFILE_MASK_TRIMMED_RVA22U64)
    return 20;

  if ((extensions & RISCV_PROFILE_MASK_TRIMMED_RVA23U64)
      != RISCV_PROFILE_MASK_TRIMMED_RVA23U64)
    return 22;

  return 23;
#else
  /* Extensions without explicit detection (yet): Zic64b.  */
  if ((extensions & RISCV_PROFILE_MASK_RVA22U64)
      != RISCV_PROFILE_MASK_RVA22U64
      || zicbop_block_size < RISCV_PROFILE_ZIC64B_MIN_BLOCK_SIZE
      || zicboz_block_size < RISCV_PROFILE_ZIC64B_MIN_BLOCK_SIZE)
    return 20;

  if ((extensions & RISCV_PROFILE_MASK_RVA23U64)
      != RISCV_PROFILE_MASK_RVA23U64)
    return 22;

  return 23;
#endif
}

#endif /* _RISCV_PROFILE_DETECTION_H */
