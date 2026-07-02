/* Data structure for RISC-V CPU features.
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

#ifndef _CPU_FEATURES_RISCV_H
#define _CPU_FEATURES_RISCV_H

#include <stdint.h>
#include <stdbool.h>
#include <sys/hwprobe.h>

#define RISCV_VENDOR_SIFIVE    0x489
#define RISCV_VENDOR_SPACEMIT  0x710
#define RISCV_VENDOR_THEAD     0x5b7

#define RISCV_MARCHID_SPACEMIT_X60 0x8000000058000001ULL
#define RISCV_MARCHID_SPACEMIT_C908 RISCV_MARCHID_SPACEMIT_X60

enum riscv_cpu_tune
{
  RISCV_CPU_TUNE_NONE = 0,
  /* Use the SpaceMiT X60 tuned implementations.  */
  RISCV_CPU_TUNE_SPACEMIT_X60,
  /* Use the C908 tuned implementations.  */
  RISCV_CPU_TUNE_C908,
};

struct cpu_features
{
  uint64_t mvendorid;
  uint64_t marchid;
  uint64_t mimpid;
  uint64_t extensions;
  uint64_t cpuperf;
  unsigned int zicboz_block_size;
  unsigned int zicbop_block_size;
  unsigned int profile_level;
  uint32_t tune;
};

/* Get a pointer to the CPU features structure.  */
extern const struct cpu_features *_dl_riscv_get_cpu_features (void)
  __attribute__ ((pure));

#define __get_cpu_features() _dl_riscv_get_cpu_features ()

#if defined (_LIBC) && !IS_IN (nonlib)
# define _dl_riscv_get_cpu_features() (&GLRO (dl_riscv_cpu_features))
#endif

#endif /* _CPU_FEATURES_RISCV_H */
