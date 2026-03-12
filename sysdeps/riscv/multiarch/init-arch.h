/* Initialize CPU feature data for RISC-V multiarch.
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

#include <ldsodefs.h>
#include <cpu-features.h>
#include <profile-detection.h>

#define INIT_ARCH()							\
  uint64_t __attribute__ ((unused)) mvendorid =			\
    GLRO(dl_riscv_cpu_features).mvendorid;				\
  uint64_t __attribute__ ((unused)) marchid =				\
    GLRO(dl_riscv_cpu_features).marchid;				\
  uint64_t __attribute__ ((unused)) mimpid =				\
    GLRO(dl_riscv_cpu_features).mimpid;					\
  uint64_t __attribute__ ((unused)) extensions =			\
    GLRO(dl_riscv_cpu_features).extensions;				\
  uint64_t __attribute__ ((unused)) cpuperf =				\
    GLRO(dl_riscv_cpu_features).cpuperf;				\
  unsigned int __attribute__ ((unused)) zicboz_block_size =		\
    GLRO(dl_riscv_cpu_features).zicboz_block_size;			\
  unsigned int __attribute__ ((unused)) zicbop_block_size =		\
    GLRO(dl_riscv_cpu_features).zicbop_block_size;			\
  unsigned int __attribute__ ((unused)) profile_level =			\
    GLRO(dl_riscv_cpu_features).profile_level;				\
  uint32_t __attribute__ ((unused)) tune_flags =			\
    GLRO(dl_riscv_cpu_features).tune_flags;				\
  bool __attribute__ ((unused)) zicclsm =				\
    riscv_hwprobe_have_zicclsm (cpuperf);				\
  bool __attribute__ ((unused)) fast_unaligned =			\
    riscv_hwprobe_fast_unaligned (cpuperf);

/* Helper macros for common extension checks.  */
#define HAS_VECTOR()	(extensions & RISCV_HWPROBE_IMA_V)
#define HAS_ZBA()	(extensions & RISCV_HWPROBE_EXT_ZBA)
#define HAS_ZBB()	(extensions & RISCV_HWPROBE_EXT_ZBB)
#define HAS_ZBS()	(extensions & RISCV_HWPROBE_EXT_ZBS)
#define HAS_ZICBOZ()	(extensions & RISCV_HWPROBE_EXT_ZICBOZ)
#define HAS_ZVBB()	(extensions & RISCV_HWPROBE_EXT_ZVBB)
