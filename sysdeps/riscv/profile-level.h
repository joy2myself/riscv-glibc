/* Header defining the RISC-V Profile level
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   In addition to the permissions in the GNU Lesser General Public
   License, the Free Software Foundation gives you unlimited
   permission to link the compiled version of this file with other
   programs, and to distribute those programs without any restriction
   coming from the use of this file.  (The Lesser General Public
   License restrictions do apply in other respects; for example, they
   cover modification of the file, and distribution when not linked
   into another program.)

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#ifndef _PROFILE_LEVEL_H
#define _PROFILE_LEVEL_H

#define MINIMUM_RISCV_PROFILE_LEVEL 0

#if defined __riscv_i && defined __riscv_m && defined __riscv_a \
    && defined __riscv_f && defined __riscv_d && defined __riscv_c
# define __RISCV_GC 1
#else
# define __RISCV_GC 0
#endif

#if __RISCV_GC && defined __riscv_zicsr && defined __riscv_zicntr \
    && defined __riscv_ziccif && defined __riscv_ziccrse \
    && defined __riscv_ziccamoa && defined __riscv_zicclsm \
    && (defined __riscv_za128rs || defined __riscv_za64rs)
# define __RISCV_PROFILE_RVA20 1
# undef MINIMUM_RISCV_PROFILE_LEVEL
# define MINIMUM_RISCV_PROFILE_LEVEL 20
#else
# define __RISCV_PROFILE_RVA20 0
#endif

#if __RISCV_PROFILE_RVA20 && defined __riscv_zic64b \
    && defined __riscv_za64rs && defined __riscv_zihintpause \
    && defined __riscv_zba && defined __riscv_zbb \
    && defined __riscv_zbs && defined __riscv_zihpm \
    && defined __riscv_zicbom && defined __riscv_zicbop \
    && defined __riscv_zicboz && defined __riscv_zfhmin \
    && defined __riscv_zkt
# define __RISCV_PROFILE_RVA22 1
# undef MINIMUM_RISCV_PROFILE_LEVEL
# define MINIMUM_RISCV_PROFILE_LEVEL 22
#else
# define __RISCV_PROFILE_RVA22 0
#endif

#if __RISCV_PROFILE_RVA22 && defined __riscv_v \
    && defined __riscv_zvfhmin && defined __riscv_zvbb \
    && defined __riscv_zvkt && defined __riscv_zihintntl \
    && defined __riscv_zicond && defined __riscv_zimop \
    && defined __riscv_zcmop && defined __riscv_zcb \
    && defined __riscv_zfa && defined __riscv_zawrs
# define __RISCV_PROFILE_RVA23 1
# undef MINIMUM_RISCV_PROFILE_LEVEL
# define MINIMUM_RISCV_PROFILE_LEVEL 23
#else
# define __RISCV_PROFILE_RVA23 0
#endif

/* Temporary debug mode for validating downward-trimmed RVA23 builds on
   systems whose hardware/toolchain cannot expose the full profile baseline.
   This must only be used together with a baseline -march that has already
   removed unsupported instructions, otherwise the resulting binaries can
   still trap with SIGILL at startup.  */
#ifdef RISCV_DEBUG_TRIMMED_RVA23
# undef __RISCV_PROFILE_RVA20
# undef __RISCV_PROFILE_RVA22
# undef __RISCV_PROFILE_RVA23
# undef MINIMUM_RISCV_PROFILE_LEVEL
# define __RISCV_PROFILE_RVA20 1
# define __RISCV_PROFILE_RVA22 1
# define __RISCV_PROFILE_RVA23 1
# define MINIMUM_RISCV_PROFILE_LEVEL 23
#endif

/* Profile level >= 23 guaranteed includes.  */
#define V_RISCV_PROFILE_LEVEL 23
#define ZVFHMIN_RISCV_PROFILE_LEVEL 23
#define ZVBB_RISCV_PROFILE_LEVEL 23
#define ZVKT_RISCV_PROFILE_LEVEL 23
#define ZIHINTNTL_RISCV_PROFILE_LEVEL 23
#define ZICOND_RISCV_PROFILE_LEVEL 23
#define ZIMOP_RISCV_PROFILE_LEVEL 23
#define ZCMOP_RISCV_PROFILE_LEVEL 23
#define ZCB_RISCV_PROFILE_LEVEL 23
#define ZFA_RISCV_PROFILE_LEVEL 23
#define ZAWRS_RISCV_PROFILE_LEVEL 23

/* Profile level >= 22 guaranteed includes.  */
#define HPM_RISCV_PROFILE_LEVEL 22
#define ZA64RS_RISCV_PROFILE_LEVEL 22
#define ZIHINTPAUSE_RISCV_PROFILE_LEVEL 22
#define ZBA_RISCV_PROFILE_LEVEL 22
#define ZBB_RISCV_PROFILE_LEVEL 22
#define ZBS_RISCV_PROFILE_LEVEL 22
#define ZIC64B_RISCV_PROFILE_LEVEL 22
#define ZICBOM_RISCV_PROFILE_LEVEL 22
#define ZICBOP_RISCV_PROFILE_LEVEL 22
#define ZICBOZ_RISCV_PROFILE_LEVEL 22
#define ZFHMIN_RISCV_PROFILE_LEVEL 22
#define ZKT_RISCV_PROFILE_LEVEL 22

/* Profile level >= 20 guaranteed includes.  */
#define FD_RISCV_PROFILE_LEVEL 20
#define A_RISCV_PROFILE_LEVEL 20
#define C_RISCV_PROFILE_LEVEL 20
#define ZICSR_RISCV_PROFILE_LEVEL 20
#define ZICNTR_RISCV_PROFILE_LEVEL 20
#define ZICCIF_RISCV_PROFILE_LEVEL 20
#define ZICCRSE_RISCV_PROFILE_LEVEL 20
#define ZICCAMOA_RISCV_PROFILE_LEVEL 20
#define ZA128RS_RISCV_PROFILE_LEVEL 20
#define ZICCLSM_RISCV_PROFILE_LEVEL 20

#define INIT_RISCV_PROFILE_LEVEL 0
#define RVA20_RISCV_PROFILE_LEVEL 20
#define RVA22_RISCV_PROFILE_LEVEL 22
#define RVA23_RISCV_PROFILE_LEVEL 23

#define PROFILE_SHOULD_BUILD(profile) \
  (IS_IN (libc) && MINIMUM_RISCV_PROFILE_LEVEL <= profile##_RISCV_PROFILE_LEVEL)

#endif
