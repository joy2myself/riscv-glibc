/* RISC-V architecture probe interface
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
   License along with the GNU C Library.  If not, see
   <https://www.gnu.org/licenses/>.  */

#ifndef _SYS_HWPROBE_H
#define _SYS_HWPROBE_H 1

#include <features.h>
#include <sched.h>
#include <stddef.h>
#include <errno.h>
#ifdef __has_include
# if __has_include (<asm/hwprobe.h>)
#  include <asm/hwprobe.h>
# endif
#endif

/* Define a (probably stale) version of the interface if the Linux headers
   aren't present.  */
#ifndef RISCV_HWPROBE_KEY_MVENDORID
struct riscv_hwprobe {
	signed long long int key;
	unsigned long long int value;
};

#define RISCV_HWPROBE_KEY_MVENDORID 0
#define RISCV_HWPROBE_KEY_MARCHID 1
#define RISCV_HWPROBE_KEY_MIMPID 2
#define RISCV_HWPROBE_KEY_BASE_BEHAVIOR 3
#define  RISCV_HWPROBE_BASE_BEHAVIOR_IMA (1 << 0)
#define RISCV_HWPROBE_KEY_IMA_EXT_0 4
#define  RISCV_HWPROBE_IMA_FD (1 << 0)
#define  RISCV_HWPROBE_IMA_C (1 << 1)
#define  RISCV_HWPROBE_IMA_V (1 << 2)
#define  RISCV_HWPROBE_EXT_ZBA (1 << 3)
#define  RISCV_HWPROBE_EXT_ZBB (1 << 4)
#define  RISCV_HWPROBE_EXT_ZBS (1 << 5)
#define  RISCV_HWPROBE_EXT_ZICBOZ (1 << 6)
#define RISCV_HWPROBE_KEY_CPUPERF_0 5
#define  RISCV_HWPROBE_MISALIGNED_UNKNOWN (0 << 0)
#define  RISCV_HWPROBE_MISALIGNED_EMULATED (1 << 0)
#define  RISCV_HWPROBE_MISALIGNED_SLOW (2 << 0)
#define  RISCV_HWPROBE_MISALIGNED_FAST (3 << 0)
#define  RISCV_HWPROBE_MISALIGNED_UNSUPPORTED (4 << 0)
#define  RISCV_HWPROBE_MISALIGNED_MASK (7 << 0)
#define RISCV_HWPROBE_KEY_ZICBOZ_BLOCK_SIZE 6

#endif /* RISCV_HWPROBE_KEY_MVENDORID */

#ifndef RISCV_HWPROBE_KEY_ZICBOZ_BLOCK_SIZE
#define RISCV_HWPROBE_KEY_ZICBOZ_BLOCK_SIZE 6
#endif

#ifndef RISCV_HWPROBE_EXT_ZICBOZ
#define RISCV_HWPROBE_EXT_ZICBOZ (1 << 6)
#endif

#ifndef RISCV_HWPROBE_EXT_ZBC
#define RISCV_HWPROBE_EXT_ZBC (1 << 7)
#endif

#ifndef RISCV_HWPROBE_EXT_ZBKB
#define RISCV_HWPROBE_EXT_ZBKB (1 << 8)
#endif

#ifndef RISCV_HWPROBE_EXT_ZBKC
#define RISCV_HWPROBE_EXT_ZBKC (1 << 9)
#endif

#ifndef RISCV_HWPROBE_EXT_ZBKX
#define RISCV_HWPROBE_EXT_ZBKX (1 << 10)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKND
#define RISCV_HWPROBE_EXT_ZKND (1 << 11)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKNE
#define RISCV_HWPROBE_EXT_ZKNE (1 << 12)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKNH
#define RISCV_HWPROBE_EXT_ZKNH (1 << 13)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKSED
#define RISCV_HWPROBE_EXT_ZKSED (1 << 14)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKSH
#define RISCV_HWPROBE_EXT_ZKSH (1 << 15)
#endif

#ifndef RISCV_HWPROBE_EXT_ZKT
#define RISCV_HWPROBE_EXT_ZKT (1 << 16)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVBB
#define RISCV_HWPROBE_EXT_ZVBB (1 << 17)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVBC
#define RISCV_HWPROBE_EXT_ZVBC (1 << 18)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKB
#define RISCV_HWPROBE_EXT_ZVKB (1 << 19)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKG
#define RISCV_HWPROBE_EXT_ZVKG (1 << 20)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKNED
#define RISCV_HWPROBE_EXT_ZVKNED (1 << 21)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKNHA
#define RISCV_HWPROBE_EXT_ZVKNHA (1 << 22)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKNHB
#define RISCV_HWPROBE_EXT_ZVKNHB (1 << 23)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKSED
#define RISCV_HWPROBE_EXT_ZVKSED (1 << 24)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKSH
#define RISCV_HWPROBE_EXT_ZVKSH (1 << 25)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVKT
#define RISCV_HWPROBE_EXT_ZVKT (1 << 26)
#endif

#ifndef RISCV_HWPROBE_EXT_ZFH
#define RISCV_HWPROBE_EXT_ZFH (1 << 27)
#endif

#ifndef RISCV_HWPROBE_EXT_ZFHMIN
#define RISCV_HWPROBE_EXT_ZFHMIN (1 << 28)
#endif

#ifndef RISCV_HWPROBE_EXT_ZIHINTNTL
#define RISCV_HWPROBE_EXT_ZIHINTNTL (1 << 29)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVFH
#define RISCV_HWPROBE_EXT_ZVFH (1 << 30)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVFHMIN
#define RISCV_HWPROBE_EXT_ZVFHMIN (1ULL << 31)
#endif

#ifndef RISCV_HWPROBE_EXT_ZFA
#define RISCV_HWPROBE_EXT_ZFA (1ULL << 32)
#endif

#ifndef RISCV_HWPROBE_EXT_ZTSO
#define RISCV_HWPROBE_EXT_ZTSO (1ULL << 33)
#endif

#ifndef RISCV_HWPROBE_EXT_ZACAS
#define RISCV_HWPROBE_EXT_ZACAS (1ULL << 34)
#endif

#ifndef RISCV_HWPROBE_EXT_ZICOND
#define RISCV_HWPROBE_EXT_ZICOND (1ULL << 35)
#endif

#ifndef RISCV_HWPROBE_EXT_ZIHINTPAUSE
#define RISCV_HWPROBE_EXT_ZIHINTPAUSE (1ULL << 36)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVE32X
#define RISCV_HWPROBE_EXT_ZVE32X (1ULL << 37)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVE32F
#define RISCV_HWPROBE_EXT_ZVE32F (1ULL << 38)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVE64X
#define RISCV_HWPROBE_EXT_ZVE64X (1ULL << 39)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVE64F
#define RISCV_HWPROBE_EXT_ZVE64F (1ULL << 40)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVE64D
#define RISCV_HWPROBE_EXT_ZVE64D (1ULL << 41)
#endif

#ifndef RISCV_HWPROBE_EXT_ZIMOP
#define RISCV_HWPROBE_EXT_ZIMOP (1ULL << 42)
#endif

#ifndef RISCV_HWPROBE_EXT_ZCA
#define RISCV_HWPROBE_EXT_ZCA (1ULL << 43)
#endif

#ifndef RISCV_HWPROBE_EXT_ZCB
#define RISCV_HWPROBE_EXT_ZCB (1ULL << 44)
#endif

#ifndef RISCV_HWPROBE_EXT_ZCD
#define RISCV_HWPROBE_EXT_ZCD (1ULL << 45)
#endif

#ifndef RISCV_HWPROBE_EXT_ZCF
#define RISCV_HWPROBE_EXT_ZCF (1ULL << 46)
#endif

#ifndef RISCV_HWPROBE_EXT_ZCMOP
#define RISCV_HWPROBE_EXT_ZCMOP (1ULL << 47)
#endif

#ifndef RISCV_HWPROBE_EXT_ZAWRS
#define RISCV_HWPROBE_EXT_ZAWRS (1ULL << 48)
#endif

#ifndef RISCV_HWPROBE_EXT_SUPM
#define RISCV_HWPROBE_EXT_SUPM (1ULL << 49)
#endif

#ifndef RISCV_HWPROBE_EXT_ZICNTR
#define RISCV_HWPROBE_EXT_ZICNTR (1ULL << 50)
#endif

#ifndef RISCV_HWPROBE_EXT_ZIHPM
#define RISCV_HWPROBE_EXT_ZIHPM (1ULL << 51)
#endif

#ifndef RISCV_HWPROBE_EXT_ZFBFMIN
#define RISCV_HWPROBE_EXT_ZFBFMIN (1ULL << 52)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVFBFMIN
#define RISCV_HWPROBE_EXT_ZVFBFMIN (1ULL << 53)
#endif

#ifndef RISCV_HWPROBE_EXT_ZVFBFWMA
#define RISCV_HWPROBE_EXT_ZVFBFWMA (1ULL << 54)
#endif

#ifndef RISCV_HWPROBE_EXT_ZICBOM
#define RISCV_HWPROBE_EXT_ZICBOM (1ULL << 55)
#endif

#ifndef RISCV_HWPROBE_EXT_ZAAMO
#define RISCV_HWPROBE_EXT_ZAAMO (1ULL << 56)
#endif

#ifndef RISCV_HWPROBE_EXT_ZALRSC
#define RISCV_HWPROBE_EXT_ZALRSC (1ULL << 57)
#endif

#ifndef RISCV_HWPROBE_EXT_ZABHA
#define RISCV_HWPROBE_EXT_ZABHA (1ULL << 58)
#endif

#ifndef RISCV_HWPROBE_EXT_ZICBOP
#define RISCV_HWPROBE_EXT_ZICBOP (1ULL << 59)
#endif

#ifndef RISCV_HWPROBE_KEY_ZICBOP_BLOCK_SIZE
#define RISCV_HWPROBE_KEY_ZICBOP_BLOCK_SIZE 15
#endif

__BEGIN_DECLS

#if defined __cplusplus || !__GNUC_PREREQ (2, 7)
# define __RISCV_HWPROBE_CPUS_TYPE cpu_set_t *
#else
/* The fourth argument to __riscv_hwprobe should be a null pointer or a
   pointer to a cpu_set_t (either the fixed-size type or allocated with
   CPU_ALLOC).  However, early versions of this header file used the
   argument type unsigned long int *.  The transparent union allows
   the argument to be either cpu_set_t * or unsigned long int * for
   compatibility.  The older header file requiring unsigned long int *
   can be identified by the lack of the __RISCV_HWPROBE_CPUS_TYPE macro.
   In C++ and with compilers that do not support transparent unions, the
   argument type must be cpu_set_t *.  */
typedef union {
	cpu_set_t *__cs;
	unsigned long int *__ul;
} __RISCV_HWPROBE_CPUS_TYPE __attribute__ ((__transparent_union__));
# define __RISCV_HWPROBE_CPUS_TYPE __RISCV_HWPROBE_CPUS_TYPE
#endif

extern int __riscv_hwprobe (struct riscv_hwprobe *__pairs,
			    size_t __pair_count, size_t __cpusetsize,
			    __RISCV_HWPROBE_CPUS_TYPE __cpus,
			    unsigned int __flags)
     __THROW __nonnull ((1)) __attr_access ((__read_write__, 1, 2));

/* A pointer to the __riscv_hwprobe function is passed as the second
   argument to ifunc selector routines. Include a function pointer type for
   convenience in calling the function in those settings. */
typedef int (*__riscv_hwprobe_t) (struct riscv_hwprobe *__pairs,
				  size_t __pair_count, size_t __cpusetsize,
				  __RISCV_HWPROBE_CPUS_TYPE __cpus,
				  unsigned int __flags)
     __nonnull ((1)) __attr_access ((__read_write__, 1, 2));

/* Helper function usable from ifunc selectors that probes a single key. */
static __inline int
__riscv_hwprobe_one(__riscv_hwprobe_t hwprobe_func,
                    long long int key,
                    unsigned long long int *value)
{
  struct riscv_hwprobe pair;
  int rc;

  /* Earlier versions of glibc pass NULL as the second ifunc parameter. Other C
     libraries on non-Linux systems may pass +1 as this function pointer to
     indicate no support. Users copying this function to exotic worlds
     (non-Linux non-glibc) may want to do additional validity checks here. */
  if (hwprobe_func == NULL)
    return ENOSYS;

  pair.key = key;
  rc = hwprobe_func (&pair, 1, 0, NULL, 0);
  if (rc != 0)
    return rc;

  if (pair.key < 0)
    return ENOENT;

  *value = pair.value;
  return 0;
}

__END_DECLS

#endif /* sys/hwprobe.h */
