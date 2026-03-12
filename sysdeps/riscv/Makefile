ifeq ($(subdir),misc)
sysdep_headers += sys/asm.h
endif

ifeq ($(subdir),resolv)
CFLAGS-ns_rr_cursor_next.c += -fno-tree-vectorize
CFLAGS-res-noaaaa.c += -fno-tree-vectorize
CFLAGS-res_send.c += -fno-tree-vectorize
CFLAGS-dns-host.c += -fno-tree-vectorize
CFLAGS-dns-network.c += -fno-tree-vectorize
endif

ifeq ($(subdir),nss)
CFLAGS-files-hosts.c += -fno-tree-vectorize
CFLAGS-getaddrinfo.c += -fno-tree-vectorize
CFLAGS-getnameinfo.c += -fno-tree-vectorize
endif

ifeq ($(subdir),elf)
gen-as-const-headers += dl-link.sym
sysdep-dl-routines += \
  dl-get-cpu-features \
  # sysdep-dl-routines

ifneq ($(have-tunables),no)
tests += tst-hwcap-tunables
tst-hwcap-tunables-ARGS = -- $(host-test-program-cmd)
endif

endif

# RISC-V's assembler also needs to know about PIC as it changes the definition
# of some assembler macros.
ASFLAGS-.os += $(pic-ccflag)

ifeq (no,$(riscv-r-align))
ASFLAGS-.os += -Wa,-mno-relax
ASFLAGS-.o += -Wa,-mno-relax
sysdep-CFLAGS += -mno-relax
endif
