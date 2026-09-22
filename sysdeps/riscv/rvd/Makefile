ifeq ($(subdir),mathvec)

libmvec-sysdep_routines = \
	$(addprefix $(libmvec-vprefix)v_d_, $(libmvec-double-funcs)) \
	$(addsuffix _data, $(addprefix v_d_, $(libmvec-double-funcs-with-data))) \
	$(addprefix v_f_, $(libmvec-float-funcs)) \
	$(addsuffix _data, $(addprefix v_f_, $(libmvec-float-funcs-with-data))) \

define riscv64-vector-cflags-template
CFLAGS-$(1).c += -march=rv64gcv -Wno-maybe-uninitialized -Wno-undef
endef

define riscv32-vector-cflags-template
CFLAGS-$(1).c += -march=rv32gcv -Wno-maybe-uninitialized -Wno-undef
endef

ifeq ($(config-machine),riscv64)
$(foreach f,$(libmvec-sysdep_routines), $(eval $(call riscv64-vector-cflags-template,$(f))))
endif

ifeq ($(config-machine),riscv32)
$(foreach f,$(libmvec-sysdep_routines), $(eval $(call riscv32-vector-cflags-template,$(f))))
endif

endif

ifneq ($(multi-arch),no)
ifeq ($(subdir)$(build-mathvec),mathyes)

ifneq ($(have-test-clang),yes)

tests += \
  $(libmvec-abi-func-tests)     \

endif

endif

endif
