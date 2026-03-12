ifeq ($(subdir),string)
sysdep_routines += \
  memcpy \
  memcpy-generic \
  memcpy-vector \
  memcpy_noalignment \
  memset \
  memset-generic \
  memset-vector \
  strcat \
  strcat-generic \
  strcat-vector \
  strcmp \
  strcmp-generic \
  strcmp-vector \
  strcpy \
  strcpy-generic \
  strcpy-spacemit-x60 \
  strcpy-vector \
  strlen \
  strlen-generic \
  strlen-vector \
  strncmp \
  strncmp-generic \
  strncmp-vector \
  # sysdep_routines

CFLAGS-memcpy_noalignment.c += -mno-strict-align
# Called during static initialization
CFLAGS-memset-generic.c += $(no-stack-protector)
CFLAGS-memcpy-generic.c += $(no-stack-protector)
endif
