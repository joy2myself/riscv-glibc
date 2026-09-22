/* Tables and constants for the single-precision vector powf function.
   The values are from glibc's scalar powf log2/exp2 source family.

   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <vecmath_config.h>

const struct v_powf_data __v_powf_data = {
  .log2_invc = {
    0x1.661ec79f8f3bep+0, 0x1.571ed4aaf883dp+0,
    0x1.49539f0f010bp+0, 0x1.3c995b0b80385p+0,
    0x1.30d190c8864a5p+0, 0x1.25e227b0b8eap+0,
    0x1.1bb4a4a1a343fp+0, 0x1.12358f08ae5bap+0,
    0x1.0953f419900a7p+0, 0x1p+0,
    0x1.e608cfd9a47acp-1, 0x1.ca4b31f026aap-1,
    0x1.b2036576afce6p-1, 0x1.9c2d163a1aa2dp-1,
    0x1.886e6037841edp-1, 0x1.767dcf5534862p-1,
  },
  .log2_logc = {
    -0x1.efec65b963019p-2, -0x1.b0b6832d4fca4p-2,
    -0x1.7418b0a1fb77bp-2, -0x1.39de91a6dcf7bp-2,
    -0x1.01d9bf3f2b631p-2, -0x1.97c1d1b3b7afp-3,
    -0x1.2f9e393af3c9fp-3, -0x1.960cbbf788d5cp-4,
    -0x1.a6f9db6475fcep-5, 0x0p+0,
    0x1.338ca9f24f53dp-4, 0x1.476a9543891bap-3,
    0x1.e840b4ac4e4d2p-3, 0x1.40645f0c6651cp-2,
    0x1.88e9c2c1b9ff8p-2, 0x1.ce0a44eb17bccp-2,
  },
  .exp2_tab = {
    0x3ff0000000000000ULL, 0x3fefd9b0d3158574ULL,
    0x3fefb5586cf9890fULL, 0x3fef9301d0125b51ULL,
    0x3fef72b83c7d517bULL, 0x3fef54873168b9aaULL,
    0x3fef387a6e756238ULL, 0x3fef1e9df51fdee1ULL,
    0x3fef06fe0a31b715ULL, 0x3feef1a7373aa9cbULL,
    0x3feedea64c123422ULL, 0x3feece086061892dULL,
    0x3feebfdad5362a27ULL, 0x3feeb42b569d4f82ULL,
    0x3feeab07dd485429ULL, 0x3feea47eb03a5585ULL,
    0x3feea09e667f3bcdULL, 0x3fee9f75e8ec5f74ULL,
    0x3feea11473eb0187ULL, 0x3feea589994cce13ULL,
    0x3feeace5422aa0dbULL, 0x3feeb737b0cdc5e5ULL,
    0x3feec49182a3f090ULL, 0x3feed503b23e255dULL,
    0x3feee89f995ad3adULL, 0x3feeff76f2fb5e47ULL,
    0x3fef199bdd85529cULL, 0x3fef3720dcef9069ULL,
    0x3fef5818dcfba487ULL, 0x3fef7c97337b9b5fULL,
    0x3fefa4afa2a490daULL, 0x3fefd0765b6e4540ULL,
  },
};
