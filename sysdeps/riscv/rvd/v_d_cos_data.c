/* Large-argument reduction table for the double-precision cos implementation.
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

const struct v_cos_data __v_cos_data = {
  .two_over_pi = {
    0x1.45f306dc9c882p+499,  0x1.4a7f09d5f47d4p+446,
    0x1.a6ee06db14ad0p+393,  -0x1.b0ef1bef806bcp+342,
    0x1.8eaf7aef1586cp+290,  0x1.c91b8e909374cp+238,
    -0x1.ff9b6d115f630p+184, 0x1.921cfe1deb1d0p+132,
    -0x1.3b5963045df74p+82,  0x1.7d4baed1213a8p+30,
    -0x1.8e3f652e82070p-22,  0x1.3991d63983530p-76,
    0x1.cfa4e422fc5e0p-127, -0x1.036be27003b40p-179,
    -0x1.0fd33f8086800p-239, -0x1.dce94beb25c20p-285,
    0x1.b4d9fb3c9f2c4p-334,  -0x1.922c2e7026588p-386,
    0x1.7fa8b5d49eeb0p-438,  0x1.faf97c5ecf41cp-490,
    0x1.cfbc529497538p-543,  -0x1.012813b81ca8cp-594,
    0x1.0ac06608df900p-649,  -0x1.251503cc10f7cp-698,
    -0x1.942f27895871cp-750, 0x1.615ee61b08660p-804,
    -0x1.99ea83ad7e5f0p-854, 0x1.1bffb1009ae60p-909,
  },
};
