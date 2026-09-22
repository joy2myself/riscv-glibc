/* Indexed coefficient tables for the double-precision vector lgamma kernel.
   Copyright (C) 2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either version 2.1
   of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  */

#include <vecmath_config.h>

const struct v_lgamma_data __v_lgamma_data = {
  .neg_bnd_to_2_c = { { 0x1.83fe966af535fp+0, -0x1.775909a36a68cp-55 },
        { 0x1.36eebb002f55dp-1, -0x1.8d4b2124a39f8p-55 },
        { 0x1.694a6058a7858p-6, -0x1.1d8c8b9b4d80ep-61 },
        { 0x1.1718d7ca09e5bp-6, 0x1.83195b07fd25p-60 },
        { 0x1.7339fe04b2764p-10, -0x1.48648f4a4bf9ep-64 },
        { 0x1.8d32f682aa0bdp-11, -0x1.90953dadfba01p-66 },
        { 0x1.809f04ee6e0fap-14, -0x1.6100a3d177e7cp-68 },
        { 0x1.48eaa81657361p-15, 0x1.42b86f83f623dp-70 },
        { 0x1.9297adb2def5ap-18, -0x1.0c295492288fdp-73 },
        { 0x1.286fb8cbaebb5p-19, -0x1.ebf1f6a584b62p-77 },
        { 0x1.a92e0a5de4bc9p-22, -0x1.688db240a9a82p-76 },
        { 0x1.1a9d4d8c6284ep-23, 0x1.08624c7186639p-78 },
        { 0x1.c4cd2594e91c9p-26, 0x1.a225e59cfef5dp-81 },
        { 0x1.18737ec8e68aap-27, 0x1.6421ab9e4c553p-82 },
        { 0x1.e6028795be5c1p-30, -0x1.b4bd2f36bb0ddp-84 },
        { 0x1.1eacaeb800afdp-31, -0x1.cedc75cbd8cc4p-87 },
        { 0x1.06bce9e1f6b9p-33, -0x1.0fb524fda64ebp-87 },
        { 0x1.2bb110a516b79p-35, 0x1.2d7a224941f92p-90 },
        { 0x1.1e0024589b848p-37, 0x1.1a84fa416703dp-92 },
        { 0x1.3ec6ebb4ba73p-39, -0x1.eb9cd73c631d9p-94 },
        { 0x1.3936af34a015p-41, -0x1.83b442df73674p-95 },
        { 0x1.57d4ece262198p-43, 0x1.20b160178fca3p-97 },
        { 0x1.5b530d802ffffp-45, 0x1.8bb2ee582dbd1p-104 },
        { 0x1.7c323d20053dp-47, -0x1.32d5dd28cac13p-103 },
        { 0x1.6131339e2b76ap-49, 0x1.4422b26549563p-104 },
        { 0x1.b4f6aaa0d8886p-52, 0x1.796e66ee63a34p-106 } },
  .neg_3_to_bnd_c = { { -0x1.ea12da904b18cp+0, -0x1.220130f99b2cbp-54 },
        { 0x1.3267f3c265a52p-1, -0x1.1c630ff19db83p-55 },
        { -0x1.4185ac30c8bf2p-4, 0x1.f161263693e12p-59 },
        { 0x1.f504accc9f19bp-7, -0x1.eacc0226c7208p-62 },
        { -0x1.8588458207eacp-9, 0x1.4b51668f3dff4p-63 },
        { 0x1.4373f7cc709b3p-11, -0x1.2474e20e777aep-66 },
        { -0x1.12239bdd6c013p-13, 0x1.46df69aa3032cp-69 },
        { 0x1.dba65e27421c4p-16, 0x1.36bfe12004625p-71 },
        { -0x1.a2d2504d7e987p-18, 0x1.a6b6dfe9fa6f4p-75 },
        { 0x1.7581739ee6087p-20, -0x1.23100aff1ab78p-77 },
        { -0x1.506c65fad6188p-22, -0x1.17d1749eb738fp-77 },
        { 0x1.318ef724f7814p-24, -0x1.6517edf27abc2p-82 },
        { -0x1.17767260d9825p-26, -0x1.03e683dfe1d87p-82 },
        { 0x1.011e34454ade3p-28, -0x1.2a844fd6e7622p-82 },
        { -0x1.db8b9e6b1e0e1p-31, -0x1.f73094d7d4e4p-85 },
        { 0x1.b9bab1fca321p-33, -0x1.89b7151db3448p-88 },
        { -0x1.9bed46f4d0aa2p-35, -0x1.1af984f4b8e4p-90 },
        { 0x1.81780d4242119p-37, -0x1.7576147eed5dcp-91 },
        { -0x1.69d3ce15a3f87p-39, 0x1.9823fc1b6be9ep-93 },
        { 0x1.5494a34ce847p-41, -0x1.0ec96907a37b4p-95 },
        { -0x1.4160eccd5982ep-43, 0x1.2b4a92d7b0c76p-98 },
        { 0x1.2fe0a640fb1ap-45, -0x1.b1f41cf9689edp-106 },
        { -0x1.1ff8659402df5p-47, 0x1.511b4097b0fccp-104 },
        { 0x1.1357a3e9e1cf1p-49, -0x1.bda0e37acdd97p-106 },
        { -0x1.0b145c5c5b9abp-51, -0x1.b52b4002dbc0cp-105 },
        { 0x1.dc1888c7036cap-54, 0x1.2a69b093c9a8bp-109 },
        { -0x1.060b52b5d6f68p-56, 0x1.d17de397a2b1p-110 } },
};
