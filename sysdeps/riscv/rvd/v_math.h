/* Utilities for RVV libmvec routines.
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

#ifndef _V_MATH_H
#define _V_MATH_H

#include <riscv_vector.h>
#include <math_private.h>
#include <assert.h>

#include <v_math_names.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#define _NEED_UNDEF_GNU_SOURCE_V_MATH
#endif

#include <math.h>
#ifndef NAN
_Static_assert (0, "NaN not available on this architecture");
#endif

#define __PASTE2_BASE(A, B) A##B
#define __PASTE2(A, B) __PASTE2_BASE (A, B)
#define __PASTE3_BASE(A, B, C) A##B##C
#define __PASTE3(A, B, C) __PASTE3_BASE (A, B, C)
#define __PASTE4_BASE(A, B, C, D) A##B##C##D
#define __PASTE4(A, B, C, D) __PASTE4_BASE (A, B, C, D)
#define __PASTE5_BASE(A, B, C, D, E) A##B##C##D##E
#define __PASTE5(A, B, C, D, E) __PASTE5_BASE (A, B, C, D, E)

/* The source-shaped log2f comparator has a separate importer entry.  */
#ifndef RVV_MATH_IS_LOG2F_ORIGINAL_log2f_original
# define RVV_MATH_IS_LOG2F_ORIGINAL_log2f_original 1
# define RVV_MATH_IS_LOG2F_ORIGINAL_INNER(name) \
  RVV_MATH_IS_LOG2F_ORIGINAL_##name
# define RVV_MATH_IS_LOG2F_ORIGINAL(name) \
  RVV_MATH_IS_LOG2F_ORIGINAL_INNER (name)
#endif

#ifndef RVV_MATH_IS_TGAMMAF_FAMILY_tgammaf
# define RVV_MATH_IS_TGAMMAF_FAMILY_tgammaf 1
# define RVV_MATH_IS_TGAMMAF_FAMILY_tgammaf_original 1
# define RVV_MATH_IS_TGAMMAF_FAMILY_lgammaf 1
# define RVV_MATH_IS_TGAMMAF_FAMILY_lgammaf_original 1
# define RVV_MATH_IS_TGAMMAF_FAMILY_INNER(name) \
  RVV_MATH_IS_TGAMMAF_FAMILY_##name
# define RVV_MATH_IS_TGAMMAF_FAMILY(name) \
  RVV_MATH_IS_TGAMMAF_FAMILY_INNER (name)
#endif

#ifndef RVV_MATH_IS_LOG2F_FAMILY_log2f
# define RVV_MATH_IS_LOG2F_FAMILY_log2f 1
# define RVV_MATH_IS_LOG2F_FAMILY_log2f_clean_assert 1
# define RVV_MATH_IS_LOG2F_FAMILY_log2f_original 1
# define RVV_MATH_IS_LOG2F_FAMILY_log2p1f 1
# define RVV_MATH_IS_LOG2F_FAMILY_log2p1f_clean_assert 1
# define RVV_MATH_IS_LOG2F_FAMILY_log2p1f_original 1
# define RVV_MATH_IS_LOG2F_FAMILY_INNER(name) \
  RVV_MATH_IS_LOG2F_FAMILY_##name
# define RVV_MATH_IS_LOG2F_FAMILY(name) \
  RVV_MATH_IS_LOG2F_FAMILY_INNER (name)
#endif

/* The log2f-original skeleton is also shared by the later log2p1f clean
   entry, whose masked scale operation needs a compile-time selector.  */
#ifndef RVV_MATH_IS_LOG2P1F_CLEAN_log2p1f
# define RVV_MATH_IS_LOG2P1F_CLEAN_log2p1f 1
# define RVV_MATH_IS_LOG2P1F_CLEAN_log2p1f_clean_assert 1
# define RVV_MATH_IS_LOG2P1F_CLEAN_INNER(name) \
  RVV_MATH_IS_LOG2P1F_CLEAN_##name
# define RVV_MATH_IS_LOG2P1F_CLEAN(name) \
  RVV_MATH_IS_LOG2P1F_CLEAN_INNER (name)
#endif

static const uint64_t class_sNaN = 0x100;
static const uint64_t class_NaN = 0x300;
static const uint64_t class_negInf = 0x1;
static const uint64_t class_posInf = 0x80;
static const uint64_t class_Inf = 0x81;
static const uint64_t class_negZero = 0x8;
static const uint64_t class_Zero = 0x18;
static const uint64_t class_posZero = 0x10;
static const uint64_t class_posDenorm = 0x20;
static const uint64_t class_Subnormal = 0x24;
static const uint64_t class_negative = 0x7;
static const uint64_t class_finite_neg = 0x6;
static const uint64_t class_finite_pos = 0x60;

/* Compile-time constants shared by the source-shaped single-precision
   expm1-family graphs.  These values are scalar constants, rather than
   indexed data, so keep them in the math header instead of a data TU.  */
#define RVV_SLEEF_EXPM1F_R_LN2 0x1.715476p+0f
#define RVV_SLEEF_EXPM1F_L2U 0x1.62e400p-1f
#define RVV_SLEEF_EXPM1F_L2L 0x1.7f7d1cp-20f
#define RVV_SLEEF_EXPM1F_MAX 88.72283172607421875f
#define RVV_SLEEF_EXPM1F_MIN -16.635532333438687f
#define RVV_SLEEF_EXPM1F_UNDERFLOW -104.0f
#define RVV_SLEEF_EXPM1F_POLY0 0.1980960224e-3f
#define RVV_SLEEF_EXPM1F_POLY1 0.1394256484e-2f
#define RVV_SLEEF_EXPM1F_POLY2 0.8333456703e-2f
#define RVV_SLEEF_EXPM1F_POLY3 0.4166637361e-1f
#define RVV_SLEEF_EXPM1F_TAYLOR_C \
  0.166666659414234244790680580464e+0f
#define RVV_SLEEF_EXPM1F_MAX_FINITE 0x1.fffffep+127f
#define RVV_SLEEF_EXPM1F_MIN_SUBNORMAL 0x1.0p-149f
#define RVV_SLEEF_EXPM1F_HALF 0.5f
#define RVV_SLEEF_EXPM1F_TWO 2.0f

/* The expk2f large-input coshf path uses the source split multiplied by two.
   Keep the already-rounded binary32 products as constants.  */
#define RVV_SLEEF_COSHF_TWO_L2U 0x1.62e400p+0f
#define RVV_SLEEF_COSHF_TWO_L2L 0x1.7f7d1cp-19f

#define RVV_LIBMVEC_COSHF_LARGE_INPUT 88.0f
#define RVV_LIBMVEC_COSHF_TINY_INPUT 0x1p-48f
#define RVV_LIBMVEC_COSHF_OVERFLOW 89.415985f
#define RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_HI RVV_SLEEF_COSHF_TWO_L2U
#define RVV_LIBMVEC_COSHF_LARGE_NEUTRAL_LO RVV_SLEEF_COSHF_TWO_L2L
#define RVV_SLEEF_COSHF_SOURCE_THRESHOLD 89.0f

/* This is the rounded low part used by the binary32 clean-side split.  */
#define RVV_LIBMVEC_LN2_LO_F -1.904654323148236017e-09f
#define RVV_LIBMVEC_EXP2M1F_MAX 128.0f
#define RVV_LIBMVEC_EXP2M1F_MIN -25.0f
#define RVV_LIBMVEC_EXP2M1F_EXACT_INT_LIMIT 24.0f

/* Scalar constants shared by the SLEEF single-precision trigonometric
   reduction paths.  Keep the actual rempif lookup table in
   v_f_sleef_trig_data.c; these values are compile-time constants, not
   indexed data.  */
#define RVV_SLEEF_TRIG_ONE_OVER_PI 0x1.45f306p-2f
#define RVV_SLEEF_TRIG_PI_F 0x1.921fb6p+1f
#define RVV_PI_F 0x1.921fb6p+1f
#define RVV_PI_LO_F -0x1.777a5cp-24f
#define RVV_ONE_OVER_PI_F RVV_SLEEF_TRIG_ONE_OVER_PI
#define RVV_ONE_OVER_PI_LO_F 0x1.b93910p-27f
#define RVV_PIBY2_F 0x1.921fb6p+0f
#define RVV_PIBY2_LO_F -0x1.777a5cp-25f
#define RVV_PIBY4_F 0x1.921fb6p-1f
#define RVV_PIBY4_LO_F -0x1.777a5cp-26f
#define RVV_SLEEF_TRIG_RANGE_MAX2 125.0f
#define RVV_SLEEF_TRIG_TWO_OVER_PI 0x1.45f306p-1f
#define RVV_SLEEF_TRIG_NEG_PIBY2_A2_F \
  (-3.1414794921875f * 0.5f)
#define RVV_SLEEF_TRIG_NEG_PIBY2_B2_F \
  (-0.00011315941810607910156f * 0.5f)
#define RVV_SLEEF_TRIG_NEG_PIBY2_C2_F \
  (-1.9841872589410058936e-09f * 0.5f)
#define RVV_SLEEF_TRIG_REMPI_TWO23 8388608.0f
#define RVV_SLEEF_TRIG_REMPI_FOUR 4.0f
#define RVV_SLEEF_TRIG_REMPI_NEG_QUARTER -0.25f
#define RVV_SLEEF_TRIG_REMPI_NEG_FOUR -4.0f
#define RVV_SLEEF_TRIG_REMPI_TWO_PI_HI \
  (3.1415927410125732422f * 2.0f)
#define RVV_SLEEF_TRIG_REMPI_TWO_PI_LO \
  (-8.7422776573475857731e-08f * 2.0f)
#define RVV_SLEEF_TRIG_REMPI_SMALL 0.7f

/* SLEEF sincosf's paired approximation coefficients.  */
#define RVV_SINCOSF_SIN_POLY0_F -0.000195169282960705459117889f
#define RVV_SINCOSF_SIN_POLY1_F 0.00833215750753879547119141f
#define RVV_SINCOSF_SIN_POLY2_F -0.166666537523269653320312f
#define RVV_SINCOSF_COS_POLY0_F -2.71811842367242206819355e-07f
#define RVV_SINCOSF_COS_POLY1_F 2.47990446951007470488548e-05f
#define RVV_SINCOSF_COS_POLY2_F -0.00138888787478208541870117f
#define RVV_SINCOSF_COS_POLY3_F 0.0416666641831398010253906f
#define RVV_SINCOSF_COS_POLY4_F -0.5f

#define RVV_SLEEF_SINF_NEG_PI_A2 -3.1414794921875f
#define RVV_SLEEF_SINF_NEG_PI_B2 -0.00011315941810607910156f
#define RVV_SLEEF_SINF_NEG_PI_C2 -1.9841872589410058936e-09f
#define RVV_SLEEF_SINF_REMPI_HALF_PI_HI \
  (3.1415927410125732422f * -0.5f)
#define RVV_SLEEF_SINF_REMPI_HALF_PI_LO \
  (-8.7422776573475857731e-08f * -0.5f)
#define RVV_SLEEF_SINF_REMPI_EXPONENT_BIAS 25
#define RVV_SLEEF_SINF_REMPI_SCALE_LIMIT (90 - 25)
#define RVV_SLEEF_SINF_REMPI_SCALE -64
#define RVV_SLEEF_SINF_REMPI_INDEX_SHIFT 4
#define RVV_SLEEF_TRIG_NEG_PI_A2_F RVV_SLEEF_SINF_NEG_PI_A2
#define RVV_SLEEF_TRIG_NEG_PI_B2_F RVV_SLEEF_SINF_NEG_PI_B2
#define RVV_SLEEF_TRIG_NEG_PI_C2_F RVV_SLEEF_SINF_NEG_PI_C2
#define RVV_SLEEF_SINF_POLY0_F RVV_SLEEF_COSF_POLY_0
#define RVV_SLEEF_SINF_POLY1_F RVV_SLEEF_COSF_POLY_1
#define RVV_SLEEF_SINF_POLY2_F RVV_SLEEF_COSF_POLY_2
#define RVV_SLEEF_SINF_POLY3_F RVV_SLEEF_COSF_POLY_3
#define RVV_SLEEF_TRIG_SIN_POLY0_F RVV_SLEEF_COSF_POLY_0
#define RVV_SLEEF_TRIG_SIN_POLY1_F RVV_SLEEF_COSF_POLY_1
#define RVV_SLEEF_TRIG_SIN_POLY2_F RVV_SLEEF_COSF_POLY_2
#define RVV_SLEEF_TRIG_SIN_POLY3_F RVV_SLEEF_COSF_POLY_3

#define RVV_SLEEF_COSF_NEG_PI_A2_HALF \
  (-3.1414794921875f * 0.5f)
#define RVV_SLEEF_COSF_NEG_PI_B2_HALF \
  (-0.00011315941810607910156f * 0.5f)
#define RVV_SLEEF_COSF_NEG_PI_C2_HALF \
  (-1.9841872589410058936e-09f * 0.5f)
#define RVV_SLEEF_COSF_REMPI_HALF_PI_HI \
  (-3.1415927410125732422f * 0.5f)
#define RVV_SLEEF_COSF_REMPI_HALF_PI_LO \
  (-8.7422776573475857731e-08f * -0.5f)
#define RVV_SLEEF_COSF_REDUCE_NEG_HALF -0.5f

/* Fixed coefficients for the fastcosf clean and source-original paths.
   Keep float and double spellings distinct for source scalar conversions.  */
#ifndef RVV_FASTCOSF_CONSTANTS
# define RVV_FASTCOSF_CONSTANTS
# define RVV_FASTCOSF_CLEAN_POLY0_F -0.1881748176e-3f
# define RVV_FASTCOSF_CLEAN_POLY1_F +0.8323502727e-2f
# define RVV_FASTCOSF_CLEAN_POLY2_F -0.1666651368e+0f
# define RVV_FASTCOSF_SOURCE_POLY0_D -0.1881748176e-3
# define RVV_FASTCOSF_SOURCE_POLY1_D +0.8323502727e-2
# define RVV_FASTCOSF_SOURCE_POLY2_D -0.1666651368e+0
# define RVV_FASTCOSF_CLEAN_M1_PI_F RVV_SLEEF_TRIG_ONE_OVER_PI
# define RVV_FASTCOSF_CLEAN_PI_F RVV_SLEEF_TRIG_PI_F
# define RVV_FASTCOSF_SOURCE_M1_PI_D 0.318309886183790671537767526745028724
# define RVV_FASTCOSF_SOURCE_PI_D 3.141592653589793238462643383279502884
# define RVV_FASTCOSF_FAST_RANGE_F 30.0f
# define RVV_FASTCOSF_XCOSF_NEG_PI_A_HALF (-3.140625f * 0.5f)
# define RVV_FASTCOSF_XCOSF_NEG_PI_B_HALF \
  (-0.0009670257568359375f * 0.5f)
# define RVV_FASTCOSF_XCOSF_NEG_PI_C_HALF \
  (-6.2771141529083251953e-07f * 0.5f)
# define RVV_FASTCOSF_XCOSF_NEG_PI_D_HALF \
  (-1.2154201256553420762e-10f * 0.5f)
# define RVV_FASTCOSF_XCOSF_TRIG_RANGE_MAX_F 39000.0f
#endif

/* Fixed coefficients shared by the erff clean and source-original paths.  */
#ifndef RVV_ERFF_CONSTANTS
# define RVV_ERFF_CONSTANTS
# define RVV_ERFF_POLY25_0_F -0.4360447008e-6f
# define RVV_ERFF_POLY25_1_F +0.6867515367e-5f
# define RVV_ERFF_POLY25_2_F -0.3045156700e-4f
# define RVV_ERFF_POLY25_3_F +0.9808536561e-4f
# define RVV_ERFF_POLY25_4_F +0.2395523916e-3f
# define RVV_ERFF_POLY25_5_F +0.1459901541e-3f
# define RVV_ERFF_POLY_GT_0_F -0.1130012848e-6f
# define RVV_ERFF_POLY_GT_1_F +0.4115272986e-5f
# define RVV_ERFF_POLY_GT_2_F -0.6928304356e-4f
# define RVV_ERFF_POLY_GT_3_F +0.7172692567e-3f
# define RVV_ERFF_POLY_GT_4_F -0.5131045356e-2f
# define RVV_ERFF_POLY_GT_5_F +0.2708637156e-1f
# define RVV_ERFF_DF25_0_0_F 0.0092883445322513580322f
# define RVV_ERFF_DF25_0_1_F -2.7863745897025330755e-11f
# define RVV_ERFF_DF25_1_0_F 0.042275499552488327026f
# define RVV_ERFF_DF25_1_1_F 1.3461399289988106057e-09f
# define RVV_ERFF_DF25_2_0_F 0.070523701608180999756f
# define RVV_ERFF_DF25_2_1_F -3.6616309318707365163e-09f
# define RVV_ERFF_DF_GT_0_0_F -0.11064319312572479248f
# define RVV_ERFF_DF_GT_0_1_F 3.7050452777225283007e-09f
# define RVV_ERFF_DF_GT_1_0_F -0.63192230463027954102f
# define RVV_ERFF_DF_GT_1_1_F -2.0200432585073177859e-08f
# define RVV_ERFF_DF_GT_2_0_F -1.1296638250350952148f
# define RVV_ERFF_DF_GT_2_1_F 2.5515120196453259252e-08f
# define RVV_ERFF_EXP_POLY_0_F 0.00136324646882712841033936f
# define RVV_ERFF_EXP_POLY_1_F 0.00836596917361021041870117f
# define RVV_ERFF_EXP_POLY_2_F 0.0416710823774337768554688f
# define RVV_ERFF_EXP_POLY_3_F 0.166665524244308471679688f
# define RVV_ERFF_EXP_POLY_4_F 0.499999850988388061523438f
# define RVV_ERFF_R_LN2_F RVV_SLEEF_EXPM1F_R_LN2
# define RVV_ERFF_L2U_F RVV_SLEEF_EXPM1F_L2U
# define RVV_ERFF_L2L_F RVV_SLEEF_EXPM1F_L2L
# define RVV_ERFF_TWO_OVER_SQRT_PI_F -1.1283792257308959961f
# define RVV_ERFF_TWO_OVER_SQRT_PI_LO_F 5.8635383422197591097e-08f
# define RVV_ERFF_RANGE_SPLIT_F 2.5f
# define RVV_ERFF_SATURATION_THRESHOLD_F 6.0f
# define RVV_ERFF_TINY_THRESHOLD_F 1e-4f
# define RVV_ERFF_EXP_UNDERFLOW_F RVV_SLEEF_EXPM1F_UNDERFLOW
#endif

/* Fixed coefficients shared by the erfcf clean and source-original paths.  */
#ifndef RVV_ERFCF_CONSTANTS
# define RVV_ERFCF_CONSTANTS
# define RVV_ERFCF_EXP_POLY_0_F RVV_SLEEF_EXPM1F_POLY0
# define RVV_ERFCF_EXP_POLY_1_F RVV_SLEEF_EXPM1F_POLY1
# define RVV_ERFCF_EXP_POLY_2_F RVV_SLEEF_EXPM1F_POLY2
# define RVV_ERFCF_EXP_POLY_3_F RVV_SLEEF_EXPM1F_POLY3
# define RVV_ERFCF_INTERVAL_0_F 1.0f
# define RVV_ERFCF_INTERVAL_1_F 2.2f
# define RVV_ERFCF_INTERVAL_2_F 4.3f
# define RVV_ERFCF_INTERVAL_3_F 10.1f
# define RVV_ERFCF_POLY_0_0_F -0.8638041618e-4f
# define RVV_ERFCF_POLY_0_1_F -0.6236977242e-5f
# define RVV_ERFCF_POLY_0_2_F -0.3869504035e+0f
# define RVV_ERFCF_POLY_0_3_F +0.1115344167e+1f
# define RVV_ERFCF_POLY_1_0_F +0.6000166177e-3f
# define RVV_ERFCF_POLY_1_1_F +0.5749821503e-4f
# define RVV_ERFCF_POLY_1_2_F +0.1288077235e+1f
# define RVV_ERFCF_POLY_1_3_F -0.9454904199e+0f
# define RVV_ERFCF_POLY_2_0_F -0.1665703603e-2f
# define RVV_ERFCF_POLY_2_1_F +0.6002851478e-5f
# define RVV_ERFCF_POLY_2_2_F -0.1816803217e+1f
# define RVV_ERFCF_POLY_2_3_F -0.3667259514e+0f
# define RVV_ERFCF_POLY_3_0_F +0.1795156277e-3f
# define RVV_ERFCF_POLY_3_1_F -0.2851036377e-2f
# define RVV_ERFCF_POLY_3_2_F +0.1249150872e+1f
# define RVV_ERFCF_POLY_3_3_F +0.7155663371e+0f
# define RVV_ERFCF_POLY_4_0_F +0.1914106123e-1f
# define RVV_ERFCF_POLY_4_1_F +0.2260518074e-1f
# define RVV_ERFCF_POLY_4_2_F -0.1328857988e+0f
# define RVV_ERFCF_POLY_4_3_F -0.1262947265e-1f
# define RVV_ERFCF_DF_0_0_D -0.102775359343930288081655368891e+0
# define RVV_ERFCF_DF_0_1_D -0.105247583459338632253369014063e+0
# define RVV_ERFCF_DF_0_2_D -0.482365310333045318680618892669e+0
# define RVV_ERFCF_DF_0_3_D -0.498961546254537647970305302739e+0
# define RVV_ERFCF_DF_1_0_D -0.636619483208481931303752546439e+0
# define RVV_ERFCF_DF_1_1_D -0.635609463574589034216723775292e+0
# define RVV_ERFCF_DF_1_2_D -0.134450203224533979217859332703e-2
# define RVV_ERFCF_DF_1_3_D -0.471199543422848492080722832666e-4
# define RVV_ERFCF_DF_2_0_D -0.112837917790537404939545770596e+1
# define RVV_ERFCF_DF_2_1_D -0.112855987376668622084547028949e+1
# define RVV_ERFCF_DF_2_2_D -0.572319781150472949561786101080e+0
# define RVV_ERFCF_DF_2_3_D -0.572364030327966044425932623525e+0
# define RVV_ERFCF_R_LN2_F RVV_ERFF_R_LN2_F
# define RVV_ERFCF_NEG_L2U_F (-RVV_ERFF_L2U_F)
# define RVV_ERFCF_NEG_L2L_F (-RVV_ERFF_L2L_F)
# define RVV_ERFCF_EXP_TAYLOR_C_F RVV_SLEEF_EXPM1F_TAYLOR_C
# define RVV_ERFCF_EXP_HALF_F RVV_SLEEF_EXPM1F_HALF
# define RVV_ERFCF_SOURCE_UNDERFLOW_F RVV_SLEEF_EXPM1F_UNDERFLOW
#endif
#define RVV_SLEEF_COSF_REDUCE_TWO 2.0f
#define RVV_SLEEF_COSF_POLY_0 2.6083159809786593541503e-06f
#define RVV_SLEEF_COSF_POLY_1 -0.0001981069071916863322258f
#define RVV_SLEEF_COSF_POLY_2 0.00833307858556509017944336f
#define RVV_SLEEF_COSF_POLY_3 -0.166666597127914428710938f

/* Scalar constants for the SLEEF single-precision cbrtf graph.  */
#define RVV_SLEEF_CBRTF_VILOGBK_THRESH 5.421010862427522e-20f
#define RVV_SLEEF_CBRTF_VILOGBK_SCALE 1.8446744073709552e19f
#define RVV_SLEEF_CBRTF_VILOGBK_MASK 0xffu
#define RVV_SLEEF_CBRTF_VILOGBK_SMALL_BIAS 0x7f
#define RVV_SLEEF_CBRTF_VILOGBK_LARGE_BIAS (64 + 0x7f)
#define RVV_SLEEF_CBRTF_Q2_1_HI 1.2599210739135742188f
#define RVV_SLEEF_CBRTF_Q2_1_LO -2.4018701694217270415e-08f
#define RVV_SLEEF_CBRTF_Q2_2_HI 1.5874010324478149414f
#define RVV_SLEEF_CBRTF_Q2_2_LO 1.9520385308169352356e-08f
#define RVV_SLEEF_CBRTF_ONE_THIRD (1.0f / 3.0f)
#define RVV_SLEEF_CBRTF_NEG_ONE_THIRD (-1.0f / 3.0f)
#define RVV_SLEEF_CBRTF_NEG_TWO_THIRDS (-2.0f / 3.0f)
#define RVV_SLEEF_CBRTF_QUOTIENT_OFFSET_F 6144.0f
#define RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR_F 3.0f
#define RVV_SLEEF_CBRTF_EXPONENT_OFFSET 1
#define RVV_SLEEF_CBRTF_REMAINDER_ONE 1
#define RVV_SLEEF_CBRTF_REMAINDER_TWO 2
#define RVV_SLEEF_CBRTF_QUOTIENT_OFFSET 6144
#define RVV_SLEEF_CBRTF_QUOTIENT_MUL 2731
#define RVV_SLEEF_CBRTF_QUOTIENT_SHIFT 13
#define RVV_SLEEF_CBRTF_QUOTIENT_DIVISOR 3
#define RVV_SLEEF_CBRTF_EXP_ADJUST 2048
#define RVV_SLEEF_CBRTF_POLY0_F -0.601564466953277587890625f
#define RVV_SLEEF_CBRTF_POLY1_F 2.8208892345428466796875f
#define RVV_SLEEF_CBRTF_POLY2_F -5.532182216644287109375f
#define RVV_SLEEF_CBRTF_POLY3_F 5.898262500762939453125f
#define RVV_SLEEF_CBRTF_POLY4_F -3.8095417022705078125f
#define RVV_SLEEF_CBRTF_POLY5_F 2.2241256237030029296875f

/* SLEEF's vcast_vf_f converts unsuffixed double expressions to float at the
   source operation point.  Keep that conversion tied to the surrounding
   source-shaped graph: GCC otherwise may hoist a data-backed conversion
   across the explicit-RNE vector-conversion window.  */
static inline float
RVV_SLEEF_ORIG_F64_TO_F32 (double value)
{
  float result;
  __asm__ __volatile__ ("fcvt.s.d %0, %1"
                        : "=f" (result)
                        : "f" (value)
                        : "memory");
  return result;
}

#ifdef _NEED_UNDEF_GNU_SOURCE_V_MATH
#undef _GNU_SOURCE
#undef _NEED_UNDEF_GNU_SOURCE_V_MATH
#endif

#endif

#ifndef LMUL
# error "Must assign an LMUL before including v_math.h"
#endif
#ifndef BIT_WIDTH
# error "Must assign BIT_WIDTH before including v_math.h"
#endif

/* Undefine previous LMUL-dependent macros to allow redefinition */
#undef MAKE_VTYPE
#undef MAKE_TYPE
#undef MAKE_FUNC
#undef MAKE_VBOOL
#undef VSET
#undef MAKE_REINTERPRET
#undef FLOAT
#undef VFLOAT
#undef VFLOAT_F
#undef INT
#undef VINT
#undef VINT_F
#undef UINT
#undef VUINT
#undef VUINT_F
#undef VBOOL
#undef VBOOL_F
#undef F_AS_I
#undef F_AS_I_F
#undef F_AS_U
#undef F_AS_U_F
#undef I_AS_F
#undef I_AS_F_F
#undef U_AS_F
#undef U_AS_F_F
#undef I_AS_U
#undef I_AS_U_F
#undef U_AS_I
#undef U_AS_I_F
#undef VMVI_VX
#undef VMVI_VX_F
#undef VMVU_VX
#undef VMVU_VX_F
#undef VFMV_VF
#undef VFMV_VF_F
#undef EXP_BIAS
#undef MAN_LEN
#undef EXP_BIAS_F
#undef MAN_LEN_F
#undef BIT_WIDTH_F
#undef fp_negOne_f
#undef fp_sNaN_f
#undef fp_qNaN_f

#define MAKE_VTYPE(A) __PASTE3 (A, BIT_WIDTH, __PASTE3 (m, LMUL, _t))
#define MAKE_TYPE(A) __PASTE3 (A, BIT_WIDTH, _t)
#define MAKE_FUNC(A) __PASTE3 (A, BIT_WIDTH, __PASTE2 (m, LMUL))

#if (BIT_WIDTH == 64)
# define EXP_BIAS 1023
# define MAN_LEN 52
# ifndef _V_MATH_SCALAR_64
#  define _V_MATH_SCALAR_64
typedef double scalar_t;
static const scalar_t fp_sNaN = __builtin_nans ("");
static const scalar_t fp_qNaN = __builtin_nan ("");
static const scalar_t fp_negInf = -__builtin_inf ();
static const scalar_t fp_negZero = -0.;
static const scalar_t fp_posZero = 0.0;
static const scalar_t fp_posOne = 0x1.0p0;
static const scalar_t fp_negOne = -0x1.0p0;
static const scalar_t fp_posHalf = 0x1.0p-1;
static const scalar_t fp_posInf = __builtin_inf ();
# endif
#elif (BIT_WIDTH == 32)
# define EXP_BIAS 127
# define MAN_LEN 23
# ifndef _V_MATH_SCALAR_32
#  define _V_MATH_SCALAR_32
typedef float scalar_t;
static const scalar_t fp_sNaN = __builtin_nansf ("");
static const scalar_t fp_qNaN = __builtin_nanf ("");
static const scalar_t fp_negInf = -__builtin_inff ();
static const scalar_t fp_negZero = -0.f;
static const scalar_t fp_posZero = 0.0f;
static const scalar_t fp_posOne = 0x1.0p0f;
static const scalar_t fp_negOne = -0x1.0p0f;
static const scalar_t fp_posHalf = 0x1.0p-1f;
static const scalar_t fp_posInf = __builtin_inff ();
# endif
#else
# error "requested BIT_WIDTH unsupported"
#endif

#define EXP_BIAS_F EXP_BIAS
#define MAN_LEN_F MAN_LEN
#define BIT_WIDTH_F BIT_WIDTH
#define fp_posOne_f fp_posOne
#define fp_posZero_f fp_posZero
#define fp_posInf_f fp_posInf
#define fp_negInf_f fp_negInf
#define fp_negZero_f fp_negZero
#define fp_negOne_f fp_negOne
#define fp_sNaN_f fp_sNaN
#define fp_qNaN_f fp_qNaN

/* For MAKE_VBOOL, the value is BIT_WIDTH/LMUL for the current element size. */
#if (BIT_WIDTH == 64)
# if (LMUL == 1)
#  define MAKE_VBOOL(A) __PASTE3 (A, 64, _t)
# elif (LMUL == 2)
#  define MAKE_VBOOL(A) __PASTE3 (A, 32, _t)
# elif (LMUL == 4)
#  define MAKE_VBOOL(A) __PASTE3 (A, 16, _t)
# elif (LMUL == 8)
#  define MAKE_VBOOL(A) __PASTE3 (A, 8, _t)
# else
#  error "requested LMUL unsupported"
# endif
#elif (BIT_WIDTH == 32)
# if (LMUL == 1)
#  define MAKE_VBOOL(A) __PASTE3 (A, 32, _t)
# elif (LMUL == 2)
#  define MAKE_VBOOL(A) __PASTE3 (A, 16, _t)
# elif (LMUL == 4)
#  define MAKE_VBOOL(A) __PASTE3 (A, 8, _t)
# elif (LMUL == 8)
#  define MAKE_VBOOL(A) __PASTE3 (A, 4, _t)
# else
#  error "requested LMUL unsupported"
# endif
#else
# error "requested BIT_WIDTH unsupported"
#endif

#define VSET __PASTE2 (__riscv_vsetvl_e, __PASTE3 (BIT_WIDTH, m, LMUL))
#define MAKE_REINTERPRET(A, B)                                                \
  __PASTE5 (__riscv_vreinterpret_v_, A, __PASTE4 (BIT_WIDTH, m, LMUL, _), B,  \
            __PASTE3 (BIT_WIDTH, m, LMUL))

#define FLOAT MAKE_TYPE (float)
#define VFLOAT MAKE_VTYPE (vfloat)
#define VFLOAT_F VFLOAT
#define INT MAKE_TYPE (int)
#define VINT MAKE_VTYPE (vint)
#define VINT_F VINT
#define UINT MAKE_TYPE (uint)
#define VUINT MAKE_VTYPE (vuint)
#define VUINT_F VUINT
#define VBOOL MAKE_VBOOL (vbool)
#define VBOOL_F VBOOL

#define F_AS_I MAKE_REINTERPRET (f, i)
#define F_AS_U MAKE_REINTERPRET (f, u)
#define I_AS_F MAKE_REINTERPRET (i, f)
#define U_AS_F MAKE_REINTERPRET (u, f)
#define I_AS_U MAKE_REINTERPRET (i, u)
#define U_AS_I MAKE_REINTERPRET (u, i)

#define VMVI_VX MAKE_FUNC (__riscv_vmv_v_x_i)
#define VMVU_VX MAKE_FUNC (__riscv_vmv_v_x_u)
#define VFMV_VF MAKE_FUNC (__riscv_vfmv_v_f_f)

/* Float-qualified spellings used by source-shaped SLEEF adapters.  The
   importer already selects the element width, so these are aliases rather
   than a second intrinsic family.  */
#define F_AS_I_F F_AS_I
#define F_AS_U_F F_AS_U
#define I_AS_F_F I_AS_F
#define U_AS_F_F U_AS_F
#define I_AS_U_F I_AS_U
#define U_AS_I_F U_AS_I
#define VMVI_VX_F VMVI_VX
#define VMVU_VX_F VMVU_VX
#define VFMV_VF_F VFMV_VF

#define IDENTIFY_F IDENTIFY

/* The scalar-source powf importer widens binary32 lanes to binary64 for its
   log2/exp2 core.  Keep the LMUL-specific widened types and intrinsics in the
   shared math header so the routine skeleton contains only the operation
   graph.  LMUL=8 is handled by the importer through two LMUL=4 calls.  */
#if (BIT_WIDTH == 32)
# if (LMUL == 1)
#  define RVV_POWF_VDOUBLE vfloat64m2_t
#  define RVV_POWF_VU64 vuint64m2_t
#  define RVV_POWF_VU32_IDX vuint32m1_t
#  define RVV_POWF_VFWCVT __riscv_vfwcvt_f_f_v_f64m2
#  define RVV_POWF_VFNCVT __riscv_vfncvt_f_f_w_f32m1
#  define RVV_POWF_VDFMV_VF __riscv_vfmv_v_f_f64m2
#  define RVV_POWF_VGATHER_F64 __riscv_vluxei32_v_f64m2
#  define RVV_POWF_VGATHER_U64 __riscv_vluxei64_v_u64m2
#  define RVV_POWF_F64_AS_U64 __riscv_vreinterpret_v_f64m2_u64m2
#  define RVV_POWF_U64_AS_F64 __riscv_vreinterpret_v_u64m2_f64m2
# elif (LMUL == 2)
#  define RVV_POWF_VDOUBLE vfloat64m4_t
#  define RVV_POWF_VU64 vuint64m4_t
#  define RVV_POWF_VU32_IDX vuint32m2_t
#  define RVV_POWF_VFWCVT __riscv_vfwcvt_f_f_v_f64m4
#  define RVV_POWF_VFNCVT __riscv_vfncvt_f_f_w_f32m2
#  define RVV_POWF_VDFMV_VF __riscv_vfmv_v_f_f64m4
#  define RVV_POWF_VGATHER_F64 __riscv_vluxei32_v_f64m4
#  define RVV_POWF_VGATHER_U64 __riscv_vluxei64_v_u64m4
#  define RVV_POWF_F64_AS_U64 __riscv_vreinterpret_v_f64m4_u64m4
#  define RVV_POWF_U64_AS_F64 __riscv_vreinterpret_v_u64m4_f64m4
# elif (LMUL == 4)
#  define RVV_POWF_VDOUBLE vfloat64m8_t
#  define RVV_POWF_VU64 vuint64m8_t
#  define RVV_POWF_VU32_IDX vuint32m4_t
#  define RVV_POWF_VFWCVT __riscv_vfwcvt_f_f_v_f64m8
#  define RVV_POWF_VFNCVT __riscv_vfncvt_f_f_w_f32m4
#  define RVV_POWF_VDFMV_VF __riscv_vfmv_v_f_f64m8
#  define RVV_POWF_VGATHER_F64 __riscv_vluxei32_v_f64m8
#  define RVV_POWF_VGATHER_U64 __riscv_vluxei64_v_u64m8
#  define RVV_POWF_F64_AS_U64 __riscv_vreinterpret_v_f64m8_u64m8
#  define RVV_POWF_U64_AS_F64 __riscv_vreinterpret_v_u64m8_f64m8
# endif
#endif

/* The source-shaped acosf helper is embedded by the acospif comparator.
   Select its public body from the importer function name so the helper
   skeleton has no routine-local enable/disable macro.  */
#ifndef RVV_MATH_IS_ACOSF_HELPERS_ONLY_acospif_original
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY_acospif_original 1
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_ACOSF_HELPERS_ONLY_##name
# define RVV_MATH_IS_ACOSF_HELPERS_ONLY(name) \
  RVV_MATH_IS_ACOSF_HELPERS_ONLY_INNER (name)
#endif

#ifndef RVV_LIBMVEC_TRIG_DOUBLE_HELPERS
#define RVV_LIBMVEC_TRIG_DOUBLE_HELPERS

#define PIBY2_INV 0x1.45f306dc9c883p-1
#define PIBY2_HI 0x1.921fb54442d18p+0
#define PIBY2_MID 0x1.1a62633145c07p-54
#define PIBY2_LO -0x1.f1976b7ed8fbcp-110
#define RVV_PI_HI_D 0x1.921fb54442d18p+1
#define RVV_PI_MID_D 0x1.1a62633145c07p-53
#define RVV_PI_Q60 0x3243f6a8885a308dLL
#define RVV_PIBY2_Q60 0x1921fb54442d1847LL
#define RVV_PIBY2_Q61 0x3243f6a8885a308dLL
#define RVV_ONE_OVER_PI_HI_D 0x1.45f306dc9c883p-2
#define RVV_ONE_OVER_PI_LO_D -0x1.6b01ec5417056p-56
#define RVV_ONE_OVER_PI_Q63 0x28be60db9391054aLL
#define RVV_ATAN2_PI_075_HI_D 0x1.2d97c7f3321d2p+1
#define RVV_ATAN2_PIBY4_HI_D 0x1.921fb54442d18p-1
#define RVV_ATAN2_PI_HI_D 0x1.921fb54442d18p+1
#define RVV_ATAN2_PI_LO_D 0x1.1a62633145c07p-53
#define RVV_ATAN2_PI_Q61 0x6487ed5110b4611aLL
#define RVV_ATAN2_HALF_Q61 0x1000000000000000LL
#define RVV_ATAN2_ONE_Q61 0x2000000000000000LL
#define RVV_ATAN2_ONE_OVER_4_D 0x1.0p-2
#define RVV_ATAN2_THREE_OVER_4_D 0x1.8p-1
#define RVV_LN2_INV_D 0x1.71547652b82fep+0
#define RVV_LN2_HI_D 0x1.62e42fefa39efp-1
#define RVV_LN2_LO_D 0x1.abc9e3b39803fp-56
#define RVV_LN2_INV_LO_D 0x1.777d0ffda0d24p-56
#define RVV_LOG10_2_HI_D 0x1.34413509f79ffp-2
#define RVV_LOG10_2_LO_D -0x1.9dc1da994fd00p-59
#define RVV_LOG10_E_HI_D 0x1.bcb7b1526e50ep-2
#define RVV_LOG10_E_LO_D 0x1.95355baaafad3p-57
#define RVV_EXP_X_MAX_D 0x1.65p+9
#define RVV_EXP_X_MIN_D -0x1.77p+9
#define RVV_LOG2_HI_D 0x1.62e42fefa4000p-1
#define RVV_LOG2_LO_D -0x1.8432a1b0e2634p-43
#define RVV_LOG2_BY2_HI_D 0x1.62e42fefa4000p-2
#define RVV_LOG2_BY2_LO_D -0x1.8432a1b0e2634p-44
#define RVV_ONE_Q60 0x1000000000000000LL
#define RVV_ONE_OVER_3_D 0x1.5555555555555p-2
#define RVV_THIRD_Q62 0x1555555555555555LL
#define RVV_ONE_Q62 0x4000000000000000LL
#define RVV_CBRT_2_Q62 0x50a28be635ca2b89LL
#define RVV_CBRT_4_Q62 0x6597fa94f5b8f20bLL
#define RVV_COSPI_R_PRIME_SCALE -0x1.0p-1
#define RVV_COSPI_POLY_RIGHT_SEED_A -0x1.8f5dd75850673p-37
#define RVV_COSPI_POLY_RIGHT_SEED_B 0x1.5d8b5ae12066ap-33
#define RVV_COSPI_POLY_RIGHT_OUTER_A -0x1.27e4f72551e3dp-22
#define RVV_COSPI_POLY_RIGHT_OUTER_B 0x1.71de35553ddb6p-19
#define RVV_COSPI_POLY_RIGHT_INNER_A 0x1.1ee950032f74cp-29
#define RVV_COSPI_POLY_RIGHT_INNER_B -0x1.ae5e4b94836f8p-26
#define RVV_COSPI_POLY_LEFT_SEED_A 0x1.a01a019b77545p-16
#define RVV_COSPI_POLY_LEFT_SEED_B -0x1.a01a019be932ap-13
#define RVV_COSPI_POLY_LEFT_OUTER_A 0x1.5555555555546p-5
#define RVV_COSPI_POLY_LEFT_OUTER_B -0x1.5555555555548p-3
#define RVV_COSPI_POLY_LEFT_INNER_A -0x1.6c16c16c1450cp-10
#define RVV_COSPI_POLY_LEFT_INNER_B 0x1.111111110f730p-7
#define RVV_COSPI_SMALL_EXP_OFFSET 5
#define RVV_COSPI_LARGE_EXP_OFFSET 53
#define RVV_HYPOT_SOURCE_MIN_D 0x1.0p-1022
#define RVV_HYPOT_SOURCE_SCALE_D 0x1.0p54
#define RVV_HYPOT_CLEAN_SCALE_UP_LIMIT_D 0x1.0p-968
#define RVV_HYPOT_CLEAN_SCALE_DOWN_LIMIT_D 0x1.0p+970
#define RVV_HYPOT_CLEAN_SCALE_UP_D 0x1.0p54
#define RVV_HYPOT_CLEAN_SCALE_DOWN_D 0x1.0p-54
#define RVV_HYPOTF_CLEAN_SCALE_LIMIT_F 0x1.0p-80f
#define RVV_HYPOTF_CLEAN_LARGE_LIMIT_F 0x1.0p80f
#define RVV_HYPOTF_CLEAN_SCALE_UP_F 0x1.0p48f
#define RVV_HYPOTF_CLEAN_SCALE_DOWN_F 0x1.0p-48f
#define RVV_HYPOTF_CLEAN_EXPONENT_GAP 16
#define RVV_HYPOTF_CLEAN_HALF_ULP_F 0x1.0p-24f
#define RVV_HYPOTF_SOURCE_FLT_MIN_F 0x1.0p-126f
#define RVV_HYPOTF_SOURCE_SCALE_F 0x1.0p24f
#define RVV_COSHF_X_MAX_D 0x1.65p9
#define RVV_COSHF_SMALL_EXP_OFFSET 5
#define RVV_COSHF_SCALE_LIMIT 60
#define RVV_COSHF_LARGE_N 50
#define RVV_COSH_P_EVEN0_D 0x1.555555555555ap-3
#define RVV_COSH_P_EVEN1_D 0x1.111111110ef6ap-7
#define RVV_COSH_P_EVEN2_D 0x1.a01a01b32b633p-13
#define RVV_COSH_P_EVEN3_D 0x1.71ddef82f4beep-19
#define RVV_COSH_P_EVEN4_D 0x1.af6eacd796f0bp-26
#define RVV_COSH_P_ODD0_D 0x1.5555555553aefp-5
#define RVV_COSH_P_ODD1_D 0x1.6c16c17a09506p-10
#define RVV_COSH_P_ODD2_D 0x1.a019b37a2b3dfp-16
#define RVV_COSH_P_ODD3_D 0x1.289788d8bdadfp-22
#define RVV_ERF_T_D 0x1.0p+1
#define RVV_ERF_ALPHA_HI_D 0x1.20dd750429b6dp+0
#define RVV_ERF_ALPHA_LO_D 0x1.1ae3a914fed80p-56
#define RVV_ERF_POLY_SCALE_D 0x1.0p+78
#define RVV_ERF_R_SCALE_D 0x1.0p+60
#define RVV_ERF_P_FIXED0_I 0x120dd750429b6d0fLL
#define RVV_ERF_P_FIXED1_I -0x1812746b0379e00cLL
#define RVV_ERF_P_FIXED2_I 0x1ce2f21a04292b5fLL
#define RVV_ERF_P_FIXED3_I -0x1b82ce31281b38e1LL
#define RVV_ERF_P_FIXED4_I 0x1565bcd0dd0bcd58LL
#define RVV_ERF_P_FIXED5_I -0xe016d9f815a019dLL
#define RVV_ERF_P_FIXED6_I 0x7e68c976c0ebcdcLL
#define RVV_ERF_P_FIXED7_I -0x3e9a49c76e6ee9aLL
#define RVV_ERF_P_FIXED8_I 0x1b9e64a589f8da9LL
#define RVV_ERF_P_TAIL0_D -0x1.5f70cd90f1878p-23
#define RVV_ERF_P_TAIL1_D 0x1.fca2b5f17c85ap-27
#define RVV_ERF_P_TAIL2_D -0x1.514eafaeffc30p-30
#define RVV_ERF_P_TAIL3_D 0x1.9b3583b6b826dp-34
#define RVV_ERF_P_TAIL4_D -0x1.c97ffcf4f4e22p-38
#define RVV_ERF_P_TAIL5_D 0x1.c2f4a46d3297dp-42
#define RVV_ERF_P_TAIL6_D -0x1.6ef3c7000b58bp-46
#define RVV_ERF_P_TAIL7_D 0x1.ac36453182837p-51
#define RVV_ERF_P_TAIL8_D -0x1.0482966738f0ep-56
#define RVV_ERF_Q_FIXED0_I 0xffff87b6641370fLL
#define RVV_ERF_Q_FIXED1_I -0x9062a79f9b29022LL
#define RVV_ERF_Q_FIXED2_I -0x11dc7e40e4efb77dLL
#define RVV_ERF_Q_FIXED3_I -0x1dd1004e1f59ed4LL
#define RVV_ERF_Q_FIXED4_I 0x1980c051527d41e7LL
#define RVV_ERF_Q_FIXED5_I 0x902cddcb829790bLL
#define RVV_ERF_Q_FIXED6_I -0x33d6f572cdbfa228LL
#define RVV_ERF_Q_FIXED7_I 0x425f9974bef87221LL
#define RVV_ERF_Q_FIXED8_I -0x5363e91dfca5d4dfLL
#define RVV_ERF_Q_TAIL0_D 0x1.b5eea4ad8cdbfp-16
#define RVV_ERF_Q_TAIL1_D -0x1.ded0a34468c8cp-18
#define RVV_ERF_Q_TAIL2_D 0x1.af4968b4d634ap-20
#define RVV_ERF_Q_TAIL3_D -0x1.51de51c57f11ap-22
#define RVV_ERF_Q_TAIL4_D 0x1.cbbf535e64b65p-25
#define RVV_ERF_Q_TAIL5_D -0x1.025a03d4fdf7bp-27
#define RVV_ERF_Q_TAIL6_D 0x1.c735f1e16e8cdp-31
#define RVV_ERF_Q_TAIL7_D -0x1.2de00f5eeee49p-34
#define RVV_ERF_Q_TAIL8_D 0x1.219bdcb68d070p-38
#define RVV_ERF_Q_TAIL9_D -0x1.7b5fc54357bcfp-43
#define RVV_ERF_Q_TAIL10_D 0x1.301ac8caec6e3p-48
#define RVV_ERF_Q_TAIL11_D -0x1.c3232aa28d427p-55
#define RVV_TWO_D 0x1.0p1
#define RVV_Q63_SCALE_D 0x1.0p63
#define RVV_Q60_SCALE_D 0x1.0p-60
#define RVV_ERFC_NEG_A_SCALED_D -0x1.ep+64
#define RVV_ERFC_B_FOR_TRANS_D 0x1.04p+2
#define RVV_ERFC_MIN_CLIP_D 0x1.0p-60
#define RVV_ERFC_MAX_CLIP_D 0x1.cp4
#define RVV_Q63_INV_SCALE_D 0x1.0p-63
#define RVV_ERFC_EXP_P_RIGHT_BASE_I 0x1acc4c50c4
#define RVV_ERFC_EXP_P_LEFT_BASE_I 0x888888889fe9c4
#define RVV_ERFC_P0_I 0x4f33682d757709e8LL
#define RVV_ERFC_P1_I -0x95970864bc25c71LL
#define RVV_ERFC_P2_I 0xa377a56796fd6fLL
#define RVV_ERFC_P3_I 0x5ea2d221c412d2dLL
#define RVV_ERFC_P4_I -0x8f0caa24847e2a3LL
#define RVV_ERFC_P5_I 0x8ac6781d49af506LL
#define RVV_ERFC_P6_I -0x67476ebb9bc1f58LL
#define RVV_ERFC_P7_I 0x3d0ed00f93b86cbLL
#define RVV_ERFC_P8_I -0x1c36fb9d9556ac0LL
#define RVV_ERFC_P9_I 0x96c3f45eaad23bLL
#define RVV_ERFC_P10_I -0x1a6d434ab9ada1LL
#define RVV_ERFC_P11_I -0x4dd9356c9c3f8LL
#define RVV_ERFC_P12_I 0x4bb31b11d0a1aLL
#define RVV_ERFC_P13_I -0xf2d325083d5bLL
#define RVV_ERFC_P14_I -0x52720383749fLL
#define RVV_ERFC_P15_I 0x33f7f3f6cb7dLL
#define RVV_ERFC_P16_I 0x4ed13a394fLL
#define RVV_ERFC_P17_I -0x770e9d9af50LL
#define RVV_ERFC_P18_I 0x108f3f3cf59LL
#define RVV_ERFC_P19_I 0x101b7f3c485LL
#define RVV_ERFC_P20_I -0x3ab6fb75adLL
#define RVV_ERFC_P21_I -0x237088721cLL
#define RVV_ERFC_P22_I 0x6ed93407eLL
#define RVV_ERFC_P23_I 0x3dbfb2c72LL
#define RVV_ERFC_EXP_P_RIGHT0_I 0x16c16c185646e2LL
#define RVV_ERFC_EXP_P_RIGHT1_I 0x3403401a3f740LL
#define RVV_ERFC_EXP_P_RIGHT2_I 0x680665cc2958LL
#define RVV_ERFC_EXP_P_RIGHT3_I 0xb8efdcde680LL
#define RVV_ERFC_EXP_P_RIGHT4_I 0x128afc94c08LL
#define RVV_ERFC_EXP_P_LEFT0_I 0x4000000000000000LL
#define RVV_ERFC_EXP_P_LEFT1_I 0x40000000000000edLL
#define RVV_ERFC_EXP_P_LEFT2_I 0x2000000000001659LL
#define RVV_ERFC_EXP_P_LEFT3_I 0xaaaaaaaaaaa201bLL
#define RVV_ERFC_EXP_P_LEFT4_I 0x2aaaaaaaaa03367LL
#define RVV_ERFINV_P_LEFT0_I -0x48dbe9f5b3eabaaLL
#define RVV_ERFINV_P_LEFT1_I -0xb35279f1a626ae5LL
#define RVV_ERFINV_P_LEFT2_I 0x33789911873d184aLL
#define RVV_ERFINV_P_LEFT3_I -0x1bf9138fc77c0fbfLL
#define RVV_ERFINV_P_LEFT4_I -0x2d4ec43bc48403d4LL
#define RVV_ERFINV_P_LEFT5_I 0x2d61deb53842cca1LL
#define RVV_ERFINV_P_LEFT6_I 0x23324eca6b3ff02LL
#define RVV_ERFINV_P_LEFT7_I -0xd4ec1d31542c4fcLL
#define RVV_ERFINV_P_LEFT8_I 0x2ecf3c60308b0f2LL
#define RVV_ERFINV_P_LEFT9_I 0x7c917b3378071eLL
#define RVV_ERFINV_P_LEFT10_I -0x1e09b597f226caLL
#define RVV_ERFINV_Q_LEFT0_I -0x52366e5b14c0970LL
#define RVV_ERFINV_Q_LEFT1_I -0xca57e95abcc599bLL
#define RVV_ERFINV_Q_LEFT2_I 0x3b6c91ec67f5759cLL
#define RVV_ERFINV_Q_LEFT3_I -0x1c40d5daa3be22bcLL
#define RVV_ERFINV_Q_LEFT4_I -0x41f11eb5d837386cLL
#define RVV_ERFINV_Q_LEFT5_I 0x3c6ce478fcd75c9aLL
#define RVV_ERFINV_Q_LEFT6_I 0xbb1cd7270cfba1dLL
#define RVV_ERFINV_Q_LEFT7_I -0x1988a4116498f1afLL
#define RVV_ERFINV_Q_LEFT8_I 0x44dc3042f103d20LL
#define RVV_ERFINV_Q_LEFT9_I 0x2390e683d02edf3LL
#define RVV_ERFINV_Q_LEFT10_I -0x8ec66f2a7e410cLL
#define RVV_ERFINV_P_RIGHT0_I 0x17a0bb69321dfLL
#define RVV_ERFINV_P_RIGHT1_I 0x402eb416ae6015LL
#define RVV_ERFINV_P_RIGHT2_I 0x2973eb18028ce34LL
#define RVV_ERFINV_P_RIGHT3_I 0x8034a7ece1d5370LL
#define RVV_ERFINV_P_RIGHT4_I 0xa76c08a74dae273LL
#define RVV_ERFINV_P_RIGHT5_I 0x11dd3876b83dd078LL
#define RVV_ERFINV_P_RIGHT6_I 0xfdd7693c3b77653LL
#define RVV_ERFINV_P_RIGHT7_I 0xb33d66152b3c223LL
#define RVV_ERFINV_P_RIGHT8_I 0x5a564c28c6a41a9LL
#define RVV_ERFINV_P_RIGHT9_I 0x1190449fe630213LL
#define RVV_ERFINV_P_RIGHT10_I -0x659c784274e1LL
#define RVV_ERFINV_Q_RIGHT0_I 0x17a09aabf9ceeLL
#define RVV_ERFINV_Q_RIGHT1_I 0x4030b9059ffcadLL
#define RVV_ERFINV_Q_RIGHT2_I 0x29b26b0d87f7855LL
#define RVV_ERFINV_Q_RIGHT3_I 0x87572a13d3fa2ddLL
#define RVV_ERFINV_Q_RIGHT4_I 0xd7a728b5620ac3cLL
#define RVV_ERFINV_Q_RIGHT5_I 0x1754392b473fd439LL
#define RVV_ERFINV_Q_RIGHT6_I 0x1791b9a091a816c2LL
#define RVV_ERFINV_Q_RIGHT7_I 0x167f71db9e13b075LL
#define RVV_ERFINV_Q_RIGHT8_I 0xcb9f5f3e5e618a4LL
#define RVV_ERFINV_Q_RIGHT9_I 0x68271fae767c68eLL
#define RVV_ERFINV_Q_RIGHT10_I 0x13745c4fa224b25LL
#define RVV_ERFINV_LOG_POLY_RIGHT0_D 0x1.74681ff881228p-14
#define RVV_ERFINV_LOG_POLY_RIGHT1_D 0x1.39751be23e4a3p-16
#define RVV_ERFINV_LOG_POLY_RIGHT2_D 0x1.30a893993e73dp-18
#define RVV_ERFINV_LOG_POLY_LEFT0_D 0x1.999999996ce82p-7
#define RVV_ERFINV_LOG_POLY_LEFT1_D 0x1.249249501b1adp-9
#define RVV_ERFINV_LOG_POLY_LEFT2_D 0x1.c71c47e7189f6p-12
#define RVV_ERFINV_LOG_POLY_LEFT3_D 0x1.55555555555dbp-4
#define RVV_ERFCINV_P_RIGHT0_I RVV_ERFINV_P_LEFT0_I
#define RVV_ERFCINV_P_RIGHT1_I RVV_ERFINV_P_LEFT1_I
#define RVV_ERFCINV_P_RIGHT2_I RVV_ERFINV_P_LEFT2_I
#define RVV_ERFCINV_P_RIGHT3_I RVV_ERFINV_P_LEFT3_I
#define RVV_ERFCINV_P_RIGHT4_I RVV_ERFINV_P_LEFT4_I
#define RVV_ERFCINV_P_RIGHT5_I RVV_ERFINV_P_LEFT5_I
#define RVV_ERFCINV_P_RIGHT6_I RVV_ERFINV_P_LEFT6_I
#define RVV_ERFCINV_P_RIGHT7_I RVV_ERFINV_P_LEFT7_I
#define RVV_ERFCINV_P_RIGHT8_I RVV_ERFINV_P_LEFT8_I
#define RVV_ERFCINV_P_RIGHT9_I RVV_ERFINV_P_LEFT9_I
#define RVV_ERFCINV_P_RIGHT10_I RVV_ERFINV_P_LEFT10_I
#define RVV_ERFCINV_Q_RIGHT0_I RVV_ERFINV_Q_LEFT0_I
#define RVV_ERFCINV_Q_RIGHT1_I RVV_ERFINV_Q_LEFT1_I
#define RVV_ERFCINV_Q_RIGHT2_I RVV_ERFINV_Q_LEFT2_I
#define RVV_ERFCINV_Q_RIGHT3_I RVV_ERFINV_Q_LEFT3_I
#define RVV_ERFCINV_Q_RIGHT4_I RVV_ERFINV_Q_LEFT4_I
#define RVV_ERFCINV_Q_RIGHT5_I RVV_ERFINV_Q_LEFT5_I
#define RVV_ERFCINV_Q_RIGHT6_I RVV_ERFINV_Q_LEFT6_I
#define RVV_ERFCINV_Q_RIGHT7_I RVV_ERFINV_Q_LEFT7_I
#define RVV_ERFCINV_Q_RIGHT8_I RVV_ERFINV_Q_LEFT8_I
#define RVV_ERFCINV_Q_RIGHT9_I RVV_ERFINV_Q_LEFT9_I
#define RVV_ERFCINV_Q_RIGHT10_I RVV_ERFINV_Q_LEFT10_I
#define RVV_ERFCINV_P_LEFT0_I RVV_ERFINV_P_RIGHT0_I
#define RVV_ERFCINV_P_LEFT1_I RVV_ERFINV_P_RIGHT1_I
#define RVV_ERFCINV_P_LEFT2_I RVV_ERFINV_P_RIGHT2_I
#define RVV_ERFCINV_P_LEFT3_I RVV_ERFINV_P_RIGHT3_I
#define RVV_ERFCINV_P_LEFT4_I RVV_ERFINV_P_RIGHT4_I
#define RVV_ERFCINV_P_LEFT5_I RVV_ERFINV_P_RIGHT5_I
#define RVV_ERFCINV_P_LEFT6_I RVV_ERFINV_P_RIGHT6_I
#define RVV_ERFCINV_P_LEFT7_I RVV_ERFINV_P_RIGHT7_I
#define RVV_ERFCINV_P_LEFT8_I RVV_ERFINV_P_RIGHT8_I
#define RVV_ERFCINV_P_LEFT9_I RVV_ERFINV_P_RIGHT9_I
#define RVV_ERFCINV_P_LEFT10_I RVV_ERFINV_P_RIGHT10_I
#define RVV_ERFCINV_Q_LEFT0_I RVV_ERFINV_Q_RIGHT0_I
#define RVV_ERFCINV_Q_LEFT1_I RVV_ERFINV_Q_RIGHT1_I
#define RVV_ERFCINV_Q_LEFT2_I RVV_ERFINV_Q_RIGHT2_I
#define RVV_ERFCINV_Q_LEFT3_I RVV_ERFINV_Q_RIGHT3_I
#define RVV_ERFCINV_Q_LEFT4_I RVV_ERFINV_Q_RIGHT4_I
#define RVV_ERFCINV_Q_LEFT5_I RVV_ERFINV_Q_RIGHT5_I
#define RVV_ERFCINV_Q_LEFT6_I RVV_ERFINV_Q_RIGHT6_I
#define RVV_ERFCINV_Q_LEFT7_I RVV_ERFINV_Q_RIGHT7_I
#define RVV_ERFCINV_Q_LEFT8_I RVV_ERFINV_Q_RIGHT8_I
#define RVV_ERFCINV_Q_LEFT9_I RVV_ERFINV_Q_RIGHT9_I
#define RVV_ERFCINV_Q_LEFT10_I RVV_ERFINV_Q_RIGHT10_I
#define RVV_ERFCINV_LOG_POLY_RIGHT0_D RVV_ERFINV_LOG_POLY_RIGHT0_D
#define RVV_ERFCINV_LOG_POLY_RIGHT1_D RVV_ERFINV_LOG_POLY_RIGHT1_D
#define RVV_ERFCINV_LOG_POLY_RIGHT2_D RVV_ERFINV_LOG_POLY_RIGHT2_D
#define RVV_ERFCINV_LOG_POLY_LEFT0_D RVV_ERFINV_LOG_POLY_LEFT0_D
#define RVV_ERFCINV_LOG_POLY_LEFT1_D RVV_ERFINV_LOG_POLY_LEFT1_D
#define RVV_ERFCINV_LOG_POLY_LEFT2_D RVV_ERFINV_LOG_POLY_LEFT2_D
#define RVV_ERFCINV_LOG_POLY_LEFT3_D RVV_ERFINV_LOG_POLY_LEFT3_D
#define RVV_ERFCINV_P_TINY_HI0_I -0x8593442eLL
#define RVV_ERFCINV_P_TINY_HI1_I -0x7f3dc156b1LL
#define RVV_ERFCINV_P_TINY_HI2_I -0x20dbbc67b1b8LL
#define RVV_ERFCINV_P_TINY_HI3_I -0x30c75b8264d44LL
#define RVV_ERFCINV_P_TINY_HI4_I -0x1ab06f39e5addeLL
#define RVV_ERFCINV_P_TINY_HI5_I -0x3d698fc2964a9cLL
#define RVV_ERFCINV_P_TINY_HI6_I 0x7f7d6a748d0c2fLL
#define RVV_ERFCINV_P_TINY_HI7_I 0x18b100818c77848LL
#define RVV_ERFCINV_P_TINY_HI8_I 0x4be1ed5d031774LL
#define RVV_ERFCINV_P_TINY_HI9_I 0x156a7e2e54e260LL
#define RVV_ERFCINV_P_TINY_LO0_I -0x4e7245b3LL
#define RVV_ERFCINV_P_TINY_LO1_I -0x1f0300096LL
#define RVV_ERFCINV_P_TINY_LO2_I -0xbc59b742LL
#define RVV_ERFCINV_P_TINY_LO3_I -0x18a421ab9LL
#define RVV_ERFCINV_P_TINY_LO4_I -0x180f2a477LL
#define RVV_ERFCINV_P_TINY_LO5_I 0xc3d4ab0bLL
#define RVV_ERFCINV_P_TINY_LO6_I 0x1729754e9LL
#define RVV_ERFCINV_P_TINY_LO7_I 0x1aca73439LL
#define RVV_ERFCINV_P_TINY_LO8_I -0x3b6c5afbLL
#define RVV_ERFCINV_P_TINY_LO9_I 0x1a0c336beLL
#define RVV_ERFCINV_Q_TINY_HI0_I -0x85933cdaLL
#define RVV_ERFCINV_Q_TINY_HI1_I -0x7f3de4b69fLL
#define RVV_ERFCINV_Q_TINY_HI2_I -0x20dd8dc1da27LL
#define RVV_ERFCINV_Q_TINY_HI3_I -0x30dc92d1cd231LL
#define RVV_ERFCINV_Q_TINY_HI4_I -0x1af5fcee397d58LL
#define RVV_ERFCINV_Q_TINY_HI5_I -0x42639eeec1d051LL
#define RVV_ERFCINV_Q_TINY_HI6_I 0x6182b99f6ca998LL
#define RVV_ERFCINV_Q_TINY_HI7_I 0x17a6848dc07624aLL
#define RVV_ERFCINV_Q_TINY_HI8_I 0x105ecd6aac52b12LL
#define RVV_ERFCINV_Q_TINY_HI9_I 0xbc944d4f8282afLL
#define RVV_ERFCINV_Q_TINY_LO0_I -0xb5b39d61LL
#define RVV_ERFCINV_Q_TINY_LO1_I -0x151d1cd35LL
#define RVV_ERFCINV_Q_TINY_LO2_I -0x1706945d7LL
#define RVV_ERFCINV_Q_TINY_LO3_I 0xabde03f9LL
#define RVV_ERFCINV_Q_TINY_LO4_I -0xc3530d28LL
#define RVV_ERFCINV_Q_TINY_LO5_I 0x662b41ecLL
#define RVV_ERFCINV_Q_TINY_LO6_I 0x938a5e35LL
#define RVV_ERFCINV_Q_TINY_LO7_I 0x8a0484b7LL
#define RVV_ERFCINV_Q_TINY_LO8_I 0x1d1e38258LL
#define RVV_ERFCINV_Q_TINY_LO9_I -0x155b50b48LL
#define RVV_ERFINV_DELTA_P0_LEFT_D -0x1.ec7dc41c17860p-2
#define RVV_ERFINV_DELTA_Q0_LEFT_D -0x1.29a0161e99446p-3
#define RVV_ERFINV_DELTA_P0_RIGHT_D -0x1.d622f4cbe0eeep-2
#define RVV_ERFINV_DELTA_Q0_RIGHT_D 0x1.f7e7557a34ae6p-2
/* erfcinv uses the complementary branch labels for these same split
   correction constants.  */
#define RVV_ERFCINV_DELTA_P0_RIGHT_D RVV_ERFINV_DELTA_P0_LEFT_D
#define RVV_ERFCINV_DELTA_Q0_RIGHT_D RVV_ERFINV_DELTA_Q0_LEFT_D
#define RVV_ERFCINV_DELTA_P0_LEFT_D RVV_ERFINV_DELTA_P0_RIGHT_D
#define RVV_ERFCINV_DELTA_Q0_LEFT_D RVV_ERFINV_DELTA_Q0_RIGHT_D
#define RVV_ERFCINV_N_ADJUST 60
#define RVV_ERFCINV_LEFT_THRESHOLD_D 0x1.2p-2
#define RVV_ERFCINV_TINY_THRESHOLD_D 0x1.0p-52
#define RVV_ERFCINV_LEFT_SEED_D 0x1.0p-3
#define RVV_ERFCINV_SCALE_X_D 0x1.0p60
#define RVV_ERFCINV_TINY_T_SCALE_D 0x1.0p64
#define RVV_ERFCINV_TINY_T_SCALE_INV_D 0x1.0p-64
#define RVV_ERFCINV_TINY_QUARANTINE_Y_D 0x1.0p4
#define RVV_ERFCINV_TINY_DF_SCALE_D 0x1.0p-34
#define RVV_ERFINV_RTPI_BY2_HI_D 0x1.c5bf891b4ef6ap-1
#define RVV_ERFINV_RTPI_BY2_LO_D 0x1.4f38760a41abbp-54
#define RVV_ERFINV_RIGHT_THRESHOLD_D 0x1.7p-1
#define RVV_ERFINV_RIGHT_SAFE_X_D 0x1.8p-1
#define RVV_ERFINV_SMALL_EXPONENT 30
#define RVV_ERFINV_N_ADJUST 0
#define RVV_LN10_HI_D 0x1.26bb1bbb55516p+1
#define RVV_LN10_LO_D -0x1.f48ad494ea3e9p-53
#define RVV_EXP10_P_INV_STD_D 0x1.a934f0979a371p+1
#define RVV_EXP10_P_HI_STD_D RVV_LOG10_2_HI_D
#define RVV_EXP10_P_LO_STD_D -0x1.9dc1da994fd21p-59
#define RVV_EXP10_X_MAX_D 0x1.36p+8
#define RVV_EXP10_X_MIN_D -0x1.46p+8
#define RVV_EXP_POLY_RIGHT0_D 0x1.71df804f1baa1p-19
#define RVV_EXP_POLY_RIGHT1_D 0x1.28aa3ea739296p-22
#define RVV_EXP_POLY_RIGHT2_D 0x1.acf86201fd199p-26
#define RVV_EXP_POLY_MID0_D 0x1.6c16c1825c970p-10
#define RVV_EXP_POLY_MID1_D 0x1.a01a00fe6f730p-13
#define RVV_EXP_POLY_MID2_D 0x1.a0199e1789c72p-16
#define RVV_EXP_POLY_LEFT0_D 0x1.55555555554d2p-3
#define RVV_EXP_POLY_LEFT1_D 0x1.5555555551307p-5
#define RVV_EXP_POLY_LEFT2_D 0x1.11111111309a4p-7
#define RVV_EXP_POLY_FINAL_D 0x1.0000000000007p-1
#define RVV_HALF_D 0x1.0p-1
#define RVV_EXPM1_X_MIN_D -0x1.5p+5
#define RVV_EXPM1_TINY_BOUND_D 0x1.0p-54
#define RVV_EXP10M1_TABLE_BOUND_D 0x1.a308a4198c9d7p-4
#define RVV_EXP10M1_SPECIAL_BOUND_D 306.0
#define RVV_EXP10M1_SCALE_BOUND_D 163840.0
#define RVV_EXP10M1_LOG10_2_D 0x1.a934f0979a371p+8
#define RVV_EXP10M1_LOG10_2_BY_N_HI_D 0x1.34413509f79ffp-9
#define RVV_EXP10M1_LOG10_2_BY_N_LO_D -0x1.9dc1da994fd21p-66
#define RVV_EXP10M1_ROUND2ZERO_D -0x1.34413509f79ffp-10
#define RVV_EXP10M1_TINY_BOUND_D 0x1.0p-1022
#define RVV_EXP10M1_USE_EXPM1_BOUND_D 0x1.0p-3
#define RVV_EXP10M1_ROUND_SHIFT_D 0x1.8p+52
#define RVV_EXP10M1_SPECIAL_OFFSET 0x6000000000000000ULL
#define RVV_EXP10M1_SPECIAL_BIAS1 0x7000000000000000ULL
#define RVV_EXP10M1_SPECIAL_BIAS2 0x3010000000000000ULL
#define RVV_EXP10M1_NEGZERO_BITS 0x8000000000000000ULL
#define RVV_EXP10M1_TABLE_INDEX_MASK 127
#define RVV_EXP10M1_SCALE_SHIFT 45
#define RVV_EXP10M1_SMALL_INDEX_MASK 0x3f
#define RVV_EXP10M1_SMALL_INDEX_OFFSET 24
#define RVV_EXP10M1_N_CLIP_MAX 64
#define RVV_EXP10M1_N_LE53 53
#define RVV_EXP10M1_N_LE_NEG54 -54
#define RVV_EXP2_OVERFLOW_BOUND_D 0x1.0p10
#define RVV_EXP2_UNDERFLOW_BOUND_D -0x1.f4p10
#define RVV_EXP2_POLY_HI0_D 0.7073164598085707425e-8
#define RVV_EXP2_POLY_HI1_D 0.4434359082926529454e-9
#define RVV_EXP2_POLY_MID_HI_A0_D 0.1321543872511327615e-5
#define RVV_EXP2_POLY_MID_HI_A1_D 0.1017819260921760451e-6
#define RVV_EXP2_POLY_MID_HI_B0_D 0.1540353045101147808e-3
#define RVV_EXP2_POLY_MID_HI_B1_D 0.1525273353517584730e-4
#define RVV_EXP2_POLY_MID_LO_A0_D 0.9618129107597600536e-2
#define RVV_EXP2_POLY_MID_LO_A1_D 0.1333355814670499073e-2
#define RVV_EXP2_POLY_MID_LO_B0_D 0.2402265069591012214e+0
#define RVV_EXP2_POLY_MID_LO_B1_D 0.5550410866482046596e-1
#define RVV_EXPM1_POLY_RIGHT0_D 0x1.71ddf7aef0679p-19
#define RVV_EXPM1_POLY_RIGHT1_D 0x1.27e4e210af311p-22
#define RVV_EXPM1_POLY_RIGHT2_D 0x1.af5ff637cd647p-26
#define RVV_EXPM1_POLY_RIGHT3_D 0x1.1f6562eae5ba9p-29
#define RVV_EXPM1_POLY_MID0_D 0x1.6c16c16c166f3p-10
#define RVV_EXPM1_POLY_MID1_D 0x1.a01a01b0207e3p-13
#define RVV_EXPM1_POLY_MID2_D 0x1.a01a01a4af90ap-16
#define RVV_EXPM1_POLY_LEFT0_D 0x1.5555555555559p-3
#define RVV_EXPM1_POLY_LEFT1_D 0x1.5555555555556p-5
#define RVV_EXPM1_POLY_LEFT2_D 0x1.111111110f62ap-7
#define RVV_EXP10M1_OUTER_POLY0_D 0x1.6c16c16c16c17p-10
#define RVV_EXP10M1_OUTER_POLY1_D 0x1.1111111111111p-7
#define RVV_EXP10M1_OUTER_POLY2_D 0x1.5555555555555p-5
#define RVV_EXP10M1_OUTER_POLY3_D 0x1.5555555555555p-3
#define RVV_EXP10M1_OUTER_POLY4_D 0x1.0000000000000p-1

/* Scalar constants for the double-precision exp2m1 kernel.  Its exp2 and
   scalem1 lookup tables remain indexed data in v_d_exp2m1_data.c.  */
#ifndef RVV_EXP2M1_CONSTANTS
# define RVV_EXP2M1_CONSTANTS
# define RVV_EXP2M1_SPECIAL_BOUND_D 1022.0
# define RVV_EXP2M1_SCALE_BOUND_D 1280.0
# define RVV_EXP2M1_TABLE_BOUND_D 0x1.5bfffffffffffp-2
# define RVV_EXP2M1_SPECIAL_OFFSET RVV_EXP10M1_SPECIAL_OFFSET
# define RVV_EXP2M1_SPECIAL_BIAS RVV_EXP10M1_SPECIAL_BIAS1
# define RVV_EXP2M1_SPECIAL_BIAS2 RVV_EXP10M1_SPECIAL_BIAS2
# define RVV_EXP2M1_SHIFT_D (0x1.8p52 / V_EXP2M1_N)
# define RVV_EXP2M1_LOG2_HI_D RVV_LN2_HI_D
# define RVV_EXP2M1_LOG2_LO_D RVV_LN2_LO_D
# define RVV_EXP2M1_C1_D 0x1.ebfbdff82c58ep-3
# define RVV_EXP2M1_C2_D 0x1.c6b08d71f5804p-5
# define RVV_EXP2M1_C3_D 0x1.3b2ab6fee7509p-7
# define RVV_EXP2M1_C4_D 0x1.5d1d37eb33b15p-10
# define RVV_EXP2M1_C5_D 0x1.423f35f371d9ap-13
# define RVV_EXP2M1_C6_D 0x1.e7d57ad9a5f93p-5
# define RVV_EXP2M1_INV_LN2_D RVV_LN2_INV_D
#endif
#define RVV_EXP2F_CLEAN_MAX_INPUT_F 128.0f
#define RVV_EXP2F_CLEAN_MIN_INPUT_F -150.0f
#define RVV_EXP2F_CLEAN_MAX_FINITE_F 0x1.fffffep127f
#define RVV_EXP2F_CLEAN_MIN_NORMAL_F 0x1.0p-126f
#define RVV_EXP2F_SOURCE_MAX_INPUT_F 128.0f
#define RVV_EXP2F_SOURCE_MIN_INPUT_F -150.0f
#define RVV_EXP2F_CLEAN_POLY0_F 0.1535920892e-3f
#define RVV_EXP2F_CLEAN_POLY1_F 0.1339262701e-2f
#define RVV_EXP2F_CLEAN_POLY2_F 0.9618384764e-2f
#define RVV_EXP2F_CLEAN_POLY3_F 0.5550347269e-1f
#define RVV_EXP2F_CLEAN_POLY4_F 0.2402264476e+0f
#define RVV_EXP2F_CLEAN_POLY5_F 0.6931471825e+0f
#define RVV_EXP2F_SOURCE_POLY0_D +0.1535920892e-3
#define RVV_EXP2F_SOURCE_POLY1_D +0.1339262701e-2
#define RVV_EXP2F_SOURCE_POLY2_D +0.9618384764e-2
#define RVV_EXP2F_SOURCE_POLY3_D +0.5550347269e-1
#define RVV_EXP2F_SOURCE_POLY4_D +0.2402264476e+0
#define RVV_EXP2F_SOURCE_POLY5_D +0.6931471825e+0
#define RVV_EXPF_CLEAN_R_LN2_F 0x1.715476p+0f
#define RVV_EXPF_CLEAN_L2U_F 0x1.62e400p-1f
#define RVV_EXPF_CLEAN_L2L_F 0x1.7f7d1cp-20f
#define RVV_EXPF_CLEAN_MAX_INPUT_F 100.0f
#define RVV_EXPF_CLEAN_MIN_INPUT_F -104.0f
#define RVV_EXPF_CLEAN_MAX_FINITE_F 0x1.fffffep127f
#define RVV_EXPF_CLEAN_MIN_NORMAL_F 0x1.0p-126f
#define RVV_EXPF_SOURCE_R_LN2_F \
  1.442695040888963407359924681001892137426645954152985934135449406931f
#define RVV_EXPF_SOURCE_L2U_F 0.693145751953125f
#define RVV_EXPF_SOURCE_L2L_F 1.428606765330187045e-06f
#define RVV_EXPF_SOURCE_MAX_INPUT_F 100.0f
#define RVV_EXPF_SOURCE_MIN_INPUT_F -104.0f
#define RVV_EXPF_CLEAN_POLY0_F 0.000198527617612853646278381f
#define RVV_EXPF_CLEAN_POLY1_F 0.00139304355252534151077271f
#define RVV_EXPF_CLEAN_POLY2_F 0.00833336077630519866943359f
#define RVV_EXPF_CLEAN_POLY3_F 0.0416664853692054748535156f
#define RVV_EXPF_CLEAN_POLY4_F 0.166666671633720397949219f
#define RVV_EXPF_CLEAN_POLY5_F 0.5f
#define RVV_EXPF_SOURCE_POLY0_D +0.000198527617612853646278381
#define RVV_EXPF_SOURCE_POLY1_D +0.00139304355252534151077271
#define RVV_EXPF_SOURCE_POLY2_D +0.00833336077630519866943359
#define RVV_EXPF_SOURCE_POLY3_D +0.0416664853692054748535156
#define RVV_EXPF_SOURCE_POLY4_D +0.166666671633720397949219
#define RVV_EXPF_SOURCE_POLY5_D +0.5
#define RVV_LOG1PF_BOUND_F 0x1.2ced32p+126f
#define RVV_LOG1PF_INV_075_F 0x1.555556p0f
#define RVV_LOG1PF_SCALE_F 0x1.0p64f
#define RVV_LOG1PF_FLT_MIN_F 0x1.0p-126f
#define RVV_SLEEF_LOG_U1_LN2_HI_F 0.69314718246459960938f
#define RVV_SLEEF_LOG_U1_LN2_LO_F -1.904654323148236017e-09f
#define RVV_SLEEF_LOG_U1_POLY0_F 0.3027294874e+0f
#define RVV_SLEEF_LOG_U1_POLY1_F 0.3996108174e+0f
#define RVV_SLEEF_LOG_U1_POLY2_F 0.6666694880e+0f
#define RVV_LOG1PF_LN2_HI_F RVV_SLEEF_LOG_U1_LN2_HI_F
#define RVV_LOG1PF_LN2_LO_F RVV_SLEEF_LOG_U1_LN2_LO_F
#define RVV_LOGP1F_BOUND_F RVV_LOG1PF_BOUND_F
#define RVV_LOGP1F_INV_075_F RVV_LOG1PF_INV_075_F
#define RVV_LOGP1F_SCALE_UP_F RVV_LOG1PF_SCALE_F
#define RVV_LOGP1F_SCALE_DOWN_F 0x1.0p-4f
#define RVV_LOGP1F_LN2_HI_F RVV_LOG1PF_LN2_HI_F
#define RVV_LOGP1F_LN2_LO_F RVV_LOG1PF_LN2_LO_F
#define RVV_LOG2F_CLEAN_BOUND_F RVV_LOG1PF_BOUND_F
#define RVV_LOG2F_CLEAN_SCALE_DOWN_F RVV_LOGP1F_SCALE_DOWN_F
#define RVV_LOG2F_CLEAN_LARGE_EXP_ADJUST 4
#define RVV_LOG2F_CLEAN_LARGE_EXP_ADJUST_F 4.0f
#define RVV_LOG2F_FLT_MIN_F RVV_LOG1PF_FLT_MIN_F
#define RVV_LOG2F_SCALE_F RVV_LOG1PF_SCALE_F
#define RVV_LOG2F_INV_075_F 1.3333333333333333333333333333333333f
#define RVV_LOG2F_SUBNORMAL_EXP_ADJUST 64
#define RVV_LOG2F_X_HI_F 2.8853900432586669922f
#define RVV_LOG2F_X_LO_F 3.2734474483568488616e-08f
#define RVV_LOG2F_POLY0_F 0.4374550283e+0f
#define RVV_LOG2F_POLY1_F 0.5764790177e+0f
#define RVV_LOG2F_POLY2_F 0.9618012905120f
#define RVV_LOG10F_CLEAN_LOG10_2_HI_F 0.30103001f
#define RVV_LOG10F_CLEAN_LOG10_2_LO_F -1.4320989e-08f
#define RVV_LOG10F_CLEAN_TWO_LOG10_E_HI_F 0.868588984f
#define RVV_LOG10F_CLEAN_TWO_LOG10_E_LO_F -2.1707573e-08f
#define RVV_LOG10F_CLEAN_SUBNORMAL_SCALE_F 0x1.0p64f
#define RVV_LOG10F_CLEAN_INV_075_F 0x1.555556p0f
#define RVV_LOG10F_SOURCE_LOG10_2_HI 0.30103001
#define RVV_LOG10F_SOURCE_LOG10_2_LO -1.432098889e-08
#define RVV_LOG10F_SOURCE_TWO_LOG10_E_HI 0.868588984
#define RVV_LOG10F_SOURCE_TWO_LOG10_E_LO -2.170757285e-08
#define RVV_LOG10F_SOURCE_SUBNORMAL_SCALE_F 0x1.0p64f
#define RVV_LOG10F_SOURCE_FLT_MIN_F 0x1.0p-126f
#define RVV_LOG10F_SOURCE_INV_075 (1.0 / 0.75)
#define RVV_LOG10F_SOURCE_SUBNORMAL_EXP_ADJUST 64
#define RVV_LOG10F_CLEAN_POLY0_F 0.1314289868e+0f
#define RVV_LOG10F_CLEAN_POLY1_F 0.1735493541e+0f
#define RVV_LOG10F_CLEAN_POLY2_F 0.2895309627e+0f
#define RVV_LOG10F_SOURCE_POLY0_D 0.1314289868
#define RVV_LOG10F_SOURCE_POLY1_D 0.1735493541
#define RVV_LOG10F_SOURCE_POLY2_D 0.2895309627
#define RVV_EXP10F_CLEAN_LOG10_2_F 0x1.a934f0p+1f
#define RVV_EXP10F_CLEAN_L10U_F 0x1.344000p-2f
#define RVV_EXP10F_CLEAN_L10L_F 0x1.3509f8p-18f
#define RVV_EXP10F_CLEAN_LN10_HI_F 2.3025851249694824219f
#define RVV_EXP10F_CLEAN_LN10_LO_F -3.1705172516493593157e-08f
#define RVV_EXP10F_CLEAN_MAX_INPUT_F 38.5318394191036238941387f
#define RVV_EXP10F_CLEAN_MIN_INPUT_F -50.0f
#define RVV_EXP10F_SOURCE_LOG10_2 3.3219280948873623478703194294893901759
#define RVV_EXP10F_SOURCE_L10U_F 0.301025390625f
#define RVV_EXP10F_SOURCE_L10L_F 4.605038981e-06f
#define RVV_EXP10F_SOURCE_LN10_HI 2.3025851249694824219
#define RVV_EXP10F_SOURCE_LN10_LO -3.1705172516493593157e-08
#define RVV_EXP10F_SOURCE_MAX_INPUT_F 38.5318394191036238941387f
#define RVV_EXP10F_SOURCE_MIN_INPUT_F -50.0f
#define RVV_EXP10F_CLEAN_POLY0_F 0.6802555919e-1f
#define RVV_EXP10F_CLEAN_POLY1_F 0.2078080326e+0f
#define RVV_EXP10F_CLEAN_POLY2_F 0.5393903852e+0f
#define RVV_EXP10F_CLEAN_POLY3_F 0.1171245337e+1f
#define RVV_EXP10F_CLEAN_POLY4_F 0.2034678698e+1f
#define RVV_EXP10F_CLEAN_POLY5_F 0.2650949001e+1f
#define RVV_EXP10F_SOURCE_POLY0_D +0.6802555919e-1
#define RVV_EXP10F_SOURCE_POLY1_D +0.2078080326e+0
#define RVV_EXP10F_SOURCE_POLY2_D +0.5393903852e+0
#define RVV_EXP10F_SOURCE_POLY3_D +0.1171245337e+1
#define RVV_EXP10F_SOURCE_POLY4_D +0.2034678698e+1
#define RVV_EXP10F_SOURCE_POLY5_D +0.2650949001e+1
#define RVV_EXP10M1F_CLEAN_LN10_HI_F RVV_EXP10F_CLEAN_LN10_HI_F
#define RVV_EXP10M1F_CLEAN_LN10_LO_F RVV_EXP10F_CLEAN_LN10_LO_F
#define RVV_EXP10M1F_SOURCE_R_LN2_F RVV_EXPF_CLEAN_R_LN2_F
#define RVV_EXP10M1F_SOURCE_L2U_F RVV_EXPF_CLEAN_L2U_F
#define RVV_EXP10M1F_SOURCE_L2L_F RVV_EXPF_CLEAN_L2L_F
#define RVV_EXP10M1F_SOURCE_EXPM1_MAX_F 88.72283172607421875f
#define RVV_EXP10M1F_SOURCE_UNDERFLOW_F -104.0f
#define RVV_EXP10M1F_SOURCE_EXPM1_TINY_F 0x1.0p-24f
#define RVV_EXP10M1F_CLEAN_MIN_INPUT_F -0x1.e1a5e2p+2f
#define RVV_EXP10M1F_CLEAN_MAX_INPUT_F 0x1.344134p+5f
#define RVV_EXP10M1F_CLEAN_MIN_NORMAL_F RVV_EXP2F_CLEAN_MIN_NORMAL_F
#define RVV_EXP10M1F_EXPM1_TAYLOR_C_F 0.166666659414234244790680580464e+0f
#define RVV_EXP10M1F_EXPM1_POLY0_F RVV_SLEEF_EXPM1F_POLY0
#define RVV_EXP10M1F_EXPM1_POLY1_F RVV_SLEEF_EXPM1F_POLY1
#define RVV_EXP10M1F_EXPM1_POLY2_F RVV_SLEEF_EXPM1F_POLY2
#define RVV_EXP10M1F_EXPM1_POLY3_F RVV_SLEEF_EXPM1F_POLY3
#define RVV_COSPIF_TINY_EXPONENT 99u
#define RVV_COSPIF_LARGE_EXPONENT 151u
#define RVV_COSPIF_TRIG_RANGE_MAX4_F 8.0e6f
#define RVV_COSPIF_POLY_COS0_F 0.3093842054e-6f
#define RVV_COSPIF_POLY_SIN0_F -0.2430611801e-7f
#define RVV_COSPIF_POLY_COS1_F -0.3657307388e-4f
#define RVV_COSPIF_POLY_SIN1_F 0.3590577080e-5f
#define RVV_COSPIF_POLY_COS2_F 0.2490393585e-2f
#define RVV_COSPIF_POLY_SIN2_F -0.3259917721e-3f
#define RVV_COSPIF_DF_FIRST_COS_HI_F -0.080745510756969451904f
#define RVV_COSPIF_DF_FIRST_COS_LO_F -1.3373665339076936258e-09f
#define RVV_COSPIF_DF_FIRST_SIN_HI_F 0.015854343771934509277f
#define RVV_COSPIF_DF_FIRST_SIN_LO_F 4.4940051354032242811e-10f
#define RVV_COSPIF_DF_SECOND_COS_HI_F 0.78539818525314331055f
#define RVV_COSPIF_DF_SECOND_COS_LO_F -2.1857338617566484855e-08f
#define RVV_COSPIF_DF_SECOND_SIN_HI_F -0.30842512845993041992f
#define RVV_COSPIF_DF_SECOND_SIN_LO_F -9.0728339030733922277e-09f
#define RVV_LGAMMA_OVERFLOW_X_D 0x1.754d9278b51a7p+1014
#define RVV_LGAMMA_NEG_BND_TO_2_SC_D 0x1p+3
#define RVV_LGAMMA_NEG_3_TO_BND_SC_D 0x1p+4
#define RVV_POW_LOG_LN2_HI_D 0x1.62e42fefa3800p-1
#define RVV_POW_LOG_LN2_LO_D 0x1.ef35793c76730p-45
#define RVV_POW_LOG_OFF 0x3fe6955500000000ULL
#define RVV_POW_EXP_SIGN_BIAS (0x800ULL << V_POW_EXP_TABLE_BITS)
#define RVV_POW_EXP_SMALL_TOP 0x3c9U
#define RVV_POW_EXP_THRES_TOP (0x408U - 0x3c9U)
#define RVV_POW_EXP_HUGE_TOP 0x409U
#define RVV_POW_EXP_INV_LN2_N_D (RVV_LN2_INV_D * (double) V_POW_EXP_N)
#define RVV_POW_EXP_NEG_LN2_HI_N_D -0x1.62e42fefa0000p-8
#define RVV_POW_EXP_NEG_LN2_LO_N_D -0x1.cf79abc9e3b3ap-47
#define RVV_POW_EXP_C2_D 0x1.ffffffffffdbdp-2
#define RVV_POW_EXP_C3_D 0x1.555555555543cp-3
#define RVV_POW_EXP_C4_D 0x1.55555cf172b91p-5
#define RVV_POW_EXP_C5_D 0x1.1111167a4d017p-7
#define RVV_POW_LOG_POLY0_D -0x1p-1
#define RVV_POW_LOG_POLY1_D -0x1.555555555556p-1
#define RVV_POW_LOG_POLY2_D 0x1.0000000000006p-1
#define RVV_POW_LOG_POLY3_D 0x1.999999959554ep-1
#define RVV_POW_LOG_POLY4_D -0x1.555555529a47ap-1
#define RVV_POW_LOG_POLY5_D -0x1.2495b9b4845e9p0
#define RVV_POW_LOG_POLY6_D 0x1.0002b8b263fc3p0

/* Fixed coefficients for the single-precision gamma helpers.  Preserve
   the source float/double types and the explicit conversion points.  */
/* Shared SLEEF single-precision logk2f polynomial.  */
#define RVV_SLEEF_LOGK2F_POLY0_F 0.2392828464508056640625f
#define RVV_SLEEF_LOGK2F_POLY1_F 0.28518211841583251953125f
#define RVV_SLEEF_LOGK2F_POLY2_F 0.400005877017974853515625f
#define RVV_SLEEF_LOGK2F_POLY3_F 0.666666686534881591796875f

#ifndef RVV_TGAMMAF_CONSTANTS
# define RVV_TGAMMAF_CONSTANTS
# define RVV_TGAMMAF_SOURCE_LOG_POLY_0_F RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_TGAMMAF_SOURCE_LOG_POLY_1_F RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_TGAMMAF_SOURCE_LOG_POLY_2_F RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_TGAMMAF_SOURCE_LOG_POLY_3_F RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_0_0_F RVV_COSPIF_POLY_SIN0_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_0_1_F RVV_COSPIF_POLY_COS0_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_1_0_F RVV_COSPIF_POLY_SIN1_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_1_1_F RVV_COSPIF_POLY_COS1_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_2_0_F RVV_COSPIF_POLY_SIN2_F
# define RVV_TGAMMAF_SOURCE_SINPIF_POLY_2_1_F RVV_COSPIF_POLY_COS2_F
# define RVV_TGAMMAF_SOURCE_SINPIF_C_0_0_D +0.015854343771934509277
# define RVV_TGAMMAF_SOURCE_SINPIF_C_0_1_D -0.080745510756969451904
# define RVV_TGAMMAF_SOURCE_SINPIF_C_1_0_D +4.4940051354032242811e-10
# define RVV_TGAMMAF_SOURCE_SINPIF_C_1_1_D -1.3373665339076936258e-09
# define RVV_TGAMMAF_SOURCE_SINPIF_P_0_0_D -0.30842512845993041992
# define RVV_TGAMMAF_SOURCE_SINPIF_P_0_1_D +0.78539818525314331055
# define RVV_TGAMMAF_SOURCE_SINPIF_P_1_0_D -9.0728339030733922277e-09
# define RVV_TGAMMAF_SOURCE_SINPIF_P_1_1_D -2.1857338617566484855e-08
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_0_D +0.000839498720672087279971000786
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_1_D -5.17179090826059219329394422e-05
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_2_D -0.000592166437353693882857342347
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_3_D +6.97281375836585777403743539e-05
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_4_D +0.000784039221720066627493314301
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_5_D -0.000229472093621399176949318732
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_6_D -0.002681327160493827160473958490
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_7_D +0.003472222222222222222175164840
# define RVV_TGAMMAF_SOURCE_GAMMA_O2_8_D +0.083333333333333333335592087900
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_0_F +0.9435157776e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_1_F +0.8670063615e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_2_F +0.4826702476e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_3_F -0.8855129778e-1f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_4_F +0.1013825238e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_5_F -0.1493408978e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_6_F +0.1697509140e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_7_F -0.2072454542e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O0_8_F +0.2705872357e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_0_F +0.1102489550e-3f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_1_F +0.8160019934e-4f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_2_F +0.1528468856e-3f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_3_F -0.2355068718e-3f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_4_F +0.4962242092e-3f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_5_F -0.1193488017e-2f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_6_F +0.2891599433e-2f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_7_F -0.7385451812e-2f
# define RVV_TGAMMAF_SOURCE_GAMMA_O1_8_F +0.2058077045e-1f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_0_0_F -0.400686534596170958447352690395e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_0_1_F -0.673523028297382446749257758235e-1f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_1_0_F +0.822466960142643054450325495997e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_1_1_F +0.322467033928981157743538726901e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_2_0_F -0.577215665946766039837398973297e+0f
# define RVV_TGAMMAF_SOURCE_GAMMA_C_2_1_F +0.422784335087484338986941629852e+0f
# define RVV_TGAMMAF_LOG_L2U_F RVV_LOG1PF_LN2_HI_F
# define RVV_TGAMMAF_LOG_L2L_F RVV_LOG1PF_LN2_LO_F
# define RVV_TGAMMAF_VILOGB_LIMIT_F RVV_SLEEF_CBRTF_VILOGBK_THRESH
# define RVV_TGAMMAF_VILOGB_SCALE_F RVV_SLEEF_CBRTF_VILOGBK_SCALE
# define RVV_TGAMMAF_TINY_F 1e-30f
# define RVV_TGAMMAF_LOG_SCALE_F (1.0f / 0.75f)
# define RVV_TGAMMAF_SIN_SCALE_F (1.0f / 4096.0f)
# define RVV_TGAMMAF_SIN_UNSCALE_F 4096.0f
# define RVV_TGAMMAF_TINY_SCALE_F 0x1.0p+60f
# define RVV_TGAMMAF_HALF_F RVV_SLEEF_EXPM1F_HALF
# define RVV_TGAMMAF_LOG_SQRT_2PI_D 0.91893853320467278056
# define RVV_TGAMMAF_LOG_PI_D 1.1447298858494001639
# define RVV_TGAMMAF_LOG2_60_D 41.58883083359671856503
# define RVV_TGAMMAF_GAMMA_O0_HI_F 1.2f
# define RVV_TGAMMAF_GAMMA_O2_LO_F 2.3f
# define RVV_TGAMMAF_GAMMA_O2_HI_F 7.0f
#endif

/* Scalar constants for the single-precision powf kernel.  The log2/exp2
   lookup tables remain indexed data in v_f_powf_data.c.  */
#ifndef RVV_POWF_CONSTANTS
# define RVV_POWF_CONSTANTS
# define RVV_POWF_LOG2_OFF 0x3f330000u
# define RVV_POWF_LOG2_TOP_MASK 0xff800000u
# define RVV_POWF_MIN_NORM 0x00800000u
# define RVV_POWF_SIGN_MASK 0x80000000u
# define RVV_POWF_ABS_MASK 0x7fffffffu
# define RVV_POWF_ONE_BITS 0x3f800000u
# define RVV_POWF_SUBNORM_EXPONENT_BIAS (23u << 23)
# define RVV_POWF_LOG2_TABLE_MASK 15u
# define RVV_POWF_EXP2_TABLE_MASK 31u
# define RVV_POWF_FLOAT_MANTISSA_BITS 23u
# define RVV_POWF_DOUBLE_MANTISSA_BITS 52u
# define RVV_POWF_LOG2_TABLE_BITS 4u
# define RVV_POWF_EXP2_TABLE_BITS 5u
# define RVV_POWF_LOG2_POLY0_D 0x1.27616c9496e0bp-2
# define RVV_POWF_LOG2_POLY1_D -0x1.71969a075c67ap-2
# define RVV_POWF_LOG2_POLY2_D 0x1.ec70a6ca7baddp-2
# define RVV_POWF_LOG2_POLY3_D -0x1.7154748bef6c8p-1
# define RVV_POWF_LOG2_POLY4_D 0x1.71547652ab82bp0
# define RVV_POWF_INT_BOUND_F 16777216.0f
# define RVV_POWF_SUBNORM_SCALE_F 0x1.0p23f
# define RVV_POWF_MIN_SUBNORMAL_F RVV_SLEEF_EXPM1F_MIN_SUBNORMAL
# define RVV_POWF_MAX_FINITE_F RVV_SLEEF_EXPM1F_MAX_FINITE
# define RVV_POWF_EXP2_SHIFT_D 0x1.8p+47
# define RVV_POWF_EXP2_C0_D 0x1.c6af84b912394p-5
# define RVV_POWF_EXP2_C1_D 0x1.ebfce50fac4f3p-3
# define RVV_POWF_EXP2_C2_D 0x1.62e42ff0c52d6p-1
# define RVV_POWF_YLOGX_UNDER_D -150.0
# define RVV_POWF_YLOGX_OVER_D 0x1.fffffffd1d571p+6
#endif

#define EXCEPTION_HANDLING_TRIG(vx, expo_x, special_args, vy_special, vlen)  \
  do                                                                          \
    {                                                                         \
      VUINT vclass = __riscv_vfclass ((vx), (vlen));                          \
      VBOOL nan_inf;                                                          \
      IDENTIFY (vclass, class_NaN | class_Inf, nan_inf, (vlen));             \
      VBOOL small_x                                                           \
        = __riscv_vmsleu ((expo_x), EXP_BIAS - MAN_LEN - 5, (vlen));         \
      (special_args) = __riscv_vmor (nan_inf, small_x, (vlen));               \
      if (__riscv_vcpop ((special_args), (vlen)) > 0)                         \
        {                                                                     \
          VBOOL id_mask;                                                      \
          IDENTIFY (vclass, class_Inf, id_mask, (vlen));                      \
          (vy_special)                                                        \
            = __riscv_vfmadd (small_x, (vx), 0x1.0p-60, (vx), (vlen));        \
          VFLOAT vy_nan_inf                                                   \
            = __riscv_vfmerge ((vx), fp_sNaN, id_mask, (vlen));               \
          vy_nan_inf                                                          \
            = __riscv_vfadd (nan_inf, vy_nan_inf, vy_nan_inf, (vlen));        \
          (vy_special)                                                        \
            = __riscv_vmerge ((vy_special), vy_nan_inf, nan_inf, (vlen));     \
          (vx) = __riscv_vfmerge ((vx), fp_posZero, (special_args), (vlen));  \
        }                                                                     \
    }                                                                         \
  while (0)

/* The 2^500 * (2/pi) table is owned by v_d_cos_data.c.  */

/* Trigonometric argument reduction for |x| >= 2^24.  */
#define LARGE_ARGUMENT_REDUCTION_Piby2(vx, vlen, x_large, n_xlarge, r_xlarge, \
                                       r_delta_xlarge)                        \
  do                                                                          \
    {                                                                         \
      VBOOL x_not_large = __riscv_vmnot ((x_large), vlen);                    \
      VFLOAT vx_work = __riscv_vfmerge ((vx), 0x1.0p30, x_not_large, vlen);   \
      VINT lsb_x = U_AS_I (__riscv_vsrl (F_AS_U (vx_work), 52, vlen));        \
      lsb_x = __riscv_vand (lsb_x, 0x7ff, vlen);                              \
      lsb_x = __riscv_vsub (lsb_x, 1075, vlen);                               \
      VUINT expo_mask = VMVU_VX (0x7ff, vlen);                                \
      expo_mask = __riscv_vsll (expo_mask, 52, vlen);                         \
      VFLOAT scale = U_AS_F (__riscv_vand (F_AS_U (vx_work), expo_mask,       \
                                           vlen));                            \
      scale = __riscv_vfmul (scale, 0x1.0p-500, vlen);                        \
      expo_mask = __riscv_vnot (expo_mask, vlen);                             \
      VUINT expo_1 = VMVU_VX (0x3ff, vlen);                                   \
      expo_1 = __riscv_vsll (expo_1, 52, vlen);                               \
      vx_work = U_AS_F (__riscv_vand (F_AS_U (vx_work), expo_mask, vlen));    \
      vx_work = U_AS_F (__riscv_vor (F_AS_U (vx_work), expo_1, vlen));        \
      VUINT j_start = I_AS_U (__riscv_vmax (lsb_x, 2, vlen));                 \
      j_start = __riscv_vsub (j_start, 2, vlen);                              \
      j_start = __riscv_vmul (j_start, 20165, vlen);                          \
      j_start = __riscv_vsrl (j_start, 20, vlen);                             \
      VUINT ind = __riscv_vsll (j_start, 3, vlen);                            \
      VFLOAT two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind,    \
                                           vlen);                             \
      VUINT peg_expo = VMVU_VX (1524, vlen);                                  \
      peg_expo = __riscv_vnmsac (peg_expo, 52, j_start, vlen);                \
      VFLOAT peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      VFLOAT s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);              \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      VFLOAT ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);               \
      VFLOAT prod_0 = s;                                                      \
      VFLOAT prod_1 = ss;                                                     \
      prod_0 = __riscv_vfmul (prod_0, scale, vlen);                           \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_1 = __riscv_vfadd (prod_1, s, vlen);                               \
      VFLOAT prod_2 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_1 = __riscv_vfmul (prod_1, scale, vlen);                           \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_2 = __riscv_vfadd (prod_2, s, vlen);                               \
      VFLOAT prod_3 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_2 = __riscv_vfmul (prod_2, scale, vlen);                           \
      VBOOL ignore_prod_0 = __riscv_vmfge (scale, 0x1.0p-446, vlen);          \
      prod_0 = __riscv_vfmerge (prod_0, 0.0, ignore_prod_0, vlen);            \
      VFLOAT flt_n = __riscv_vfmul (prod_0, 0x1.0p-12, vlen);                 \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_0 = __riscv_vfnmsac (prod_0, 0x1.0p12, flt_n, vlen);               \
      flt_n = __riscv_vfmul (prod_1, 0x1.0p-12, vlen);                        \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_1 = __riscv_vfnmsac (prod_1, 0x1.0p12, flt_n, vlen);               \
      flt_n = __riscv_vfadd (prod_1, prod_2, vlen);                           \
      flt_n = __riscv_vfadd (flt_n, prod_0, vlen);                            \
      (n_xlarge) = __riscv_vfcvt_x (flt_n, vlen);                             \
      flt_n = __riscv_vfcvt_f ((n_xlarge), vlen);                             \
      prod_0 = __riscv_vfsub (prod_0, flt_n, vlen);                           \
      VFLOAT r_hi_local = __riscv_vfadd (prod_0, prod_1, vlen);               \
      VFLOAT r_lo_local = __riscv_vfsub (prod_0, r_hi_local, vlen);           \
      r_lo_local = __riscv_vfadd (r_lo_local, prod_1, vlen);                  \
      VFLOAT tmp_1 = __riscv_vfadd (r_hi_local, prod_2, vlen);                \
      VFLOAT tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                 \
      tmp_2 = __riscv_vfadd (tmp_2, prod_2, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      ss = __riscv_vfmsub (vx_work, two_by_pi, s, vlen);                      \
      prod_3 = __riscv_vfadd (prod_3, s, vlen);                               \
      VFLOAT prod_4 = I_AS_F (__riscv_vor (F_AS_I (ss), F_AS_I (ss), vlen));  \
      prod_3 = __riscv_vfmul (prod_3, scale, vlen);                           \
      tmp_1 = __riscv_vfadd (r_hi_local, prod_3, vlen);                       \
      tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                        \
      tmp_2 = __riscv_vfadd (tmp_2, prod_3, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      ind = __riscv_vadd (ind, 8, vlen);                                      \
      two_by_pi = __riscv_vluxei64 (__v_cos_data.two_over_pi, ind, vlen);     \
      peg_expo = __riscv_vsub (peg_expo, 52, vlen);                           \
      peg = U_AS_F (__riscv_vsll (peg_expo, 52, vlen));                       \
      peg = __riscv_vfsgnj (peg, two_by_pi, vlen);                            \
      peg = __riscv_vfsgnjx (peg, vx_work, vlen);                             \
      s = __riscv_vfmadd (vx_work, two_by_pi, peg, vlen);                     \
      s = __riscv_vfsub (s, peg, vlen);                                       \
      prod_4 = __riscv_vfadd (prod_4, s, vlen);                               \
      prod_4 = __riscv_vfmul (prod_4, scale, vlen);                           \
      tmp_1 = __riscv_vfadd (r_hi_local, prod_4, vlen);                       \
      tmp_2 = __riscv_vfsub (r_hi_local, tmp_1, vlen);                        \
      tmp_2 = __riscv_vfadd (tmp_2, prod_4, vlen);                            \
      r_hi_local = tmp_1;                                                     \
      r_lo_local = __riscv_vfadd (r_lo_local, tmp_2, vlen);                   \
      (r_xlarge) = __riscv_vfmul (r_hi_local, PIBY2_HI, vlen);                \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmsub (r_hi_local, PIBY2_HI, (r_xlarge), vlen);            \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmacc ((r_delta_xlarge), PIBY2_MID, r_hi_local, vlen);     \
      (r_delta_xlarge)                                                        \
        = __riscv_vfmacc ((r_delta_xlarge), PIBY2_HI, r_lo_local, vlen);      \
    }                                                                         \
  while (0)

#endif /* RVV_LIBMVEC_TRIG_DOUBLE_HELPERS */

/* Shared importer-qualified helper names for the single-precision erfcf
   source-shaped finite graphs.  The helper bodies stay in the two erfcf
   skeletons; only their LMUL/SIMDLEN naming and the source RNE conversion
   spelling are common here.  */
#undef RVV_ERFCF_HELPER_NAME2
#undef RVV_ERFCF_HELPER_NAME
#undef RVV_ERFCF_HELPER
#define RVV_ERFCF_HELPER_NAME2(name, lmul, simdlen) \
  rvv_erfcf_##name##_##lmul##_##simdlen
#define RVV_ERFCF_HELPER_NAME(name, lmul, simdlen) \
  RVV_ERFCF_HELPER_NAME2 (name, lmul, simdlen)
#define RVV_ERFCF_HELPER(name) \
  RVV_ERFCF_HELPER_NAME (name, LMUL, SIMDLEN)

#undef RVV_SLEEF_ORIG_ERFCF_VFCVT_X_RNE
#undef RVV_SLEEF_ORIG_ERFCF_HELPER_NAME2
#undef RVV_SLEEF_ORIG_ERFCF_HELPER_NAME
#undef RVV_SLEEF_ORIG_ERFCF_HELPER
#define RVV_SLEEF_ORIG_ERFCF_VFCVT_X_RNE(value, vlen) \
  RVV_VFCVT_X_F_RNE ((value), (vlen))
#define RVV_SLEEF_ORIG_ERFCF_HELPER_NAME2(name, lmul, simdlen) \
  rvv_sleef_orig_erfcf_##name##_##lmul##_##simdlen
#define RVV_SLEEF_ORIG_ERFCF_HELPER_NAME(name, lmul, simdlen) \
  RVV_SLEEF_ORIG_ERFCF_HELPER_NAME2 (name, lmul, simdlen)
#define RVV_SLEEF_ORIG_ERFCF_HELPER(name) \
  RVV_SLEEF_ORIG_ERFCF_HELPER_NAME (name, LMUL, SIMDLEN)

/* Shared importer-qualified helper names for the single-precision erff
   source-shaped finite graphs.  The helper bodies remain in the two erff
   skeletons; only their LMUL/SIMDLEN naming and source RNE conversion are
   common here.  */
#undef RVV_ERFF_POLY2DF_VF2
#undef RVV_ERFF_POLY2DF_B_VF2
#undef RVV_ERFF_POLY4DF_VF2
#undef RVV_ERFF_EXPKF_VF2
#undef RVV_ERFF_VLDEXP_VF
#define RVV_ERFF_POLY2DF_VF2 \
  __PASTE4 (erff_poly2df_vf2_, LMUL, _, SIMDLEN)
#define RVV_ERFF_POLY2DF_B_VF2 \
  __PASTE4 (erff_poly2df_b_vf2_, LMUL, _, SIMDLEN)
#define RVV_ERFF_POLY4DF_VF2 \
  __PASTE4 (erff_poly4df_vf2_, LMUL, _, SIMDLEN)
#define RVV_ERFF_EXPKF_VF2 \
  __PASTE4 (erff_expkf_vf2_, LMUL, _, SIMDLEN)
#define RVV_ERFF_VLDEXP_VF \
  __PASTE4 (erff_vldexp_vf_, LMUL, _, SIMDLEN)

#undef RVV_SLEEF_ORIG_ERFF_VFCVT_X_RNE
#undef RVV_SLEEF_ORIG_ERFF_HELPER_NAME2
#undef RVV_SLEEF_ORIG_ERFF_HELPER_NAME
#undef RVV_SLEEF_ORIG_ERFF_HELPER
#define RVV_SLEEF_ORIG_ERFF_VFCVT_X_RNE(value, vlen) \
  RVV_VFCVT_X_F_RNE ((value), (vlen))
#define RVV_SLEEF_ORIG_ERFF_HELPER_NAME2(name, lmul, simdlen) \
  rvv_sleef_orig_erff_##name##_##lmul##_##simdlen
#define RVV_SLEEF_ORIG_ERFF_HELPER_NAME(name, lmul, simdlen) \
  RVV_SLEEF_ORIG_ERFF_HELPER_NAME2 (name, lmul, simdlen)
#define RVV_SLEEF_ORIG_ERFF_HELPER(name) \
  RVV_SLEEF_ORIG_ERFF_HELPER_NAME (name, LMUL, SIMDLEN)

/* Shared importer-qualified helper names for the single-precision exp10f
   source-shaped and clean finite graphs.  The helper bodies remain in their
   respective skeletons; only LMUL naming and source RNE spelling are shared
   here.  */
#undef RVV_LIBMVEC_EXP10F_VPOW2I
#undef RVV_LIBMVEC_EXP10F_VLDEXP2
#define RVV_LIBMVEC_EXP10F_VPOW2I \
  __PASTE2 (rvv_libmvec_exp10f_vpow2i_m, LMUL)
#define RVV_LIBMVEC_EXP10F_VLDEXP2 \
  __PASTE2 (rvv_libmvec_exp10f_vldexp2_m, LMUL)

#undef RVV_SLEEF_ORIG_EXP10F_VFCVT_X_RNE
#undef RVV_SLEEF_ORIG_EXP10F_VPOW2I
#undef RVV_SLEEF_ORIG_EXP10F_VLDEXP2
#define RVV_SLEEF_ORIG_EXP10F_VFCVT_X_RNE(value, vlen) \
  RVV_VFCVT_X_F_RNE ((value), (vlen))
#define RVV_SLEEF_ORIG_EXP10F_VPOW2I \
  __PASTE2 (rvv_sleef_orig_exp10f_vpow2i_m, LMUL)
#define RVV_SLEEF_ORIG_EXP10F_VLDEXP2 \
  __PASTE2 (rvv_sleef_orig_exp10f_vldexp2_m, LMUL)

/* Source-shaped single-precision erff uses the uncompensated scalar-plus
   double-float addition spelling from SLEEF.  Keep it separate from the
   compensated clean DFADD2 helper.  */
#ifndef RVV_ERFF_SOURCE_DF_HELPERS
# define RVV_ERFF_SOURCE_DF_HELPERS
# define DFADD_VF2_VF_VF2(a, b_hi, b_lo, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _a = (a);                                                        \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a, _b_hi, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _a, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_b_hi, _v, (vlen)), _b_lo, (vlen));        \
    }                                                                         \
  while (0)
#endif

/* The source-shaped asinf helper is embedded by the asinpif comparator.
   Select its public body from the importer function name so the helper
   skeleton has no routine-local enable/disable macro.  */
#ifndef RVV_MATH_IS_ASINF_HELPERS_ONLY_asinpif_original
# define RVV_MATH_IS_ASINF_HELPERS_ONLY_asinpif_original 1
# define RVV_MATH_IS_ASINF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_ASINF_HELPERS_ONLY_##name
# define RVV_MATH_IS_ASINF_HELPERS_ONLY(name) \
  RVV_MATH_IS_ASINF_HELPERS_ONLY_INNER (name)
#endif

/* The source-shaped cosf helper is embedded by the fastcosf comparator and
   by the clean fastcosf implementation.  Select its public body from the
   importer function name so nested skeletons need no control defines.  */
#ifndef RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf
# define RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf 1
# define RVV_MATH_IS_COSF_HELPERS_ONLY_fastcosf_original 1
# define RVV_MATH_IS_COSF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_COSF_HELPERS_ONLY_##name
# define RVV_MATH_IS_COSF_HELPERS_ONLY(name) \
  RVV_MATH_IS_COSF_HELPERS_ONLY_INNER (name)
#endif

/* The fastcosf source skeleton is embedded by the clean fastcosf skeleton.
   Its source-original public body is therefore emitted only for the
   fastcosf_original importer.  */
#ifndef RVV_MATH_IS_FASTCOSF_HELPERS_ONLY_fastcosf
# define RVV_MATH_IS_FASTCOSF_HELPERS_ONLY_fastcosf 1
# define RVV_MATH_IS_FASTCOSF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_FASTCOSF_HELPERS_ONLY_##name
# define RVV_MATH_IS_FASTCOSF_HELPERS_ONLY(name) \
  RVV_MATH_IS_FASTCOSF_HELPERS_ONLY_INNER (name)
#endif

/* The source-shaped log1pf skeleton is embedded by clean log1pf.  Select
   only its finite helper bodies for that importer; the source-original
   importer receives the public comparison body.  */
#ifndef RVV_MATH_IS_LOG1PF_HELPERS_ONLY_log1pf
# define RVV_MATH_IS_LOG1PF_HELPERS_ONLY_log1pf 1
# define RVV_MATH_IS_LOG1PF_HELPERS_ONLY_log1pf_original 0
# define RVV_MATH_IS_LOG1PF_HELPERS_ONLY_INNER(name) \
  RVV_MATH_IS_LOG1PF_HELPERS_ONLY_##name
# define RVV_MATH_IS_LOG1PF_HELPERS_ONLY(name) \
  RVV_MATH_IS_LOG1PF_HELPERS_ONLY_INNER (name)
#endif

/* The private xcosf fallback needs the clean quarantine only when it is used
   as a provider by clean fastcosf.  The source-original fastcosf retains the
   source's late conversion and special-value ordering.  */
#ifndef RVV_MATH_IS_FASTCOSF_CLEAN_EH_fastcosf
# define RVV_MATH_IS_FASTCOSF_CLEAN_EH_fastcosf 1
# define RVV_MATH_IS_FASTCOSF_CLEAN_EH_INNER(name) \
  RVV_MATH_IS_FASTCOSF_CLEAN_EH_##name
# define RVV_MATH_IS_FASTCOSF_CLEAN_EH(name) \
  RVV_MATH_IS_FASTCOSF_CLEAN_EH_INNER (name)
#endif

/* Shared source-shaped acosf names and constants.  The clean finite kernel
   and the independent SLEEF comparator use distinct prefixes, but both
   source operation graphs are instantiated by the same importer LMUL.  */
#ifndef RVV_ACOSF_COMMON_DEFINITIONS
# define RVV_ACOSF_COMMON_DEFINITIONS
# define RVV_ACOSF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ACOSF_PIBY2_HI_F 0x1.921fb6p+0f
# define RVV_ACOSF_PIBY2_LO_F -0x1.777a5cp-25f
# define RVV_ACOSF_PI_HI_F RVV_SLEEF_TRIG_PI_F
# define RVV_ACOSF_PI_LO_F -0x1.777a5cp-24f
# define RVV_ACOSF_C0 0.1666677296e+0f
# define RVV_ACOSF_C1 0.7495029271e-1f
# define RVV_ACOSF_C2 0.4547423869e-1f
# define RVV_ACOSF_C3 0.2424046025e-1f
# define RVV_ACOSF_C4 0.4197454825e-1f

/* The pinned SLEEF source passes the five polynomial literals as double
   expressions to vcast_vf_f(float).  Keep them unsuffixed here so the
   conversion remains observable under a non-RNE caller FRM.  */
# define RVV_SLEEF_ORIG_ACOSF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_SLEEF_ORIG_ACOSF_C0 0.1666677296e+0
# define RVV_SLEEF_ORIG_ACOSF_C1 0.7495029271e-1
# define RVV_SLEEF_ORIG_ACOSF_C2 0.4547423869e-1
# define RVV_SLEEF_ORIG_ACOSF_C3 0.2424046025e-1
# define RVV_SLEEF_ORIG_ACOSF_C4 0.4197454825e-1
# define RVV_SLEEF_ORIG_ACOSF_PI2_HI 3.1415927410125732422f / 2
# define RVV_SLEEF_ORIG_ACOSF_PI2_LO -8.7422776573475857731e-08f / 2
# define RVV_SLEEF_ORIG_ACOSF_PI_HI 3.1415927410125732422f
# define RVV_SLEEF_ORIG_ACOSF_PI_LO -8.7422776573475857731e-08f
#endif

/* Shared source-shaped asinf names and constants.  The clean importer and
   the independent SLEEF comparator use separate prefixes for their helper
   symbols, while both mangle through the current importer LMUL.  The
   `_original` polynomial literals intentionally remain unsuffixed: pinned
   SLEEF passes them as double expressions to vcast_vf_f(float), so their
   conversion remains source-shaped under a non-RNE caller FRM.  */
#ifndef RVV_ASINF_COMMON_DEFINITIONS
# define RVV_ASINF_COMMON_DEFINITIONS
# define RVV_ASINF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ASINF_PIBY4_HI 3.1415927410125732422f / 4.0f
# define RVV_ASINF_PIBY4_LO -8.7422776573475857731e-08f / 4.0f
# define RVV_ASINF_C0 0.1666677296e+0f
# define RVV_ASINF_C1 0.7495029271e-1f
# define RVV_ASINF_C2 0.4547423869e-1f
# define RVV_ASINF_C3 0.2424046025e-1f
# define RVV_ASINF_C4 0.4197454825e-1f
# define RVV_SLEEF_ORIG_ASINF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_SLEEF_ORIG_ASINF_C0 0.1666677296e+0
# define RVV_SLEEF_ORIG_ASINF_C1 0.7495029271e-1
# define RVV_SLEEF_ORIG_ASINF_C2 0.4547423869e-1
# define RVV_SLEEF_ORIG_ASINF_C3 0.2424046025e-1
# define RVV_SLEEF_ORIG_ASINF_C4 0.4197454825e-1
# define RVV_SLEEF_ORIG_ASINF_PI4_HI 3.1415927410125732422f / 4.0f
# define RVV_SLEEF_ORIG_ASINF_PI4_LO \
  -8.7422776573475857731e-08f / 4.0f
#endif

/* Source-shaped helper names used by the asinpif composite.  Keep this
   importer-dependent mangle in the common math header; the routine skeleton
   itself should contain only the finite operation graph and data references.  */
#ifndef RVV_ASINPIF_COMMON_DEFINITIONS
# define RVV_ASINPIF_COMMON_DEFINITIONS
# define RVV_ASINPIF_SOURCE_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
# define RVV_ASINPIF_ONE_OVER_PI_HI_F RVV_SLEEF_TRIG_ONE_OVER_PI
# define RVV_ASINPIF_ONE_OVER_PI_LO_F 0x1.b93910p-27f
# define RVV_ASINPIF_PI4_HI_F RVV_ASINF_PIBY4_HI
# define RVV_ASINPIF_PI4_LO_F RVV_ASINF_PIBY4_LO
# define RVV_ASINPIF_HALF_F 0.5f
# define RVV_ASINPIF_TWO_F 2.0f
# define RVV_ASINPIF_TINY_SQUARE_LIMIT_F 0x1p-63f
# define RVV_ASINPIF_POLY0_F 0.1666677296e+0f
# define RVV_ASINPIF_POLY1_F 0.7495029271e-1f
# define RVV_ASINPIF_POLY2_F 0.4547423869e-1f
# define RVV_ASINPIF_POLY3_F 0.2424046025e-1f
# define RVV_ASINPIF_POLY4_F 0.4197454825e-1f
#endif

/* Target-owned source-shaped SLEEF DF/FMA helpers.  These macros are
   kept in the common math header so skeletons do not carry a private
   helper header.  The operation graphs are unchanged from the former
   former private SLEEF helper translation unit.  */
#ifndef RVV_SLEEF_ORIGINAL_DF_HELPERS
#define RVV_SLEEF_ORIGINAL_DF_HELPERS

#define RVV_SLEEF_ORIG_MANGLE_INNER(base, lmul) base##lmul
#define RVV_SLEEF_ORIG_MANGLE(base, lmul) \
  RVV_SLEEF_ORIG_MANGLE_INNER (base, lmul)

#define RVV_SLEEF_ORIG_DFNORMALIZE_VF2(a_hi, a_lo, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _a_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _a_lo, (vlen));       \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFNEG_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      (r_hi) = VFNEG_V (_a_hi, (vlen));                                      \
      (r_lo) = VFNEG_V (_a_lo, (vlen));                                      \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF_VF(x, y, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x, _y, (vlen));                                  \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                      \
          VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)),                   \
          VFSUB_VV (_y, _v, (vlen)), (vlen));                                 \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD_VF_VF(x, y, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x, _y, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_x, _s, (vlen)), _y, (vlen));              \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)),                   \
          VFSUB_VV (_y_hi, _v, (vlen)), (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)         \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x_hi, _y, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_y, _v, (vlen)), (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _x_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD2_VF2_VF_S(x_hi, x_lo, y, r_hi, r_lo, vlen)       \
  RVV_SLEEF_ORIG_DFADD2_VF2_VF ((x_hi), (x_lo), VFMV_VF ((y), (vlen)),       \
                                (r_hi), (r_lo), (vlen))

#define RVV_SLEEF_ORIG_DFADD_VF2_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)      \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFADD_VV (_x_hi, _y, (vlen));                               \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x_hi, _s, (vlen)), _y, (vlen));        \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _x_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x, _s, (vlen)), _y_hi, (vlen));        \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFADD_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x_hi, _s, (vlen)), _y_hi, (vlen));     \
      _t = VFADD_VV (_t, _x_lo, (vlen));                                      \
      _t = VFADD_VV (_t, _y_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

/* Source dfadd2_vf2_vf2_vf2: compensated pair-plus-pair addition.  */
#define RVV_SLEEF_ORIG_DFADD2_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,     \
                                      vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _v = VFSUB_VV (_s, _x_hi, (vlen));                               \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_x_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_y_hi, _v, (vlen)), (vlen));                              \
      VFLOAT _lo = VFADD_VV (_x_lo, _y_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _lo, (vlen));                                    \
    }                                                                         \
  while (0)

/* Source dfsub_vf2_vf2_vf2: pair subtraction with the source's sequential
   high-part subtraction and low-part correction order.  */
#define RVV_SLEEF_ORIG_DFSUB_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFSUB_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _t = VFSUB_VV (_x_hi, _s, (vlen));                               \
      _t = VFSUB_VV (_t, _y_hi, (vlen));                                      \
      _t = VFADD_VV (_t, _x_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFSUB_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFSCALE_VF2_VF(x_hi, x_lo, scale, r_hi, r_lo, vlen)   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _scale = (scale);                                                \
      (r_hi) = VFMUL_VV (_x_hi, _scale, (vlen));                              \
      (r_lo) = VFMUL_VV (_x_lo, _scale, (vlen));                              \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFMUL_VF_VF(x, y, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFMUL_VV (_x, _y, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMSUB_VV (_x, _y, _s, (vlen));                                \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFMUL_VF2_VF(x_hi, x_lo, y, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y = (y);                                                        \
      VFLOAT _s = VFMUL_VV (_x_hi, _y, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_x_lo, _y, VFMSUB_VV (_x_hi, _y, _s, (vlen)),       \
                          (vlen));                                           \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFMUL_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r_hi, r_lo,      \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_x_hi, _y_hi, (vlen));                            \
      VFLOAT _residual = VFMSUB_VV (_x_hi, _y_hi, _s, (vlen));                \
      VFLOAT _cross_hi = VFMADD_VV (_x_lo, _y_hi, _residual, (vlen));         \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_x_hi, _y_lo, _cross_hi, (vlen));                   \
    }                                                                         \
  while (0)

/* Source dfmul_vf_vf2_vf2: return the collapsed scalar product of two
   double-float pairs with the source's nested FMA association.  */
#define RVV_SLEEF_ORIG_DFMUL_VF2_VF2_SCALAR(x_hi, x_lo, y_hi, y_lo, result,   \
                                            vlen)                              \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _low = VFMADD_VV (_x_lo, _y_hi,                                 \
                               VFMUL_VV (_x_hi, _y_lo, (vlen)), (vlen));      \
      (result) = VFMADD_VV (_x_hi, _y_hi, _low, (vlen));                      \
    }                                                                         \
  while (0)

/* Source dfrec_vf2_vf in the FMA branch: reciprocal correction for a scalar
   high part retains the source's nested negative-product sequence.  */
#define RVV_SLEEF_ORIG_DFREC_VF2_VF(x, r_hi, r_lo, vlen)                     \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _x, (vlen));              \
      VFLOAT _u = VFNMSUB_VV (_x, _t, VFMV_VF (1.0f, (vlen)), (vlen));         \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_t, _u, (vlen));                                     \
    }                                                                         \
  while (0)

/* Source dfrec_vf2_vf2 in the FMA branch: reciprocal correction retains the
   source's nested negative-product sequence.  */
#define RVV_SLEEF_ORIG_DFREC_VF2_VF2(x_hi, x_lo, r_hi, r_lo, vlen)            \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _x_hi, (vlen));            \
      VFLOAT _u = VFNMSUB_VV (_x_hi, _t, VFMV_VF (1.0f, (vlen)), (vlen));      \
      _u = VFNMSAC_VV (_u, _x_lo, _t, (vlen));                                \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_t, _u, (vlen));                                     \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFDIV_VF2_VF2(n_hi, n_lo, d_hi, d_lo, r_hi, r_lo,     \
                                     vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _n_hi = (n_hi);                                                  \
      VFLOAT _n_lo = (n_lo);                                                  \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _s = VFMUL_VV (_n_hi, _t, (vlen));                               \
      VFLOAT _u = VFMSUB_VV (_t, _n_hi, _s, (vlen));                          \
      VFLOAT _v = VFNMSUB_VV (_d_hi, _t, VFMV_VF (1.0f, (vlen)), (vlen));     \
      _v = VFNMSUB_VV (_d_lo, _t, _v, (vlen));                                \
      VFLOAT _nlo_t = VFMADD_VV (_n_lo, _t, _u, (vlen));                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (_s, _v, _nlo_t, (vlen));                            \
    }                                                                         \
  while (0)

#define RVV_SLEEF_ORIG_DFSQU_VF2(x_hi, x_lo, r_hi, r_lo, vlen)                \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                  \
      VFLOAT _x_lo = (x_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_x_hi, _x_hi, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMADD_VV (VFADD_VV (_x_hi, _x_hi, (vlen)), _x_lo,             \
                          VFMSUB_VV (_x_hi, _x_hi, _s, (vlen)), (vlen));      \
    }                                                                         \
  while (0)

/* Source dfsqrt_vf2_vf2 in common/df.h.  The pinned RVV source does not
   define ENABLE_RECSQRT_SP, so this is the source's sqrt/reciprocal branch:
   sqrt(high+low), square it as a source DF pair, add it to the input pair,
   multiply by the source DF reciprocal, then scale by one half.  */
#define RVV_SLEEF_ORIG_DFSQRT_VF2_VF2(d_hi, d_lo, r_hi, r_lo, vlen)           \
  do                                                                          \
    {                                                                         \
      VFLOAT _d_hi = (d_hi);                                                   \
      VFLOAT _d_lo = (d_lo);                                                   \
      VFLOAT _t = __riscv_vfsqrt (VFADD_VV (_d_hi, _d_lo, (vlen)), (vlen));    \
      VFLOAT _tsq_hi, _tsq_lo;                                                \
      RVV_SLEEF_ORIG_DFMUL_VF_VF (_t, _t, _tsq_hi, _tsq_lo, (vlen));           \
      VFLOAT _sum_hi, _sum_lo;                                                \
      RVV_SLEEF_ORIG_DFADD2_VF2_VF2                                          \
        (_d_hi, _d_lo, _tsq_hi, _tsq_lo, _sum_hi, _sum_lo, (vlen));            \
      VFLOAT _rec_hi, _rec_lo;                                                \
      RVV_SLEEF_ORIG_DFREC_VF2_VF (_t, _rec_hi, _rec_lo, (vlen));              \
      VFLOAT _root_hi, _root_lo;                                              \
      RVV_SLEEF_ORIG_DFMUL_VF2_VF2                                           \
        (_sum_hi, _sum_lo, _rec_hi, _rec_lo, _root_hi, _root_lo, (vlen));      \
      RVV_SLEEF_ORIG_DFSCALE_VF2_VF                                           \
        (_root_hi, _root_lo, VFMV_VF (0.5f, (vlen)),                          \
         (r_hi), (r_lo), (vlen));                                             \
    }                                                                         \
  while (0)

#endif /* RVV_SLEEF_ORIGINAL_DF_HELPERS */

/* Shared source-original gamma helper names and widened clean-repair vector
   types.  The public tgammaf entry is not enabled by this block; lgammaf uses
   only the gammafk/logk2f subgraph below.  */
#undef RVV_SLEEF_ORIG_TGAMMAF_MANGLE
#define RVV_SLEEF_ORIG_TGAMMAF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)

#undef RVV_LGAMMAF_REPAIR_MANGLE_INNER
#undef RVV_LGAMMAF_REPAIR_MANGLE
#undef RVV_LGAMMAF_REPAIR_CALL_INNER
#undef RVV_LGAMMAF_REPAIR_CALL
#define RVV_LGAMMAF_REPAIR_MANGLE_INNER(base, lmul) base##lmul
#define RVV_LGAMMAF_REPAIR_MANGLE(base, lmul) \
  RVV_LGAMMAF_REPAIR_MANGLE_INNER (base, lmul)
#define RVV_LGAMMAF_REPAIR_CALL_INNER(base, lmul) base##lmul
#define RVV_LGAMMAF_REPAIR_CALL(base, lmul) \
  RVV_LGAMMAF_REPAIR_CALL_INNER (base, lmul)

#undef RVV_LGAMMAF_REPAIR_D
#undef RVV_LGAMMAF_REPAIR_U
#undef RVV_LGAMMAF_REPAIR_I
#undef RVV_LGAMMAF_REPAIR_B
#undef RVV_LGAMMAF_REPAIR_VL
#undef RVV_LGAMMAF_REPAIR_NC
#undef RVV_LGAMMAF_REPAIR_MV
#undef RVV_LGAMMAF_REPAIR_MERGE
#undef RVV_LGAMMAF_REPAIR_MERGE_X
#undef RVV_LGAMMAF_REPAIR_F2U
#undef RVV_LGAMMAF_REPAIR_U2F
#undef RVV_LGAMMAF_REPAIR_U2I
#if (BIT_WIDTH == 32) && (LMUL == 1)
# define RVV_LGAMMAF_REPAIR_D vfloat64m2_t
# define RVV_LGAMMAF_REPAIR_U vuint64m2_t
# define RVV_LGAMMAF_REPAIR_I vint64m2_t
# define RVV_LGAMMAF_REPAIR_B vbool32_t
# define RVV_LGAMMAF_REPAIR_VL __riscv_vfwcvt_f_f_v_f64m2
# define RVV_LGAMMAF_REPAIR_NC __riscv_vfncvt_f_f_w_f32m1
# define RVV_LGAMMAF_REPAIR_MV __riscv_vfmv_v_f_f64m2
# define RVV_LGAMMAF_REPAIR_MERGE __riscv_vmerge_vvm_f64m2
# define RVV_LGAMMAF_REPAIR_MERGE_X __riscv_vmerge_vxm_i64m2
# define RVV_LGAMMAF_REPAIR_F2U __riscv_vreinterpret_v_f64m2_u64m2
# define RVV_LGAMMAF_REPAIR_U2F __riscv_vreinterpret_v_u64m2_f64m2
# define RVV_LGAMMAF_REPAIR_U2I __riscv_vreinterpret_v_u64m2_i64m2
#elif (BIT_WIDTH == 32) && (LMUL == 2)
# define RVV_LGAMMAF_REPAIR_D vfloat64m4_t
# define RVV_LGAMMAF_REPAIR_U vuint64m4_t
# define RVV_LGAMMAF_REPAIR_I vint64m4_t
# define RVV_LGAMMAF_REPAIR_B vbool16_t
# define RVV_LGAMMAF_REPAIR_VL __riscv_vfwcvt_f_f_v_f64m4
# define RVV_LGAMMAF_REPAIR_NC __riscv_vfncvt_f_f_w_f32m2
# define RVV_LGAMMAF_REPAIR_MV __riscv_vfmv_v_f_f64m4
# define RVV_LGAMMAF_REPAIR_MERGE __riscv_vmerge_vvm_f64m4
# define RVV_LGAMMAF_REPAIR_MERGE_X __riscv_vmerge_vxm_i64m4
# define RVV_LGAMMAF_REPAIR_F2U __riscv_vreinterpret_v_f64m4_u64m4
# define RVV_LGAMMAF_REPAIR_U2F __riscv_vreinterpret_v_u64m4_f64m4
# define RVV_LGAMMAF_REPAIR_U2I __riscv_vreinterpret_v_u64m4_i64m4
#elif (BIT_WIDTH == 32) && (LMUL == 4)
# define RVV_LGAMMAF_REPAIR_D vfloat64m8_t
# define RVV_LGAMMAF_REPAIR_U vuint64m8_t
# define RVV_LGAMMAF_REPAIR_I vint64m8_t
# define RVV_LGAMMAF_REPAIR_B vbool8_t
# define RVV_LGAMMAF_REPAIR_VL __riscv_vfwcvt_f_f_v_f64m8
# define RVV_LGAMMAF_REPAIR_NC __riscv_vfncvt_f_f_w_f32m4
# define RVV_LGAMMAF_REPAIR_MV __riscv_vfmv_v_f_f64m8
# define RVV_LGAMMAF_REPAIR_MERGE __riscv_vmerge_vvm_f64m8
# define RVV_LGAMMAF_REPAIR_MERGE_X __riscv_vmerge_vxm_i64m8
# define RVV_LGAMMAF_REPAIR_F2U __riscv_vreinterpret_v_f64m8_u64m8
# define RVV_LGAMMAF_REPAIR_U2F __riscv_vreinterpret_v_u64m8_f64m8
# define RVV_LGAMMAF_REPAIR_U2I __riscv_vreinterpret_v_u64m8_i64m8
#endif

/* Source-original SLEEF gamma helpers.  These are shared by tgammaf and the
   lgammaf source-shaped comparator, but the public entries remain in their
   respective skeletons so clean EH and source final-selection policy cannot
   be accidentally coupled.  */
#if (BIT_WIDTH == 32) && RVV_MATH_IS_TGAMMAF_FAMILY (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M1
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M1
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M2
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M2
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M4
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M4
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M8
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_M8
#   define RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL

static inline VINT
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_truncate_m)
  (VFLOAT x, size_t vlen)
{
  return __riscv_vfcvt_rtz_x (x, vlen);
}

static inline VINT
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_vilogbk_m)
  (VFLOAT d, size_t vlen)
{
  VBOOL small = __riscv_vmflt
    (d, RVV_TGAMMAF_VILOGB_LIMIT_F, vlen);
  VFLOAT scaled = VFMUL_VF (d, RVV_TGAMMAF_VILOGB_SCALE_F, vlen);
  d = VMERGE_VVM_F (d, scaled, small, vlen);

  VUINT q = __riscv_vsrl (F_AS_U (d), MAN_LEN_F, vlen);
  q = __riscv_vand (q, 0xff, vlen);
  VINT bias = VMERGE_VVM_I
    (VMVI_VX (0x7f, vlen), VMVI_VX (64 + 0x7f, vlen), small, vlen);
  return __riscv_vsub (U_AS_I (q), bias, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_vpow2i_m)
  (VINT q, size_t vlen)
{
  VINT biased = __riscv_vadd (q, EXP_BIAS_F, vlen);
  return I_AS_F (__riscv_vsll (biased, MAN_LEN_F, vlen));
}

static inline void
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m)
  (VFLOAT d_hi, VFLOAT d_lo,
   VFLOAT *r_hi, VFLOAT *r_lo, size_t vlen)
{
  VFLOAT d_adj = VFMUL_VF (d_hi, RVV_TGAMMAF_LOG_SCALE_F, vlen);
  VINT e = RVV_SLEEF_ORIG_TGAMMAF_MANGLE
    (rvv_sleef_orig_tgammaf_vilogbk_m) (d_adj, vlen);

  VINT neg_e = __riscv_vneg (e, vlen);
  VFLOAT scale = RVV_SLEEF_ORIG_TGAMMAF_MANGLE
    (rvv_sleef_orig_tgammaf_vpow2i_m) (neg_e, vlen);
  VFLOAT m_hi, m_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (d_hi, d_lo, scale, m_hi, m_lo, vlen);

  VFLOAT num_hi, num_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (m_hi, m_lo, -1.0f, num_hi, num_lo, vlen);
  VFLOAT den_hi, den_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (m_hi, m_lo, 1.0f, den_hi, den_lo, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (num_hi, num_lo, den_hi, den_lo, x_hi, x_lo, vlen);
  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFSQU_VF2 (x_hi, x_lo, x2_hi, x2_lo, vlen);

  VFLOAT poly = VFMV_VF (RVV_TGAMMAF_SOURCE_LOG_POLY_0_F, vlen);
  poly = VFMADD_VV
    (poly, x2_hi, VFMV_VF (RVV_TGAMMAF_SOURCE_LOG_POLY_1_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, x2_hi, VFMV_VF (RVV_TGAMMAF_SOURCE_LOG_POLY_2_F, vlen), vlen);
  poly = VFMADD_VV
    (poly, x2_hi, VFMV_VF (RVV_TGAMMAF_SOURCE_LOG_POLY_3_F, vlen), vlen);

  VFLOAT e_f = __riscv_vfcvt_f (e, vlen);
  VFLOAT s_hi, s_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_TGAMMAF_LOG_L2U_F, vlen),
     VFMV_VF (RVV_TGAMMAF_LOG_L2L_F, vlen), e_f, s_hi, s_lo, vlen);

  VFLOAT two_x_hi, two_x_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (x_hi, x_lo, VFMV_VF (2.0f, vlen), two_x_hi, two_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, two_x_hi, two_x_lo, s_hi, s_lo, vlen);

  VFLOAT x2x_hi, x2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x2_hi, x2_lo, x_hi, x_lo, x2x_hi, x2x_lo, vlen);
  VFLOAT tail_hi, tail_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (x2x_hi, x2x_lo, poly, tail_hi, tail_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (s_hi, s_lo, tail_hi, tail_lo, *r_hi, *r_lo, vlen);
}

static inline void
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_sinpifk_m)
  (VFLOAT d, VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VFLOAT u = VFMUL_VF (d, 4.0f, vlen);
  VINT q = RVV_SLEEF_ORIG_TGAMMAF_MANGLE
    (rvv_sleef_orig_tgammaf_truncate_m) (u, vlen);

  VUINT q_u = I_AS_U (q);
  VUINT q_sign = __riscv_vsrl (q_u, 31, vlen);
  VUINT q_sign_x1 = VXOR_VV (q_sign, VMVU_VX (1, vlen), vlen);
  q = __riscv_vadd (q, U_AS_I (q_sign_x1), vlen);
  q = __riscv_vand (q, VMVI_VX (~1, vlen), vlen);

  VBOOL o = __riscv_vmseq (__riscv_vand (q, 2, vlen), 2, vlen);
  VFLOAT q_f = __riscv_vfcvt_f (q, vlen);
  VFLOAT s = VFSUB_VV (u, q_f, vlen);
  VFLOAT t = s;
  s = VFMUL_VV (s, s, vlen);

  VFLOAT s2_hi, s2_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (t, t, s2_hi, s2_lo, vlen);

  VFLOAT uu = VFMERGE_VFM
    (VFMV_VF (RVV_TGAMMAF_SOURCE_SINPIF_POLY_0_1_F, vlen),
     RVV_TGAMMAF_SOURCE_SINPIF_POLY_0_0_F, o, vlen);
  uu = VFMADD_VV
    (uu, s, VFMERGE_VFM
       (VFMV_VF (RVV_TGAMMAF_SOURCE_SINPIF_POLY_1_1_F, vlen),
        RVV_TGAMMAF_SOURCE_SINPIF_POLY_1_0_F, o, vlen), vlen);
  uu = VFMADD_VV
    (uu, s, VFMERGE_VFM
       (VFMV_VF (RVV_TGAMMAF_SOURCE_SINPIF_POLY_2_1_F, vlen),
        RVV_TGAMMAF_SOURCE_SINPIF_POLY_2_0_F, o, vlen), vlen);

  VFLOAT x_hi, x_lo;
  VFLOAT us = VFMUL_VV (uu, s, vlen);
  VFLOAT c_hi = VFMERGE_VFM
    (VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_C_0_1_D), vlen),
     RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_C_0_0_D), o, vlen);
  VFLOAT c_lo = VFMERGE_VFM
    (VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_C_1_1_D), vlen),
     RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_C_1_0_D), o, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF_VF2
    (us, c_hi, c_lo, x_hi, x_lo, vlen);

  VFLOAT s2x_hi, s2x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s2_hi, s2_lo, x_hi, x_lo, s2x_hi, s2x_lo, vlen);
  VFLOAT p_hi = VFMERGE_VFM
    (VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_P_0_1_D), vlen),
     RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_P_0_0_D), o, vlen);
  VFLOAT p_lo = VFMERGE_VFM
    (VFMV_VF
       (RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_P_1_1_D), vlen),
     RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_SOURCE_SINPIF_P_1_0_D), o, vlen);
  VFLOAT x2_hi, x2_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (s2x_hi, s2x_lo, p_hi, p_lo, x2_hi, x2_lo, vlen);

  VFLOAT m_hi = VMERGE_VVM_F (t, s2_hi, o, vlen);
  VFLOAT m_lo = VMERGE_VVM_F (VFMV_VF (0.0f, vlen), s2_lo, o, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x2_hi, x2_lo, m_hi, m_lo, x_hi, x_lo, vlen);

  VFLOAT one_hi, one_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (x_hi, x_lo, 1.0f, one_hi, one_lo, vlen);
  x_hi = VMERGE_VVM_F (x_hi, one_hi, o, vlen);
  x_lo = VMERGE_VVM_F (x_lo, one_lo, o, vlen);

  VBOOL flip = __riscv_vmseq (__riscv_vand (q, 4, vlen), 4, vlen);
  VINT flip_bits = VMERGE_VVM_I
    (VMVI_VX (0, vlen), VMVI_VX (0x80000000, vlen), flip, vlen);
  VFLOAT flip_f = I_AS_F (flip_bits);
  x_hi = __riscv_vfsgnjx (x_hi, flip_f, vlen);
  x_lo = __riscv_vfsgnjx (x_lo, flip_f, vlen);

  *r_hi = x_hi;
  *r_lo = x_lo;
}

static inline void
RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_gammafk_m)
  (VFLOAT a, VFLOAT *out_clc_hi,
   VFLOAT *out_clc_lo, VFLOAT *out_clln_hi, VFLOAT *out_clln_lo,
   VFLOAT *out_clld_hi, VFLOAT *out_clld_lo, size_t vlen)
{
  VFLOAT clc_hi = VFMV_VF (0.0f, vlen);
  VFLOAT clc_lo = VFMV_VF (0.0f, vlen);
  VFLOAT clln_hi = VFMV_VF (1.0f, vlen);
  VFLOAT clln_lo = VFMV_VF (0.0f, vlen);
  VFLOAT clld_hi = VFMV_VF (1.0f, vlen);
  VFLOAT clld_lo = VFMV_VF (0.0f, vlen);
  VFLOAT x_hi = a;
  VFLOAT x_lo = VFMV_VF (0.0f, vlen);
  VFLOAT y_hi, y_lo, z_hi, z_lo;
  VFLOAT t, u;

  VBOOL otiny = __riscv_vmflt
    (__riscv_vfabs (a, vlen), RVV_TGAMMAF_TINY_F, vlen);
  VBOOL oref = __riscv_vmflt (a, RVV_TGAMMAF_HALF_F, vlen);

  VFLOAT one_minus_a_hi, one_minus_a_lo;
  RVV_SLEEF_ORIG_DFADD2_VF_VF
    (VFMV_VF (1.0f, vlen), VFNEG_V (a, vlen),
     one_minus_a_hi, one_minus_a_lo, vlen);
  x_hi = VMERGE_VVM_F (x_hi, one_minus_a_hi, oref, vlen);
  x_lo = VMERGE_VVM_F (x_lo, one_minus_a_lo, oref, vlen);
  x_hi = VMERGE_VVM_F (x_hi, VFMV_VF (0.0f, vlen), otiny, vlen);
  x_lo = VMERGE_VVM_F (x_lo, VFMV_VF (0.0f, vlen), otiny, vlen);

  VBOOL o0 = __riscv_vmand
    (__riscv_vmfge (x_hi, RVV_TGAMMAF_HALF_F, vlen),
     __riscv_vmfle (x_hi, RVV_TGAMMAF_GAMMA_O0_HI_F, vlen), vlen);
  VBOOL o2 = __riscv_vmfge (x_hi, RVV_TGAMMAF_GAMMA_O2_LO_F, vlen);

  VFLOAT tmp_hi, tmp_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (x_hi, x_lo, 1.0f, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (tmp_hi, tmp_lo, x_hi, x_lo, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (y_hi, y_lo, y_hi, y_lo, vlen);

  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (x_hi, x_lo, 2.0f, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (tmp_hi, tmp_lo, y_hi, y_lo, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (y_hi, y_lo, y_hi, y_lo, vlen);

  VBOOL o = __riscv_vmand
    (o2, __riscv_vmfle (x_hi, RVV_TGAMMAF_GAMMA_O2_HI_F, vlen), vlen);
  clln_hi = VMERGE_VVM_F (clln_hi, y_hi, o, vlen);
  clln_lo = VMERGE_VVM_F (clln_lo, y_lo, o, vlen);

  VFLOAT x3_hi, x3_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (x_hi, x_lo, 3.0f, x3_hi, x3_lo, vlen);
  x_hi = VMERGE_VVM_F (x_hi, x3_hi, o, vlen);
  x_lo = VMERGE_VVM_F (x_lo, x3_lo, o, vlen);

  VFLOAT t1 = VFDIV_VV (VFMV_VF (1.0f, vlen), x_hi, vlen);
  VFLOAT add = VFMERGE_VFM
    (VFMV_VF (-2.0f, vlen), -1.0f, o0, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (x_hi, x_lo, add, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (tmp_hi, tmp_lo, tmp_hi, tmp_lo, vlen);
  t = VMERGE_VVM_F (tmp_hi, t1, o2, vlen);

# define RVV_SLEEF_ORIG_TGAMMAF_COEFF(i)                                  \
  VFMERGE_VFM                                                            \
    (VFMV_VF (__PASTE3 (RVV_TGAMMAF_SOURCE_GAMMA_O1_, i, _F), vlen),      \
     RVV_SLEEF_ORIG_F64_TO_F32                                           \
       (__PASTE3 (RVV_TGAMMAF_SOURCE_GAMMA_O0_, i, _F)), o0, vlen)
# define RVV_SLEEF_ORIG_TGAMMAF_COEFF2(i)                                 \
  VFMERGE_VFM                                                            \
    (RVV_SLEEF_ORIG_TGAMMAF_COEFF (i),                                    \
     RVV_SLEEF_ORIG_F64_TO_F32                                           \
       (__PASTE3 (RVV_TGAMMAF_SOURCE_GAMMA_O2_, i, _D)), o2, vlen)

  u = RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (0);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (1), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (2), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (3), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (4), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (5), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (6), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (7), vlen);
  u = VFMADD_VV (u, t, RVV_SLEEF_ORIG_TGAMMAF_COEFF2 (8), vlen);

# undef RVV_SLEEF_ORIG_TGAMMAF_COEFF2
# undef RVV_SLEEF_ORIG_TGAMMAF_COEFF

  VFLOAT logx_hi, logx_lo;
  RVV_SLEEF_ORIG_TGAMMAF_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m)
    (x_hi, x_lo, &logx_hi, &logx_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (x_hi, x_lo, -0.5f, tmp_hi, tmp_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (tmp_hi, tmp_lo, logx_hi, logx_lo, y_hi, y_lo, vlen);
  VFLOAT neg_x_hi, neg_x_lo;
  RVV_SLEEF_ORIG_DFNEG_VF2
    (x_hi, x_lo, neg_x_hi, neg_x_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (y_hi, y_lo, neg_x_hi, neg_x_lo, y_hi, y_lo, vlen);

  float log_sqrt_hi_s
    = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_LOG_SQRT_2PI_D);
  float log_sqrt_lo_s
    = RVV_SLEEF_ORIG_F64_TO_F32
      (RVV_TGAMMAF_LOG_SQRT_2PI_D - (double) log_sqrt_hi_s);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (y_hi, y_lo, VFMV_VF (log_sqrt_hi_s, vlen),
     VFMV_VF (log_sqrt_lo_s, vlen), y_hi, y_lo, vlen);

  VFLOAT ut_hi, ut_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (u, t, ut_hi, ut_lo, vlen);
  VFLOAT c = VFMERGE_VFM
    (VFMV_VF (RVV_TGAMMAF_SOURCE_GAMMA_C_0_1_F, vlen),
     RVV_TGAMMAF_SOURCE_GAMMA_C_0_0_F, o0, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (ut_hi, ut_lo, c, z_hi, z_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (z_hi, z_lo, t, z_hi, z_lo, vlen);
  c = VFMERGE_VFM
    (VFMV_VF (RVV_TGAMMAF_SOURCE_GAMMA_C_1_1_F, vlen),
     RVV_TGAMMAF_SOURCE_GAMMA_C_1_0_F, o0, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (z_hi, z_lo, c, z_hi, z_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (z_hi, z_lo, t, z_hi, z_lo, vlen);
  c = VFMERGE_VFM
    (VFMV_VF (RVV_TGAMMAF_SOURCE_GAMMA_C_2_1_F, vlen),
     RVV_TGAMMAF_SOURCE_GAMMA_C_2_0_F, o0, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF
    (z_hi, z_lo, c, z_hi, z_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (z_hi, z_lo, t, z_hi, z_lo, vlen);

  clc_hi = VMERGE_VVM_F (z_hi, y_hi, o2, vlen);
  clc_lo = VMERGE_VVM_F (z_lo, y_lo, o2, vlen);

  VFLOAT one_ut_hi, one_ut_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (ut_hi, ut_lo, 1.0f, one_ut_hi, one_ut_lo, vlen);
  clld_hi = VMERGE_VVM_F (clld_hi, one_ut_hi, o2, vlen);
  clld_lo = VMERGE_VVM_F (clld_lo, one_ut_lo, o2, vlen);

  float logpi_hi_s = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_LOG_PI_D);
  float logpi_lo_s = RVV_SLEEF_ORIG_F64_TO_F32
    (RVV_TGAMMAF_LOG_PI_D - (double) logpi_hi_s);
  VFLOAT logpi_hi = VFMV_VF (logpi_hi_s, vlen);
  VFLOAT logpi_lo = VFMV_VF (logpi_lo_s, vlen);
  RVV_SLEEF_ORIG_DFNEG_VF2
    (clc_hi, clc_lo, neg_x_hi, neg_x_lo, vlen);
  VFLOAT clc_ref_hi, clc_ref_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (logpi_hi, logpi_lo, neg_x_hi, neg_x_lo, clc_ref_hi, clc_ref_lo, vlen);
  clc_hi = VMERGE_VVM_F (clc_hi, clc_ref_hi, oref, vlen);
  clc_lo = VMERGE_VVM_F (clc_lo, clc_ref_lo, oref, vlen);

  float log2_60_hi_s
    = RVV_SLEEF_ORIG_F64_TO_F32 (RVV_TGAMMAF_LOG2_60_D);
  float log2_60_lo_s = RVV_SLEEF_ORIG_F64_TO_F32
    (RVV_TGAMMAF_LOG2_60_D - (double) log2_60_hi_s);
  clc_hi = VMERGE_VVM_F
    (clc_hi, VFMV_VF (log2_60_hi_s, vlen), otiny, vlen);
  clc_lo = VMERGE_VVM_F
    (clc_lo, VFMV_VF (log2_60_lo_s, vlen), otiny, vlen);

  VFLOAT y_base_hi = clln_hi;
  VFLOAT y_base_lo = clln_lo;
  VFLOAT clln_sel_hi = VMERGE_VVM_F (clld_hi, clln_hi, oref, vlen);
  VFLOAT clln_sel_lo = VMERGE_VVM_F (clld_lo, clln_lo, oref, vlen);
  clln_hi = VMERGE_VVM_F
    (clln_sel_hi, VFMV_VF (1.0f, vlen), otiny, vlen);
  clln_lo = VMERGE_VVM_F
    (clln_sel_lo, VFMV_VF (0.0f, vlen), otiny, vlen);

  VFLOAT x_ref_hi = VFMV_VF (0.0f, vlen);
  VFLOAT x_ref_lo = VFMV_VF (0.0f, vlen);
  if (__riscv_vcpop (oref, vlen) > 0)
    {
      VFLOAT a_q = VFMUL_VF (a, RVV_TGAMMAF_SIN_SCALE_F, vlen);
      VINT q = RVV_SLEEF_ORIG_TGAMMAF_MANGLE
        (rvv_sleef_orig_tgammaf_truncate_m) (a_q, vlen);
      VFLOAT q_f = __riscv_vfcvt_f (q, vlen);
      VFLOAT t_rem = VFSUB_VV
        (a, VFMUL_VF (q_f, RVV_TGAMMAF_SIN_UNSCALE_F, vlen), vlen);
      VFLOAT sin_hi, sin_lo;
      RVV_SLEEF_ORIG_TGAMMAF_MANGLE
        (rvv_sleef_orig_tgammaf_sinpifk_m)
        (t_rem, &sin_hi, &sin_lo, vlen);
      RVV_SLEEF_ORIG_DFMUL_VF2_VF2
        (clld_hi, clld_lo, sin_hi, sin_lo, x_ref_hi, x_ref_lo, vlen);
    }

  VFLOAT clld_sel_hi = VMERGE_VVM_F
    (y_base_hi, x_ref_hi, oref, vlen);
  VFLOAT clld_sel_lo = VMERGE_VVM_F
    (y_base_lo, x_ref_lo, oref, vlen);
  VFLOAT a_scaled = VFMUL_VF (a, RVV_TGAMMAF_TINY_SCALE_F, vlen);
  clld_hi = VMERGE_VVM_F
    (clld_sel_hi, a_scaled, otiny, vlen);
  clld_lo = VMERGE_VVM_F
    (clld_sel_lo, VFMV_VF (0.0f, vlen), otiny, vlen);

  *out_clc_hi = clc_hi;
  *out_clc_lo = clc_lo;
  *out_clln_hi = clln_hi;
  *out_clln_lo = clln_lo;
  *out_clld_hi = clld_hi;
  *out_clld_lo = clld_lo;
}

#  undef RVV_SLEEF_ORIG_TGAMMAF_HELPERS_THIS_LMUL
# endif
#endif

/* Target-owned source-shaped log2f exponent helpers.  The clean log2f path,
   the log2f comparator, and the log2p1f component composites share these
   exact bit-level adapters.  This part of v_math.h is intentionally
   re-entered for each importer LMUL, so guard the TU-level helper definitions
   explicitly and emit all four LMUL spellings only once.  */
#ifndef RVV_LOG2F_COMMON_HELPERS
#define RVV_LOG2F_COMMON_HELPERS
#define RVV_LOG2F_VTYPE(kind, lmul) \
  __PASTE3 (kind, __PASTE3 (32, m, lmul), _t)
#define RVV_LOG2F_REINTERPRET(a, b, lmul) \
  __PASTE5 (__riscv_vreinterpret_v_, a, \
            __PASTE4 (32, m, lmul, _), b, __PASTE3 (32, m, lmul))
#define RVV_LOG2F_DEFINE_HELPERS(lmul)                                      \
static inline RVV_LOG2F_VTYPE (vint, lmul)                                  \
__PASTE2 (rvv_sleef_orig_log2f_vilogb2k_m, lmul)                             \
  (RVV_LOG2F_VTYPE (vfloat, lmul) d, size_t vlen)                            \
{                                                                            \
  /* SLEEF vilogb2k extracts the unbiased exponent of a normalized float.  */ \
  RVV_LOG2F_VTYPE (vuint, lmul) q =                                        \
    __riscv_vsrl (RVV_LOG2F_REINTERPRET (f, u, lmul) (d), 23, vlen);         \
  q = __riscv_vand (q, 0xff, vlen);                                         \
  return __riscv_vsub (RVV_LOG2F_REINTERPRET (u, i, lmul) (q), 127, vlen);    \
}                                                                            \
                                                                             \
static inline RVV_LOG2F_VTYPE (vfloat, lmul)                                 \
__PASTE2 (rvv_sleef_orig_log2f_vldexp3_m, lmul)                              \
  (RVV_LOG2F_VTYPE (vfloat, lmul) d, RVV_LOG2F_VTYPE (vint, lmul) q,         \
   size_t vlen)                                                              \
{                                                                            \
  /* SLEEF vldexp3 adds q to the exponent field of a normalized float.  */    \
  return RVV_LOG2F_REINTERPRET (i, f, lmul)                                 \
    (__riscv_vadd (RVV_LOG2F_REINTERPRET (f, i, lmul) (d),                   \
                   __riscv_vsll (q, 23, vlen), vlen));                      \
}

#if (BIT_WIDTH == 32) && RVV_MATH_IS_LOG2F_FAMILY (V_MATH_FUNC)
RVV_LOG2F_DEFINE_HELPERS (1)
RVV_LOG2F_DEFINE_HELPERS (2)
RVV_LOG2F_DEFINE_HELPERS (4)
RVV_LOG2F_DEFINE_HELPERS (8)
#endif

#undef RVV_LOG2F_DEFINE_HELPERS
#undef RVV_LOG2F_REINTERPRET
#undef RVV_LOG2F_VTYPE
#endif /* RVV_LOG2F_COMMON_HELPERS */

/* Shared asinhf constants and importer mangles.  Keep the clean and
   source-original operation graphs independently named while placing their
   public definitions in the common math header.  */
#ifndef RVV_ASINHF_COMMON_DEFINITIONS
# define RVV_ASINHF_COMMON_DEFINITIONS
# define RVV_ASINHF_SQRT_FLT_MAX 1.8446743e+19f
# define RVV_ASINHF_L2U 0.69314718246459960938f
# define RVV_ASINHF_L2L -1.904654323148236017e-09f
# define RVV_ASINHF_LOG_P0 RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_ASINHF_LOG_P1 RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_ASINHF_LOG_P2 RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_ASINHF_LOG_P3 RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_ASINHF_CLEAN_LOGK2F \
  __PASTE2 (rvv_libmvec_asinhf_logk2f_m, LMUL)
# define RVV_ASINHF_CLEAN_LOGF_LARGE \
  __PASTE2 (rvv_libmvec_asinhf_logf_large_m, LMUL)
# define RVV_SLEEF_ORIG_ASINHF_SQRT_FLT_MAX 1.8446743e+19f
# define RVV_SLEEF_ORIG_ASINHF_L2U 0.69314718246459960938f
# define RVV_SLEEF_ORIG_ASINHF_L2L -1.904654323148236017e-09f
# define RVV_SLEEF_ORIG_ASINHF_LOG_P0 RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P1 RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P2 RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_SLEEF_ORIG_ASINHF_LOG_P3 RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_SLEEF_ORIG_ASINHF_MANGLE(base, lmul) \
  RVV_SLEEF_ORIG_MANGLE (base, lmul)
#endif

/* Shared importer mangles for the single-precision atanh source-shaped
   helpers.  The operation bodies remain in the routine skeletons, while
   these LMUL-dependent names live in the common math header so the skeletons
   contain no routine-local helper macro block.  */
#ifndef RVV_ATANHF_COMMON_DEFINITIONS
# define RVV_ATANHF_COMMON_DEFINITIONS
# define RVV_ATANHF_LOG_POLY0_F RVV_SLEEF_LOGK2F_POLY0_F
# define RVV_ATANHF_LOG_POLY1_F RVV_SLEEF_LOGK2F_POLY1_F
# define RVV_ATANHF_LOG_POLY2_F RVV_SLEEF_LOGK2F_POLY2_F
# define RVV_ATANHF_LOG_POLY3_F RVV_SLEEF_LOGK2F_POLY3_F
# define RVV_LIBMVEC_ATANHF_VILOGBK \
  __PASTE2 (rvv_libmvec_atanhf_vilogbk_m, LMUL)
# define RVV_LIBMVEC_ATANHF_VPOW2I \
  __PASTE2 (rvv_libmvec_atanhf_vpow2i_m, LMUL)
# define RVV_LIBMVEC_ATANHF_LOGK2F \
  __PASTE2 (rvv_libmvec_atanhf_logk2f_m, LMUL)
#endif

/* The double-precision atan2 and atan2pi skeletons share one operation graph.
   Keep their mode-specific constants and fixed-point post-selection in this
   common header; the skeleton itself should not carry a compile-control
   define merely because the importer has a different public function name.  */
#ifndef RVV_MATH_IS_ATAN2_FAMILY_atan2
# define RVV_MATH_IS_ATAN2_FAMILY_atan2 1
# define RVV_MATH_IS_ATAN2_FAMILY_atan2pi 1
# define RVV_MATH_IS_ATAN2_FAMILY_INNER(name) \
  RVV_MATH_IS_ATAN2_FAMILY_##name
# define RVV_MATH_IS_ATAN2_FAMILY(name) \
  RVV_MATH_IS_ATAN2_FAMILY_INNER (name)
#endif

#ifndef RVV_MATH_IS_ATAN2PI_atan2pi
# define RVV_MATH_IS_ATAN2PI_atan2pi 1
# define RVV_MATH_IS_ATAN2PI_INNER(name) RVV_MATH_IS_ATAN2PI_##name
# define RVV_MATH_IS_ATAN2PI(name) RVV_MATH_IS_ATAN2PI_INNER (name)
#endif

#undef RVV_ATAN2_Y_NONZERO_X_ZERO
#undef RVV_ATAN2_Y_ZERO_X_NEGATIVE
#undef RVV_ATAN2_POSINF_POSINF
#undef RVV_ATAN2_NEGINF_INF
#undef RVV_ATAN2_NEGINF_FINITE
#undef RVV_ATAN2_FINITE_NEGATIVE_ZERO
#undef RVV_ATAN2_Y_INF_FINITE
#undef RVV_ATAN2_PI_HI
#undef RVV_ATAN2_PI_LO
#undef RVV_ATAN2_PIBY2_HI
#undef RVV_ATAN2_FIXED_HALF
#undef RVV_ATAN2_FIXED_ONE
#undef RVV_ATAN2_POST_SCALE

#if (BIT_WIDTH == 64) && RVV_MATH_IS_ATAN2_FAMILY (V_MATH_FUNC)
# if RVV_MATH_IS_ATAN2PI (V_MATH_FUNC)
#  define RVV_ATAN2_Y_NONZERO_X_ZERO (fp_posHalf)
#  define RVV_ATAN2_Y_ZERO_X_NEGATIVE (fp_posOne)
#  define RVV_ATAN2_POSINF_POSINF RVV_ATAN2_ONE_OVER_4_D
#  define RVV_ATAN2_NEGINF_INF RVV_ATAN2_THREE_OVER_4_D
#  define RVV_ATAN2_NEGINF_FINITE (fp_posOne)
#  define RVV_ATAN2_FINITE_NEGATIVE_ZERO (fp_posOne)
#  define RVV_ATAN2_Y_INF_FINITE (fp_posHalf)
#  define RVV_ATAN2_FIXED_HALF RVV_ATAN2_HALF_Q61
#  define RVV_ATAN2_FIXED_ONE RVV_ATAN2_ONE_Q61
#  define RVV_ATAN2_POST_SCALE(value, vlen) \
  __riscv_vsmul ((value), RVV_ONE_OVER_PI_Q63, 1, (vlen))
# else
#  define RVV_ATAN2_Y_NONZERO_X_ZERO PIBY2_HI
#  define RVV_ATAN2_Y_ZERO_X_NEGATIVE RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_POSINF_POSINF RVV_ATAN2_PIBY4_HI_D
#  define RVV_ATAN2_NEGINF_INF RVV_ATAN2_PI_075_HI_D
#  define RVV_ATAN2_NEGINF_FINITE RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_FINITE_NEGATIVE_ZERO RVV_ATAN2_PI_HI_D
#  define RVV_ATAN2_Y_INF_FINITE PIBY2_HI
#  define RVV_ATAN2_FIXED_HALF RVV_PIBY2_Q61
#  define RVV_ATAN2_FIXED_ONE RVV_ATAN2_PI_Q61
#  define RVV_ATAN2_POST_SCALE(value, vlen) (value)
# endif
# define RVV_ATAN2_PI_HI RVV_ATAN2_PI_HI_D
# define RVV_ATAN2_PI_LO RVV_ATAN2_PI_LO_D
# define RVV_ATAN2_PIBY2_HI PIBY2_HI
#endif

/* Clean asinhf is the first importer that needs this subset of the generic
   float-float operation family.  Keep only the operations used by asinhf
   here; later routines can extend the common family at their first use.  */
#ifndef RVV_ASINHF_CLEAN_DF_HELPERS
# define RVV_ASINHF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF_S(a_hi, a_lo, b, r_hi, r_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, VFMV_VF ((b), (vlen)), (vlen));            \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),      \
                     VFRSUB_VF (_v, (b), (vlen)), (vlen)),                    \
          _a_lo, (vlen));                                                     \
    }                                                                         \
  while (0)

# define DFADD_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _b_hi, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _b_hi, (vlen)),              \
          VFADD_VV (_a_lo, _b_lo, (vlen)), (vlen));                            \
    }                                                                         \
  while (0)

# define DFADD_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a_hi, _b, (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _b, (vlen)), _a_lo,         \
          (vlen));                                                            \
    }                                                                         \
  while (0)

# define DFSQU_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      (r_hi) = VFMUL_VV (_a_hi, _a_hi, (vlen));                               \
      (r_lo) = VFMSUB_VV (_a_hi, _a_hi, (r_hi), (vlen));                      \
      (r_lo) = VFMACC_VV ((r_lo), VFADD_VV (_a_hi, _a_hi, (vlen)), _a_lo,    \
                          (vlen));                                           \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b, _s, (vlen));                          \
      _t = VFMACC_VV (_t, _a_lo, _b, (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFMUL_VV (_a_hi, _b_hi, (vlen));                            \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b_hi, _s, (vlen));                       \
      _t = VFMACC_VV (_t, _a_hi, _b_lo, (vlen));                              \
      _t = VFMACC_VV (_t, _a_lo, _b_hi, (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFDIV_VF2_VF2_VF2(n_hi, n_lo, d_hi, d_lo, r_hi, r_lo, vlen)          \
  do                                                                          \
    {                                                                         \
      VFLOAT _n_hi = (n_hi);                                                  \
      VFLOAT _n_lo = (n_lo);                                                  \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _s = VFMUL_VV (_n_hi, _t, (vlen));                               \
      VFLOAT _u = VFNMSUB_VV (_s, _d_hi, _n_hi, (vlen));                      \
      _u = VFNMSAC_VV (_u, _s, _d_lo, (vlen));                                \
      _u = VFADD_VV (_u, _n_lo, (vlen));                                      \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMUL_VV (_u, _t, (vlen));                                     \
    }                                                                         \
  while (0)

# define DFSCALE_VF2_VF2_S(a_hi, a_lo, s, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      (r_hi) = VFMUL_VF ((a_hi), (s), (vlen));                                \
      (r_lo) = VFMUL_VF ((a_lo), (s), (vlen));                                \
    }                                                                         \
  while (0)

# define DFSCALE_VF2_VF2_VF(a_hi, a_lo, s, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _s = (s);                                                        \
      (r_hi) = VFMUL_VV ((a_hi), _s, (vlen));                                 \
      (r_lo) = VFMUL_VV ((a_lo), _s, (vlen));                                 \
    }                                                                         \
  while (0)

#endif /* RVV_ASINHF_CLEAN_DF_HELPERS */

/* The atanf family uses these two generic float-float operations from the
   common helper family.  Keep them in v_math.h so the routine skeletons do
   not need a private helper header; the remaining generic operations can be
   added here by their first real user.  */
#ifndef RVV_VFLOAT2_COMMON_HELPERS
# define RVV_VFLOAT2_COMMON_HELPERS

# define VSEL_VF2_VO(o, a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)             \
  do                                                                          \
    {                                                                         \
      (r_hi) = __riscv_vmerge ((b_hi), (a_hi), (o), (vlen));                   \
      (r_lo) = __riscv_vmerge ((b_lo), (a_lo), (o), (vlen));                   \
    }                                                                         \
  while (0)

# define DFNEG_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                           \
  do                                                                          \
    {                                                                         \
      (r_hi) = VFNEG_V ((a_hi), (vlen));                                       \
      (r_lo) = VFNEG_V ((a_lo), (vlen));                                       \
    }                                                                         \
  while (0)

#endif /* RVV_VFLOAT2_COMMON_HELPERS */

/* Float-float operations used by the clean cosf kernel.  Keep these
   importer-dependent operations in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSF_CLEAN_DF_HELPERS
# define RVV_COSF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF_VF(a, b, r_hi, r_lo, vlen)                            \
  do                                                                          \
    {                                                                         \
      VFLOAT _a = (a);                                                        \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a, _b, (vlen));                                 \
      VFLOAT _v = VFSUB_VV (_s, _a, (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                      \
          VFSUB_VV (_a, VFSUB_VV (_s, _v, (vlen)), (vlen)),                  \
          VFSUB_VV (_b, _v, (vlen)), (vlen));                                \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a_hi, _b, (vlen));                              \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                              \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_b, _v, (vlen)), (vlen));                                 \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _a_lo, (vlen));                                 \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF2_VF2(a_hi, a_lo, b_hi, b_lo, r_hi, r_lo, vlen)       \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b_hi = (b_hi);                                                  \
      VFLOAT _b_lo = (b_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _b_hi, (vlen));                           \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                              \
      VFLOAT _t = VFADD_VV (                                                   \
          VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),                \
          VFSUB_VV (_b_hi, _v, (vlen)), (vlen));                              \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, VFADD_VV (_a_lo, _b_lo, (vlen)), (vlen));       \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF_VF(a, b, r_hi, r_lo, vlen)                              \
  do                                                                          \
    {                                                                         \
      VFLOAT _a = (a);                                                        \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a, _b, (vlen));                                  \
      (r_hi) = _s;                                                            \
      (r_lo) = VFMSUB_VV (_a, _b, _s, (vlen));                                \
    }                                                                         \
  while (0)

# define DFMUL_VF2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFMUL_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _t = VFMSUB_VV (_a_hi, _b, _s, (vlen));                          \
      _t = VFMACC_VV (_t, _a_lo, _b, (vlen));                                \
      (r_hi) = _s;                                                            \
      (r_lo) = _t;                                                            \
    }                                                                         \
  while (0)

# define DFMUL_VF_VF2_VF2(x_hi, x_lo, y_hi, y_lo, r, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _y_hi = (y_hi);                                                 \
      VFLOAT _y_lo = (y_lo);                                                 \
      VFLOAT _t = VFMUL_VV (_x_hi, _y_lo, (vlen));                           \
      _t = VFMACC_VV (_t, _x_lo, _y_hi, (vlen));                             \
      (r) = VFMACC_VV (_t, _x_hi, _y_hi, (vlen));                            \
    }                                                                         \
  while (0)

# define DFNORMALIZE_VF2_VF2(a_hi, a_lo, r_hi, r_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _s = VFADD_VV (_a_hi, _a_lo, (vlen));                            \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (VFSUB_VV (_a_hi, _s, (vlen)), _a_lo, (vlen));       \
    }                                                                         \
  while (0)

#endif /* RVV_COSF_CLEAN_DF_HELPERS */

/* Float-float reciprocal used by the clean coshf kernel.  Keep this
   importer-dependent operation in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSHF_CLEAN_DF_HELPERS
# define RVV_COSHF_CLEAN_DF_HELPERS

# define DFREC_VF2_VF2(d_hi, d_lo, r_hi, r_lo, vlen)                         \
  do                                                                          \
    {                                                                         \
      VFLOAT _d_hi = (d_hi);                                                  \
      VFLOAT _d_lo = (d_lo);                                                  \
      VFLOAT _t = VFDIV_VV (VFMV_VF (1.0f, (vlen)), _d_hi, (vlen));           \
      VFLOAT _u = VFNMSUB_VV (_t, _d_hi, VFMV_VF (1.0f, (vlen)), (vlen));      \
      _u = VFNMSAC_VV (_u, _t, _d_lo, (vlen));                                \
      (r_hi) = _t;                                                            \
      (r_lo) = VFMUL_VV (_u, _t, (vlen));                                     \
    }                                                                         \
  while (0)

#endif /* RVV_COSHF_CLEAN_DF_HELPERS */

/* Float-float additions used by the clean cospif kernel.  Keep these
   importer-dependent operations in the common math header; the skeleton
   itself should expose only the algorithm and its data accesses.  */
#ifndef RVV_COSPIF_CLEAN_DF_HELPERS
# define RVV_COSPIF_CLEAN_DF_HELPERS

# define DFADD2_VF2_VF(a_hi, a_lo, b, r_hi, r_lo, vlen)                      \
  do                                                                          \
    {                                                                         \
      VFLOAT _a_hi = (a_hi);                                                  \
      VFLOAT _a_lo = (a_lo);                                                  \
      VFLOAT _b = (b);                                                        \
      VFLOAT _s = VFADD_VV (_a_hi, _b, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _a_hi, (vlen));                               \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (                                                     \
          VFADD_VV (VFSUB_VV (_a_hi, VFSUB_VV (_s, _v, (vlen)), (vlen)),      \
                    VFSUB_VV (_b, _v, (vlen)), (vlen)),                       \
          _a_lo, (vlen));                                                     \
    }                                                                         \
  while (0)

# define DFADD2_VF2_VF_VF2(x, y_hi, y_lo, r_hi, r_lo, vlen)                   \
  do                                                                          \
    {                                                                         \
      VFLOAT _x = (x);                                                        \
      VFLOAT _y_hi = (y_hi);                                                  \
      VFLOAT _y_lo = (y_lo);                                                  \
      VFLOAT _s = VFADD_VV (_x, _y_hi, (vlen));                               \
      VFLOAT _v = VFSUB_VV (_s, _x, (vlen));                                  \
      VFLOAT _t = VFADD_VV (VFSUB_VV (_x, VFSUB_VV (_s, _v, (vlen)), (vlen)), \
                            VFSUB_VV (_y_hi, _v, (vlen)), (vlen));             \
      (r_hi) = _s;                                                            \
      (r_lo) = VFADD_VV (_t, _y_lo, (vlen));                                  \
    }                                                                         \
  while (0)

#endif /* RVV_COSPIF_CLEAN_DF_HELPERS */

/* Fixed scalar coefficients and endpoint constants for the single-precision
   atan family.  */
#ifndef RVV_ATANF_CONSTANTS
# define RVV_ATANF_CONSTANTS
/* Polynomial shared by the single-precision atan kernels.  */
# define RVV_ATANF_POLY0_F -0.00176397908944636583328247f
# define RVV_ATANF_POLY1_F 0.0107900900766253471374512f
# define RVV_ATANF_POLY2_F -0.0309564601629972457885742f
# define RVV_ATANF_POLY3_F 0.0577365085482597351074219f
# define RVV_ATANF_POLY4_F -0.0838950723409652709960938f
# define RVV_ATANF_POLY5_F 0.109463557600975036621094f
# define RVV_ATANF_POLY6_F -0.142626821994781494140625f
# define RVV_ATANF_POLY7_F 0.199983194470405578613281f
# define RVV_ATANF_POLY8_F -0.333332866430282592773438f

/* Input bit patterns for atanpif endpoint repairs.  */
# define RVV_ATANPIF_FIX_NEG_0_U 0xbf7fb455u
# define RVV_ATANPIF_FIX_NEG_1_U 0xbf7e6ad9u
# define RVV_ATANPIF_FIX_NEG_2_U 0xbf78335fu
# define RVV_ATANPIF_FIX_NEG_3_U 0xbf7670b7u
# define RVV_ATANPIF_FIX_NEG_4_U 0xbf762a13u
# define RVV_ATANPIF_FIX_NEG_5_U 0xbf71578au
# define RVV_ATANPIF_FIX_NEG_6_U 0xbf70b248u
# define RVV_ATANPIF_FIX_NEG_7_U 0xbf6eda9du
# define RVV_ATANPIF_FIX_NEG_8_U 0xbf6e9333u
# define RVV_ATANPIF_FIX_NEG_9_U 0xbf6941b1u
# define RVV_ATANPIF_FIX_POS_0_U 0x3f77419du
# define RVV_ATANPIF_FIX_POS_1_U 0x3f70f56fu
# define RVV_ATANPIF_FIX_POS_2_U 0x3f72d2e7u
# define RVV_ATANPIF_FIX_POS_3_U 0x3f7349e2u
# define RVV_ATANPIF_FIX_POS_4_U 0x3f77c250u
# define RVV_ATANPIF_FIX_POS_5_U 0x3f7bf501u
# define RVV_ATANPIF_FIX_POS_6_U 0x3f7eead8u

# define RVV_ATANF_PI_F RVV_SLEEF_TRIG_PI_F
# define RVV_ATANF_PI2_F 0x1.921fb6p+0f
# define RVV_ATANF_PI4_F 0x1.921fb6p-1f
# define RVV_ATANF_PI2_HI_F 1.5707963705062866211f
# define RVV_ATANF_PI2_LO_F -4.3711388286737928865e-08f
# define RVV_ATANF_INV_PI_HI_F RVV_SLEEF_TRIG_ONE_OVER_PI
# define RVV_ATANF_INV_PI_LO_F 0x1.b93910p-27f
# define RVV_ATANF_TINY_F 2.9387372783541830947e-39f
# define RVV_ATANF_SCALE_F 16777216.0f
# define RVV_ATANF_NORM_PI_F 1.0f
# define RVV_ATANF_NORM_PI2_F 0.5f
# define RVV_ATANF_NORM_PI4_F 0.25f
# define RVV_ATANPIF_BIAS_F 0x1.0p-25f
#endif

/* LMUL-qualified helpers used by the clean single-precision atan family.
   Importer wrappers may explicitly activate the block, while direct
   skeleton users select it from V_MATH_FUNC.  Per-LMUL guards are required
   because v_math.h is re-included for every importer instantiation.  */
#ifndef RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atanf
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atanf 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atan2f 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_atan2pif 1
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_INNER(name)                  \
  RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_##name
# define RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG(name)                        \
  RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG_INNER (name)
#endif

/* Masked double-float helpers for the source-shaped log1pf finite graph.
   These are used only by clean log1pf when finite and EH-special lanes are
   mixed.  Keep the active-lane operation order identical to the unmasked
   source helpers above.  */
#ifndef RVV_SLEEF_LOG1PF_MASKED_DF
# define RVV_SLEEF_LOG1PF_MASKED_DF
# define RVV_SLEEF_LOG1PF_ADD_VV_M(mask, a, b, vlen) \
  __riscv_vfadd ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_ADD_VF_M(mask, a, b, vlen) \
  __riscv_vfadd ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_SUB_VV_M(mask, a, b, vlen) \
  __riscv_vfsub ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_SUB_VF_M(mask, a, b, vlen) \
  __riscv_vfsub ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_MUL_VV_M(mask, a, b, vlen) \
  __riscv_vfmul ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_MUL_VF_M(mask, a, b, vlen) \
  __riscv_vfmul ((mask), (a), (b), (vlen))
# define RVV_SLEEF_LOG1PF_MADD_VV_M(mask, a, b, c, vlen) \
  __riscv_vfmadd ((mask), (a), (b), (c), (vlen))
# define RVV_SLEEF_LOG1PF_MSUB_VV_M(mask, a, b, c, vlen) \
  __riscv_vfmsub ((mask), (a), (b), (c), (vlen))
# define RVV_SLEEF_LOG1PF_NMSUB_VV_M(mask, a, b, c, vlen) \
  __riscv_vfnmsub ((mask), (a), (b), (c), (vlen))
# define RVV_SLEEF_LOG1PF_DIV_VV_M(mask, a, b, vlen) \
  __riscv_vfdiv ((mask), (a), (b), (vlen))

# define RVV_SLEEF_LOG1PF_DFMUL_VF2_VF_M(mask, x_hi, x_lo, y, r_hi, r_lo, \
                                         vlen)                                \
  do                                                                         \
    {                                                                        \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _y = (y);                                                       \
      VFLOAT _s = RVV_SLEEF_LOG1PF_MUL_VV_M ((mask), _x_hi, _y, (vlen));     \
      (r_hi) = _s;                                                           \
      (r_lo) = RVV_SLEEF_LOG1PF_MADD_VV_M                                  \
        ((mask), _x_lo, _y,                                                   \
         RVV_SLEEF_LOG1PF_MSUB_VV_M ((mask), _x_hi, _y, _s, (vlen)),         \
         (vlen));                                                            \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFADD2_VF_VF_M(mask, x, y, r_hi, r_lo, vlen)       \
  do                                                                         \
    {                                                                        \
      VFLOAT _x = (x);                                                       \
      VFLOAT _y = (y);                                                       \
      VFLOAT _s = RVV_SLEEF_LOG1PF_ADD_VV_M ((mask), _x, _y, (vlen));       \
      VFLOAT _v = RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _s, _x, (vlen));       \
      (r_hi) = _s;                                                           \
      (r_lo) = RVV_SLEEF_LOG1PF_ADD_VV_M                                   \
        ((mask),                                                              \
         RVV_SLEEF_LOG1PF_SUB_VV_M                                          \
           ((mask), _x,                                                      \
            RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _s, _v, (vlen)), (vlen)),      \
         RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _y, _v, (vlen)), (vlen));        \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFADD_VF_VF_M(mask, x, y, r_hi, r_lo, vlen)        \
  do                                                                         \
    {                                                                        \
      VFLOAT _x = (x);                                                       \
      VFLOAT _y = (y);                                                       \
      VFLOAT _s = RVV_SLEEF_LOG1PF_ADD_VV_M ((mask), _x, _y, (vlen));       \
      (r_hi) = _s;                                                           \
      (r_lo) = RVV_SLEEF_LOG1PF_ADD_VV_M                                   \
        ((mask), RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _x, _s, (vlen)), _y,    \
         (vlen));                                                            \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_M(mask, x_hi, x_lo, y_hi, y_lo,     \
                                          r_hi, r_lo, vlen)                  \
  do                                                                         \
    {                                                                        \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _y_hi = (y_hi);                                                 \
      VFLOAT _y_lo = (y_lo);                                                 \
      VFLOAT _s = RVV_SLEEF_LOG1PF_ADD_VV_M                                \
        ((mask), _x_hi, _y_hi, (vlen));                                     \
      VFLOAT _t = RVV_SLEEF_LOG1PF_ADD_VV_M                                \
        ((mask), RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _x_hi, _s, (vlen)),      \
         _y_hi, (vlen));                                                     \
      _t = RVV_SLEEF_LOG1PF_ADD_VV_M ((mask), _t, _x_lo, (vlen));            \
      _t = RVV_SLEEF_LOG1PF_ADD_VV_M ((mask), _t, _y_lo, (vlen));            \
      (r_hi) = _s;                                                           \
      (r_lo) = _t;                                                           \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFADD_VF2_VF2_VF_M(mask, x_hi, x_lo, y, r_hi,     \
                                             r_lo, vlen)                     \
  do                                                                         \
    {                                                                        \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _y = (y);                                                       \
      VFLOAT _s = RVV_SLEEF_LOG1PF_ADD_VV_M                                \
        ((mask), _x_hi, _y, (vlen));                                         \
      VFLOAT _t = RVV_SLEEF_LOG1PF_ADD_VV_M                                \
        ((mask), RVV_SLEEF_LOG1PF_SUB_VV_M ((mask), _x_hi, _s, (vlen)),      \
         _y, (vlen));                                                        \
      (r_hi) = _s;                                                           \
      (r_lo) = RVV_SLEEF_LOG1PF_ADD_VV_M ((mask), _t, _x_lo, (vlen));        \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFSCALE_VF2_VF_M(mask, x_hi, x_lo, scale, r_hi,   \
                                           r_lo, vlen)                       \
  do                                                                         \
    {                                                                        \
      VFLOAT _x_hi = (x_hi);                                                 \
      VFLOAT _x_lo = (x_lo);                                                 \
      VFLOAT _scale = (scale);                                               \
      (r_hi) = RVV_SLEEF_LOG1PF_MUL_VV_M                                    \
        ((mask), _x_hi, _scale, (vlen));                                     \
      (r_lo) = RVV_SLEEF_LOG1PF_MUL_VV_M                                    \
        ((mask), _x_lo, _scale, (vlen));                                     \
    }                                                                        \
  while (0)

# define RVV_SLEEF_LOG1PF_DFDIV_VF2_VF2_M(mask, n_hi, n_lo, d_hi, d_lo,     \
                                          r_hi, r_lo, vlen)                  \
  do                                                                         \
    {                                                                        \
      VFLOAT _n_hi = (n_hi);                                                 \
      VFLOAT _n_lo = (n_lo);                                                 \
      VFLOAT _d_hi = (d_hi);                                                 \
      VFLOAT _d_lo = (d_lo);                                                 \
      VFLOAT _one = VFMV_VF (1.0f, (vlen));                                  \
      VFLOAT _t = RVV_SLEEF_LOG1PF_DIV_VV_M                                  \
        ((mask), _one, _d_hi, (vlen));                                       \
      VFLOAT _s = RVV_SLEEF_LOG1PF_MUL_VV_M                                 \
        ((mask), _n_hi, _t, (vlen));                                         \
      VFLOAT _u = RVV_SLEEF_LOG1PF_MSUB_VV_M                                \
        ((mask), _t, _n_hi, _s, (vlen));                                     \
      VFLOAT _v = RVV_SLEEF_LOG1PF_NMSUB_VV_M                               \
        ((mask), _d_hi, _t, _one, (vlen));                                   \
      _v = RVV_SLEEF_LOG1PF_NMSUB_VV_M                                      \
        ((mask), _d_lo, _t, _v, (vlen));                                     \
      VFLOAT _nlo_t = RVV_SLEEF_LOG1PF_MADD_VV_M                             \
        ((mask), _n_lo, _t, _u, (vlen));                                     \
      (r_hi) = _s;                                                           \
      (r_lo) = RVV_SLEEF_LOG1PF_MADD_VV_M                                    \
        ((mask), _s, _v, _nlo_t, (vlen));                                    \
    }                                                                        \
  while (0)
#endif

/* Shared coshf importer names.  The clean source-shaped helper and the
   source-original comparator use distinct function prefixes, but both need
   the same per-LMUL naming discipline.  Keep these names here so neither
   coshf skeleton carries a routine-local macro block.  */
#undef RVV_COSHF_SOURCE_MANGLE
#define RVV_COSHF_SOURCE_MANGLE(base) __PASTE3 (base, _m, LMUL)
#undef RVV_COSHF_EXPK2F_VF2
#define RVV_COSHF_EXPK2F_VF2 \
  __PASTE4 (coshf_expk2f_vf2_, LMUL, _, SIMDLEN)
#undef RVV_SLEEF_ORIG_COSHF_MANGLE
#define RVV_SLEEF_ORIG_COSHF_MANGLE(base) \
  RVV_SLEEF_ORIG_MANGLE (base, LMUL)
#undef RVV_SLEEF_ORIG_COSHF_VPOW2I
#define RVV_SLEEF_ORIG_COSHF_VPOW2I \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_vpow2i_m)
#undef RVV_SLEEF_ORIG_COSHF_VLDEXP2
#define RVV_SLEEF_ORIG_COSHF_VLDEXP2 \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_vldexp2_m)
#undef RVV_SLEEF_ORIG_COSHF_EXPK2F
#define RVV_SLEEF_ORIG_COSHF_EXPK2F \
  RVV_SLEEF_ORIG_COSHF_MANGLE (rvv_sleef_orig_coshf_expk2f_m)

/* Shared source-association helper names for the single-precision cospif
   finite graph.  The helper body remains in the clean skeleton, while its
   importer-qualified name is owned here.  */
#undef RVV_COSPIF_SOURCE_MANGLE
#define RVV_COSPIF_SOURCE_MANGLE(base) __PASTE3 (base, _m, LMUL)
#undef RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR
#define RVV_COSPIF_SOURCE_DFMUL_PAIR_PAIR \
  RVV_COSPIF_SOURCE_MANGLE (cospif_source_dfmul_pair_pair)

#if defined(RVV_LIBMVEC_ATAN2F_HELPERS)                                        \
  || RVV_LIBMVEC_CLEAN_ATAN_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M1
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M1
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M2
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M2
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M4
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M4
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_LIBMVEC_ATAN2F_HELPERS_M8
#   define RVV_LIBMVEC_ATAN2F_HELPERS_M8
#   define RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
static inline VINT
RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
  (VFLOAT f0, VFLOAT f1, VINT x, VINT y, size_t vlen)
{
  VBOOL m = __riscv_vmflt (f0, f1, vlen);
  return __riscv_vmerge (y, x, m, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (visinf2_f, LMUL)
  (VFLOAT d, VFLOAT m, size_t vlen)
{
  VUINT abs_d = __riscv_vand (F_AS_U (d), 0x7fffffff, vlen);
  VBOOL is_inf = __riscv_vmseq (abs_d, 0x7f800000, vlen);
  VUINT sign = __riscv_vand (F_AS_U (d), 0x80000000, vlen);
  VUINT out_bits = __riscv_vor (F_AS_U (m), sign, vlen);
  VFLOAT out = U_AS_F (out_bits);
  return __riscv_vmerge (VFMV_VF (0.0f, vlen), out, is_inf, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (atan2kf_u1, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT x_hi, VFLOAT x_lo, VFLOAT *r_hi,
   VFLOAT *r_lo, size_t vlen)
{
  VINT q = RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
    (x_hi, VFMV_VF (0.0f, vlen), VMVI_VX (-2, vlen), VMVI_VX (0, vlen),
     vlen);
  VBOOL p = __riscv_vmflt (x_hi, 0.0f, vlen);
  VUINT sign = VMVU_VX (0x80000000, vlen);
  VUINT flip = __riscv_vmerge (VMVU_VX (0, vlen), sign, p, vlen);
  x_hi = U_AS_F (__riscv_vxor (F_AS_U (x_hi), flip, vlen));
  x_lo = U_AS_F (__riscv_vxor (F_AS_U (x_lo), flip, vlen));

  q = RVV_LIBMVEC_IMPL_NAME (vsel_vi2_vf_vf_vi2_vi2, LMUL)
    (x_hi, y_hi, __riscv_vadd (q, 1, vlen), q, vlen);
  p = __riscv_vmflt (x_hi, y_hi, vlen);

  VFLOAT neg_x_hi, neg_x_lo;
  DFNEG_VF2_VF2 (x_hi, x_lo, neg_x_hi, neg_x_lo, vlen);

  VFLOAT s_hi, s_lo, t_hi, t_lo;
  VSEL_VF2_VO (p, neg_x_hi, neg_x_lo, y_hi, y_lo, s_hi, s_lo, vlen);
  VSEL_VF2_VO (p, y_hi, y_lo, x_hi, x_lo, t_hi, t_lo, vlen);

  RVV_SLEEF_ORIG_DFDIV_VF2_VF2
    (s_hi, s_lo, t_hi, t_lo, s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFSQU_VF2
    (s_hi, s_lo, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2
    (t_hi, t_lo, t_hi, t_lo, vlen);

  VFLOAT u = VFMV_VF (RVV_ATANF_POLY0_F, vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY1_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY2_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY3_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY4_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY5_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY6_F, vlen), vlen);
  u = __riscv_vfmadd (u, t_hi, VFMV_VF (RVV_ATANF_POLY7_F, vlen), vlen);

  VFLOAT u_t = __riscv_vfmul (u, t_hi, vlen);
  VFLOAT poly_hi, poly_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_ATANF_POLY8_F, vlen), u_t, poly_hi, poly_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (t_hi, t_lo, poly_hi, poly_lo, t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2
    (VFMV_VF (1.0f, vlen), t_hi, t_lo, sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (s_hi, s_lo, sum_hi, sum_lo, t_hi, t_lo, vlen);

  VFLOAT qf = __riscv_vfcvt_f (q, vlen);
  VFLOAT add_hi, add_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ATANF_PI2_HI_F, vlen),
     VFMV_VF (RVV_ATANF_PI2_LO_F, vlen), qf,
     add_hi, add_lo, vlen);

  RVV_SLEEF_ORIG_DFADD_VF2_VF2
    (add_hi, add_lo, t_hi, t_lo, *r_hi, *r_lo, vlen);
}
#  undef RVV_LIBMVEC_ATAN2F_HELPERS_THIS_LMUL
# endif
#endif

/* Source-shaped SLEEF atan-family helpers.  They are kept separately
   mangled from the clean helper block so the comparison objects retain the
   source operation graph while sharing the public math header.  */
#ifndef RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanf_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atanpif_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atan2f_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_atan2pif_original 1
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_INNER(name)                  \
  RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_##name
# define RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG(name)                        \
  RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG_INNER (name)
#endif

#if defined(RVV_SLEEF_ATAN_HELPERS)                                           \
  || RVV_LIBMVEC_SLEEF_ATAN_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M1
#   define RVV_SLEEF_ATAN_HELPERS_M1
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M2
#   define RVV_SLEEF_ATAN_HELPERS_M2
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M4
#   define RVV_SLEEF_ATAN_HELPERS_M4
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_SLEEF_ATAN_HELPERS_M8
#   define RVV_SLEEF_ATAN_HELPERS_M8
#   define RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
  (VFLOAT f0, VFLOAT f1, VINT x, VINT y, size_t vlen)
{
  VBOOL m = __riscv_vmflt (f0, f1, vlen);
  return __riscv_vmerge (y, x, m, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_visinf2_m, LMUL)
  (VFLOAT d, VFLOAT m, size_t vlen)
{
  VUINT abs_d = __riscv_vand (F_AS_U (d), 0x7fffffff, vlen);
  VBOOL is_inf = __riscv_vmseq (abs_d, 0x7f800000, vlen);
  VUINT sign = __riscv_vand (F_AS_U (d), 0x80000000, vlen);
  VUINT out_bits = __riscv_vor (F_AS_U (m), sign, vlen);
  VUINT selected = __riscv_vmerge (VMVU_VX (0, vlen), out_bits, is_inf,
                                   vlen);
  return U_AS_F (selected);
}

static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan2kf_u1_m, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT x_hi, VFLOAT x_lo,
   VFLOAT *r_hi, VFLOAT *r_lo,
   size_t vlen)
{
  VINT q = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
    (x_hi, VFMV_VF (0.0f, vlen), VMVI_VX (-2, vlen), VMVI_VX (0, vlen),
     vlen);
  VBOOL p = __riscv_vmflt (x_hi, 0.0f, vlen);
  VUINT sign = VMVU_VX (0x80000000, vlen);
  VUINT flip = __riscv_vmerge (VMVU_VX (0, vlen), sign, p, vlen);
  x_hi = U_AS_F (__riscv_vxor (F_AS_U (x_hi), flip, vlen));
  x_lo = U_AS_F (__riscv_vxor (F_AS_U (x_lo), flip, vlen));

  q = RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan_select_m, LMUL)
    (x_hi, y_hi, __riscv_vadd (q, 1, vlen), q, vlen);
  p = __riscv_vmflt (x_hi, y_hi, vlen);

  VFLOAT neg_x_hi, neg_x_lo;
  RVV_SLEEF_ORIG_DFNEG_VF2 (x_hi, x_lo, neg_x_hi, neg_x_lo, vlen);

  VFLOAT s_hi = __riscv_vmerge (y_hi, neg_x_hi, p, vlen);
  VFLOAT s_lo = __riscv_vmerge (y_lo, neg_x_lo, p, vlen);
  VFLOAT t_hi = __riscv_vmerge (x_hi, y_hi, p, vlen);
  VFLOAT t_lo = __riscv_vmerge (x_lo, y_lo, p, vlen);

  RVV_SLEEF_ORIG_DFDIV_VF2_VF2 (s_hi, s_lo, t_hi, t_lo,
                                s_hi, s_lo, vlen);
  RVV_SLEEF_ORIG_DFSQU_VF2 (s_hi, s_lo, t_hi, t_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (t_hi, t_lo, t_hi, t_lo, vlen);

  VFLOAT u = VFMV_VF (RVV_ATANF_POLY0_F, vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY1_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY2_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY3_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY4_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY5_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY6_F, vlen), vlen);
  u = VFMADD_VV (u, t_hi, VFMV_VF (RVV_ATANF_POLY7_F, vlen), vlen);

  VFLOAT u_t = VFMUL_VV (u, t_hi, vlen);
  VFLOAT poly_hi, poly_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF
    (VFMV_VF (RVV_ATANF_POLY8_F, vlen), u_t, poly_hi, poly_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (t_hi, t_lo, poly_hi, poly_lo,
                                t_hi, t_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD_VF_VF2 (VFMV_VF (1.0f, vlen), t_hi, t_lo,
                               sum_hi, sum_lo, vlen);
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2 (s_hi, s_lo, sum_hi, sum_lo,
                                t_hi, t_lo, vlen);

  VFLOAT qf = __riscv_vfcvt_f (q, vlen);
  VFLOAT add_hi, add_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (VFMV_VF (RVV_ATANF_PI2_HI_F, vlen),
     VFMV_VF (RVV_ATANF_PI2_LO_F, vlen), qf,
     add_hi, add_lo, vlen);
  RVV_SLEEF_ORIG_DFADD_VF2_VF2 (add_hi, add_lo, t_hi, t_lo,
                                *r_hi, *r_lo, vlen);
}
#  undef RVV_SLEEF_ATAN_HELPERS_THIS_LMUL
# endif
#endif

/* expint1 is the only importer that enables this target-specific fixed-point
   and polynomial macro family.  */
#define RVV_MATH_IS_EXPINT1_expint1 1
#define RVV_MATH_IS_EXPINT1_INNER(name) RVV_MATH_IS_EXPINT1_##name
#define RVV_MATH_IS_EXPINT1(name) RVV_MATH_IS_EXPINT1_INNER (name)
#if RVV_MATH_IS_EXPINT1 (V_MATH_FUNC)
# if (BIT_WIDTH == 64)
#define RVV_EXPINT1_NEG_LOG2_HI_D -0x1.62e42fefa4000p-1
#define RVV_EXPINT1_NEG_LOG2_LO_D 0x1.8432a1b0e2634p-43
#define RVV_EXPINT1_DELTA_P0_D -0x1.991be2150638dp-91
#define RVV_EXPINT1_DELTA_Q0_D 0x1.91d57a67cdde2p-91
#define RVV_EXPINT1_P_ONE_Q62 RVV_ONE_Q62
#define RVV_EXPINT1_P_HALF_Q62 0x2000000000000000LL
#define RVV_EXPINT1_LOG2_HI_D RVV_LN2_HI_D
#define RVV_EXPINT1_LOG2_LO_D RVV_LN2_LO_D
#define RVV_EXPINT1_NEG_LOG2_INV_D -0x1.71547652b82fep+0
#define RVV_EXPINT1_NEAR0_SCALE_D 0x1.0p63
#define RVV_EXPINT1_NEAR0_INV_SCALE_D 0x1.0p-62
#define RVV_EXPINT1_RAT_SCALE_P75_D 0x1.0p-75
#define RVV_EXPINT1_RAT_SCALE_P81_D 0x1.0p-81
#define RVV_EXPINT1_RAT_SCALE_Q74_D 0x1.0p-74
#define RVV_EXPINT1_RAT_SCALE_Q81_D 0x1.0p-81
#define RVV_EXPINT1_RAT_SCALE_P89_D 0x1.0p-89
#define RVV_EXPINT1_RAT_SCALE_Q89_D 0x1.0p-89
#define RVV_EXPINT1_X_LT_ONE_SUBSTITUTE_D RVV_HALF_D
#define RVV_EXPINT1_X_GE_ONE_SUBSTITUTE_D 0x1.0p1
#define RVV_EXPINT1_X_GE_ONE_MAX_D 0x1.77p+9
#define RVV_EXPINT1_DENORM_SCALE_D 0x1.0p64
#define RVV_EXPINT1_DENORM_ADJUST 64
#define RVV_EXPINT1_EXPONENT_OFFSET 0x96
#define RVV_EXPINT1_RAT_SCALE_D 0x1.0p62
#define RVV_EXPINT1_P_NEAR0_Q62_0_I -0x24f119f8df6c31e9LL
#define RVV_EXPINT1_P_NEAR0_Q62_1_I 0x3fffffffffffffdcLL
#define RVV_EXPINT1_P_NEAR0_Q62_2_I -0x0ffffffffffff993LL
#define RVV_EXPINT1_P_NEAR0_Q62_3_I 0x038e38e38e3862efLL
#define RVV_EXPINT1_P_NEAR0_Q62_4_I -0x00aaaaaaaaa516c7LL
#define RVV_EXPINT1_P_NEAR0_0_D 0x1.b4e81b4c194fap-10
#define RVV_EXPINT1_P_NEAR0_1_D -0x1.e573ac379c696p-13
#define RVV_EXPINT1_P_NEAR0_2_D 0x1.db8b66c555673p-16
#define RVV_EXPINT1_P_NEAR0_3_D -0x1.a01962ee439b6p-19
#define RVV_EXPINT1_P_NEAR0_4_D 0x1.48bd8c51ff717p-22
#define RVV_EXPINT1_P_NEAR0_5_D -0x1.d8dbdf85f5051p-26
#define RVV_EXPINT1_P_NEAR0_6_D 0x1.355966c463c2ap-29
#define RVV_EXPINT1_P_NEAR0_7_D -0x1.5f2978a4477ccp-33
#define RVV_EXPINT1_P_NEAR0_8_D 0x1.0c38e425fee47p-37
#define RVV_EXPINT1_P0_I 0x01cd7ed8aff2c99aLL
#define RVV_EXPINT1_P1_I 0x569fe822aee57cb5LL
#define RVV_EXPINT1_P2_I 0x066ef748e71155e7LL
#define RVV_EXPINT1_P3_I 0x3eef1f3e5518e60cLL
#define RVV_EXPINT1_P4_I 0x0565f4de088a3f6fLL
#define RVV_EXPINT1_P5_I 0x110eb5a49f6eb3f7LL
#define RVV_EXPINT1_P6_I 0x1eb3d58b6063612aLL
#define RVV_EXPINT1_P7_I 0x1e310833e20e2c95LL
#define RVV_EXPINT1_P8_I 0x0ef22d4e45dc890fLL
#define RVV_EXPINT1_P9_I 0x0332908f3ee32e50LL
#define RVV_EXPINT1_P10_I 0x00363991b63aebbeLL
#define RVV_EXPINT1_P11_I 0x00003fcc5b1fe05fLL
#define RVV_EXPINT1_Q0_I 0x01cd7ed8aff2c99aLL
#define RVV_EXPINT1_Q1_I 0x586d66fb5ed84240LL
#define RVV_EXPINT1_Q2_I 0x06c3c9b231105002LL
#define RVV_EXPINT1_Q3_I 0x450cdf1ba384745aLL
#define RVV_EXPINT1_Q4_I 0x0325d39f7df69cefLL
#define RVV_EXPINT1_Q5_I 0x0adb448d52de9c87LL
#define RVV_EXPINT1_Q6_I 0x162a2f5f98fba589LL
#define RVV_EXPINT1_Q7_I 0x1a29d73cbd365659LL
#define RVV_EXPINT1_Q8_I 0x10f713cf2428b2ffLL
#define RVV_EXPINT1_Q9_I 0x0582b960921c6deeLL
#define RVV_EXPINT1_Q10_I 0x00c13a848700a0a5LL
#define RVV_EXPINT1_Q11_I 0x0007802b6e574e3eLL
#define RVV_EXPINT1_P_EXP3_D 0x1.555555555555ap+59
#define RVV_EXPINT1_P_EXP4_D 0x1.5555555555533p+57
#define RVV_EXPINT1_P_EXP5_D 0x1.111111110ef1dp+55
#define RVV_EXPINT1_P_EXP6_D 0x1.6c16c16c23cabp+52
#define RVV_EXPINT1_P_EXP7_D 0x1.a01a01b2eeafdp+49
#define RVV_EXPINT1_P_EXP8_D 0x1.a01a016c97838p+46
#define RVV_EXPINT1_P_EXP9_D 0x1.71ddf0af3f3a4p+43
#define RVV_EXPINT1_P_EXP10_D 0x1.27e542d471a01p+40
#define RVV_EXPINT1_P_EXP11_D 0x1.af6bfc694314ap+36
#define RVV_EXPINT1_P_EXP12_D 0x1.1ef1a5cf633bap+33
#define RVV_EXPINT1_NEG_LOG_POLY0_D 0x1.c71c51c73bb7fp-12
#define RVV_EXPINT1_NEG_LOG_POLY1_D 0x1.74664bed42062p-14
#define RVV_EXPINT1_NEG_LOG_POLY2_D 0x1.39a071f83b771p-16
#define RVV_EXPINT1_NEG_LOG_POLY3_D 0x1.2f123764244dfp-18
#define RVV_EXPINT1_NEG_LOG_POLY4_D 0x1.5555555555594p-4
#define RVV_EXPINT1_NEG_LOG_POLY5_D 0x1.999999997f6b6p-7
#define RVV_EXPINT1_NEG_LOG_POLY6_D 0x1.2492494248a48p-9
#undef RVV_EXPINT1_PSTEP_I_SLL
#undef RVV_EXPINT1_PSTEP_I
#undef RVV_EXPINT1_SMUL_Q63
#define RVV_EXPINT1_SMUL_Q63 __PASTE2 (rvv_expint1_smul_q63_m, LMUL)
#define RVV_EXPINT1_PSTEP_I __PASTE2 (rvv_expint1_pstep_i_m, LMUL)
#define RVV_EXPINT1_PSTEP_I_SLL __PASTE2 (rvv_expint1_pstep_i_sll_m, LMUL)
static inline VINT
RVV_EXPINT1_SMUL_Q63 (VINT a, VINT b, size_t vlen)
{
  VINT hi = __riscv_vmulh (a, b, vlen);
  VUINT lo = I_AS_U (__riscv_vmul (a, b, vlen));
  VUINT carry = __riscv_vsrl (lo, 63, vlen);
  return __riscv_vadd (__riscv_vsll (hi, 1, vlen), U_AS_I (carry), vlen);
}
static inline VINT
RVV_EXPINT1_PSTEP_I (int64_t coeff, VINT x, VINT poly, size_t vlen)
{
  VINT mul = RVV_EXPINT1_SMUL_Q63 (poly, x, vlen);
  return __riscv_vsadd (mul, coeff, vlen);
}
static inline VINT
RVV_EXPINT1_PSTEP_I_SLL (int64_t coeff, VINT x, unsigned int k, VINT poly,
                         size_t vlen)
{
  VINT mul = RVV_EXPINT1_SMUL_Q63 (poly, x, vlen);
  mul = __riscv_vsll (mul, k, vlen);
  return __riscv_vsadd (mul, coeff, vlen);
}
#undef PSTEP_I
#undef PSTEP_I_SLL
#define PSTEP_I(COEFF_j, X, POLY, vlen)                                       \
  RVV_EXPINT1_PSTEP_I ((COEFF_j), (X), (POLY), (vlen))
#define PSTEP_I_SLL(COEFF_j, X, K, POLY, vlen)                                \
  RVV_EXPINT1_PSTEP_I_SLL ((COEFF_j), (X), (K), (POLY), (vlen))
#ifndef V_D_EXPINT1_SKELETON_ONCE
#define V_D_EXPINT1_SKELETON_ONCE

#define NEG_LOG2_HI RVV_EXPINT1_NEG_LOG2_HI_D
#define NEG_LOG2_LO RVV_EXPINT1_NEG_LOG2_LO_D

#define EXCEPTION_HANDLING(vx, special_args, vy_special, n_adjust, vlen)      \
  do                                                                          \
    {                                                                         \
      VUINT vclass = __riscv_vfclass ((vx), (vlen));                          \
      IDENTIFY (vclass, 0x3BF, (special_args), (vlen));                       \
      UINT nb_special_args = __riscv_vcpop ((special_args), (vlen));          \
      (n_adjust) = VMVI_VX (0, (vlen));                                       \
      if (nb_special_args > 0)                                                \
        {                                                                     \
          VBOOL id_mask;                                                      \
          /* Substitute negative arguments with sNaN.  */                     \
          IDENTIFY (vclass, class_negative, id_mask, (vlen));                 \
          (vx) = __riscv_vfmerge ((vx), fp_sNaN, id_mask, (vlen));             \
          /* Substitute -0 argument with +0.  */                              \
          IDENTIFY (vclass, class_negZero, id_mask, (vlen));                  \
          (vx) = __riscv_vfmerge ((vx), fp_posZero, id_mask, (vlen));         \
          /* Eliminate positive denorm input from special_args.  */           \
          IDENTIFY (vclass, 0x39F, (special_args), (vlen));                   \
          /* For narrowed set of special arguments, compute vx+vfrec7(vx). */ \
          (vy_special) = __riscv_vfrec7 ((special_args), (vx), (vlen));       \
          (vy_special)                                                       \
              = __riscv_vfadd ((special_args), (vy_special), (vx), (vlen));   \
          (vx) = __riscv_vfmerge ((vx), fp_posOne, (special_args), (vlen));   \
          /* Scale up input for positive denormals.  */                       \
          IDENTIFY (vclass, class_posDenorm, id_mask, (vlen));                \
          (n_adjust) = __riscv_vmerge ((n_adjust), RVV_EXPINT1_DENORM_ADJUST, id_mask, (vlen));      \
          VFLOAT vx_normalized                                                \
              = __riscv_vfmul (id_mask, (vx), RVV_EXPINT1_DENORM_SCALE_D, (vlen));              \
          (vx) = __riscv_vmerge ((vx), vx_normalized, id_mask, (vlen));       \
          (vx) = __riscv_vfmerge ((vx), fp_posOne, (special_args), (vlen));   \
        }                                                                     \
    }                                                                         \
  while (0)

/* Polynomial coefficients for expint1 on [0,1]: -log(x) + poly(x) */
/* p0, ..., p4 in Q62 format */
#define P_near0_0 RVV_EXPINT1_P_NEAR0_Q62_0_I
#define P_near0_1 RVV_EXPINT1_P_NEAR0_Q62_1_I
#define P_near0_2 RVV_EXPINT1_P_NEAR0_Q62_2_I
#define P_near0_3 RVV_EXPINT1_P_NEAR0_Q62_3_I
#define P_near0_4 RVV_EXPINT1_P_NEAR0_Q62_4_I

/* p5, ..., p13 as floating point */
#define P_near0_5 RVV_EXPINT1_P_NEAR0_0_D
#define P_near0_6 RVV_EXPINT1_P_NEAR0_1_D
#define P_near0_7 RVV_EXPINT1_P_NEAR0_2_D
#define P_near0_8 RVV_EXPINT1_P_NEAR0_3_D
#define P_near0_9 RVV_EXPINT1_P_NEAR0_4_D
#define P_near0_10 RVV_EXPINT1_P_NEAR0_5_D
#define P_near0_11 RVV_EXPINT1_P_NEAR0_6_D
#define P_near0_12 RVV_EXPINT1_P_NEAR0_7_D
#define P_near0_13 RVV_EXPINT1_P_NEAR0_8_D

/* Rational function coefficients for x >= 1: exp(-x) * y * P(y)/Q(y), y=1/x */
#define P_0 RVV_EXPINT1_P0_I
#define P_1 RVV_EXPINT1_P1_I
#define P_2 RVV_EXPINT1_P2_I
#define P_3 RVV_EXPINT1_P3_I
#define P_4 RVV_EXPINT1_P4_I
#define P_5 RVV_EXPINT1_P5_I
#define P_6 RVV_EXPINT1_P6_I
#define P_7 RVV_EXPINT1_P7_I
#define P_8 RVV_EXPINT1_P8_I
#define P_9 RVV_EXPINT1_P9_I
#define P_10 RVV_EXPINT1_P10_I
#define P_11 RVV_EXPINT1_P11_I
#define delta_P_0 RVV_EXPINT1_DELTA_P0_D

#define Q_0 RVV_EXPINT1_Q0_I
#define Q_1 RVV_EXPINT1_Q1_I
#define Q_2 RVV_EXPINT1_Q2_I
#define Q_3 RVV_EXPINT1_Q3_I
#define Q_4 RVV_EXPINT1_Q4_I
#define Q_5 RVV_EXPINT1_Q5_I
#define Q_6 RVV_EXPINT1_Q6_I
#define Q_7 RVV_EXPINT1_Q7_I
#define Q_8 RVV_EXPINT1_Q8_I
#define Q_9 RVV_EXPINT1_Q9_I
#define Q_10 RVV_EXPINT1_Q10_I
#define Q_11 RVV_EXPINT1_Q11_I
#define delta_Q_0 RVV_EXPINT1_DELTA_Q0_D

/* Exponential polynomial coefficients scaled by 2^62 */
#define P_one RVV_EXPINT1_P_ONE_Q62
#define P_half RVV_EXPINT1_P_HALF_Q62
#define P_exp_3 RVV_EXPINT1_P_EXP3_D
#define P_exp_4 RVV_EXPINT1_P_EXP4_D
#define P_exp_5 RVV_EXPINT1_P_EXP5_D
#define P_exp_6 RVV_EXPINT1_P_EXP6_D
#define P_exp_7 RVV_EXPINT1_P_EXP7_D
#define P_exp_8 RVV_EXPINT1_P_EXP8_D
#define P_exp_9 RVV_EXPINT1_P_EXP9_D
#define P_exp_10 RVV_EXPINT1_P_EXP10_D
#define P_exp_11 RVV_EXPINT1_P_EXP11_D
#define P_exp_12 RVV_EXPINT1_P_EXP12_D

#define LOG2_HI RVV_EXPINT1_LOG2_HI_D
#define LOG2_LO RVV_EXPINT1_LOG2_LO_D
#define NEG_LOG2_INV RVV_EXPINT1_NEG_LOG2_INV_D

#define EXPINT1_NEG_LOG(vx_in, n_adjust, y_hi, y_lo, vlen)                    \
  do                                                                          \
    {                                                                         \
      VFLOAT vx = (vx_in);                                                    \
      VINT n = U_AS_I (__riscv_vadd (                                         \
          __riscv_vsrl (F_AS_U (vx), MAN_LEN - 8, (vlen)), RVV_EXPINT1_EXPONENT_OFFSET, (vlen)));    \
      n = __riscv_vsra (n, 8, (vlen));                                        \
      n = __riscv_vsub (n, EXP_BIAS, (vlen));                                 \
      vx = I_AS_F (                                                           \
          __riscv_vsub (F_AS_I (vx), __riscv_vsll (n, MAN_LEN, (vlen)),        \
                        (vlen)));                                             \
      n = __riscv_vsub (n, (n_adjust), (vlen));                               \
      VFLOAT n_flt = __riscv_vfcvt_f (n, (vlen));                             \
                                                                              \
      VFLOAT numer = __riscv_vfsub (vx, fp_posOne, (vlen));                   \
      numer = __riscv_vfadd (numer, numer, (vlen));                           \
      VFLOAT denom = __riscv_vfadd (vx, fp_posOne, (vlen));                   \
      VFLOAT delta_d = __riscv_vfrsub (denom, fp_posOne, (vlen));             \
      delta_d = __riscv_vfadd (delta_d, vx, (vlen));                          \
      VFLOAT r, delta_r;                                                      \
      DIV_N1D2 (numer, denom, delta_d, r, delta_r, (vlen));                   \
                                                                              \
      VFLOAT rsq = __riscv_vfmul (r, r, (vlen));                              \
      VFLOAT rcube = __riscv_vfmul (rsq, r, (vlen));                          \
      VFLOAT r6 = __riscv_vfmul (rcube, rcube, (vlen));                       \
                                                                              \
      VFLOAT poly_right = PSTEP (                                             \
          RVV_EXPINT1_NEG_LOG_POLY0_D, rsq,                                    \
          PSTEP (RVV_EXPINT1_NEG_LOG_POLY1_D, rsq,                             \
                 PSTEP (RVV_EXPINT1_NEG_LOG_POLY2_D,                            \
                        RVV_EXPINT1_NEG_LOG_POLY3_D, rsq,                      \
                        (vlen)),                                             \
                 (vlen)),                                                    \
          (vlen));                                                           \
                                                                              \
      VFLOAT poly_left =                                                      \
          PSTEP (RVV_EXPINT1_NEG_LOG_POLY4_D, rsq,                              \
                 PSTEP (RVV_EXPINT1_NEG_LOG_POLY5_D,                            \
                        RVV_EXPINT1_NEG_LOG_POLY6_D, rsq,                       \
                        (vlen)),                                             \
                 (vlen));                                                    \
                                                                              \
      VFLOAT poly = __riscv_vfmadd (poly_right, r6, poly_left, (vlen));       \
      delta_r = __riscv_vfmsac (delta_r, NEG_LOG2_LO, n_flt, (vlen));         \
      poly = __riscv_vfnmsub (poly, rcube, delta_r, (vlen));                  \
                                                                              \
      VFLOAT A = __riscv_vfmul (n_flt, NEG_LOG2_HI, (vlen));                  \
      r = __riscv_vfsgnjx (r, fp_negOne, (vlen));                             \
      FAST2SUM (A, r, (y_hi), (y_lo), (vlen));                                \
      (y_lo) = __riscv_vfadd ((y_lo), poly, (vlen));                          \
    }                                                                         \
  while (0)

#define EXPINT_POLY_01(vx, y_hi, y_lo, vlen)                                  \
  do                                                                          \
    {                                                                         \
      VINT P_FIX, X;                                                          \
      X = __riscv_vfcvt_x (__riscv_vfmul ((vx), RVV_EXPINT1_NEAR0_SCALE_D, (vlen)), (vlen));   \
      P_FIX = PSTEP_I (                                                       \
          P_near0_0, X,                                                       \
          PSTEP_I (P_near0_1, X,                                              \
                   PSTEP_I (P_near0_2, X,                                     \
                            PSTEP_I (P_near0_3, X,                            \
                                     VMVI_VX (P_near0_4, (vlen)), (vlen)),    \
                            (vlen)),                                         \
                   (vlen)),                                                  \
          (vlen));                                                           \
      VFLOAT p_left_hi = __riscv_vfcvt_f (P_FIX, (vlen));                     \
      P_FIX = __riscv_vsub (P_FIX, __riscv_vfcvt_x (p_left_hi, (vlen)),        \
                            (vlen));                                         \
      VFLOAT p_left_lo = __riscv_vfcvt_f (P_FIX, (vlen));                     \
      p_left_lo = __riscv_vfmul (p_left_lo, RVV_EXPINT1_NEAR0_INV_SCALE_D, (vlen));               \
      p_left_hi = __riscv_vfmul (p_left_hi, RVV_EXPINT1_NEAR0_INV_SCALE_D, (vlen));               \
                                                                              \
      VFLOAT poly_right =                                                    \
          PSTEP (P_near0_9, (vx),                                             \
                 PSTEP (P_near0_10, (vx),                                     \
                        PSTEP (P_near0_11, (vx),                              \
                               PSTEP (P_near0_12, P_near0_13, (vx), (vlen)),  \
                               (vlen)),                                      \
                        (vlen)),                                             \
                 (vlen));                                                    \
                                                                              \
      VFLOAT poly_mid =                                                      \
          PSTEP (P_near0_5, (vx),                                             \
                 PSTEP (P_near0_6, (vx),                                      \
                        PSTEP (P_near0_7, P_near0_8, (vx), (vlen)), (vlen)),  \
                 (vlen));                                                    \
      VFLOAT x4 = __riscv_vfmul ((vx), (vx), (vlen));                         \
      x4 = __riscv_vfmul (x4, x4, (vlen));                                    \
      VFLOAT x5 = __riscv_vfmul ((vx), x4, (vlen));                           \
      poly_mid = __riscv_vfmacc (poly_mid, x4, poly_right, (vlen));           \
      p_left_lo = __riscv_vfmacc (p_left_lo, x5, poly_mid, (vlen));           \
      KNUTH2SUM (p_left_hi, p_left_lo, (y_hi), (y_lo), (vlen));               \
    }                                                                         \
  while (0)

#define EXPINT1_RAT_GE1(x_hi, x_lo, scale, rat_hi, rat_lo, vlen)              \
  do                                                                          \
    {                                                                         \
      VINT P75, P81, P89, Q74, Q81, Q89;                                      \
      VINT _X;                                                                \
      FLT2FIX ((x_hi), (x_lo), (scale), _X, (vlen));                          \
      P75 = PSTEP_I_SLL (P_10, _X, 1, VMVI_VX (P_11, (vlen)), (vlen));         \
      P75 = PSTEP_I_SLL (P_9, _X, 1, P75, (vlen));                            \
      P75 = PSTEP_I_SLL (P_8, _X, 1, P75, (vlen));                            \
      P75 = PSTEP_I_SLL (P_7, _X, 1, P75, (vlen));                            \
      P75 = PSTEP_I_SLL (P_6, _X, 1, P75, (vlen));                            \
      P75 = PSTEP_I_SLL (P_5, _X, 1, P75, (vlen));                            \
      P75 = PSTEP_I_SLL (P_4, _X, 1, P75, (vlen));                            \
      VFLOAT _xsq_hi, _xsq_lo;                                                \
      SQR_X2 ((x_hi), (x_lo), _xsq_hi, _xsq_lo, (vlen));                      \
      VFLOAT _p_right_hi, _p_right_lo;                                        \
      FIX2FLT (P75, RVV_EXPINT1_RAT_SCALE_P75_D, _p_right_hi, _p_right_lo, (vlen));             \
      VFLOAT _p_tmp1_hi, _p_tmp1_lo;                                          \
      PROD_X2Y2 (_xsq_hi, _xsq_lo, _p_right_hi, _p_right_lo, _p_tmp1_hi,       \
                 _p_tmp1_lo, (vlen));                                         \
                                                                              \
      P81 = PSTEP_I_SLL (P_2, _X, 1, VMVI_VX (P_3, (vlen)), (vlen));           \
      VFLOAT _p_mid_hi, _p_mid_lo;                                            \
      FIX2FLT (P81, RVV_EXPINT1_RAT_SCALE_P81_D, _p_mid_hi, _p_mid_lo, (vlen));                 \
      VFLOAT _p_tmp2_hi, _p_tmp2_lo;                                          \
      POS2SUM (_p_tmp1_hi, _p_mid_hi, _p_tmp2_hi, _p_tmp2_lo, (vlen));        \
      _p_tmp2_lo = __riscv_vfadd (_p_tmp2_lo, _p_tmp1_lo, (vlen));            \
      _p_tmp2_lo = __riscv_vfadd (_p_tmp2_lo, _p_mid_lo, (vlen));             \
                                                                              \
      PROD_X2Y2 (_xsq_hi, _xsq_lo, _p_tmp2_hi, _p_tmp2_lo, _p_tmp1_hi,         \
                 _p_tmp1_lo, (vlen));                                         \
      VFLOAT _p_left_hi, _p_left_lo;                                          \
      P89 = PSTEP_I_SLL (P_0, _X, 1, VMVI_VX (P_1, (vlen)), (vlen));           \
      FIX2FLT (P89, RVV_EXPINT1_RAT_SCALE_P89_D, _p_left_hi, _p_left_lo, (vlen));               \
      POS2SUM (_p_left_hi, _p_tmp1_hi, _p_tmp2_hi, _p_tmp2_lo, (vlen));       \
      _p_tmp2_lo = __riscv_vfadd (_p_tmp2_lo, _p_left_lo, (vlen));            \
      _p_tmp2_lo = __riscv_vfadd (_p_tmp2_lo, _p_tmp1_lo, (vlen));            \
      _p_tmp2_lo = __riscv_vfadd (_p_tmp2_lo, delta_P_0, (vlen));             \
      VFLOAT _p_hi, _p_lo;                                                    \
      FAST2SUM (_p_tmp2_hi, _p_tmp2_lo, _p_hi, _p_lo, (vlen));                \
                                                                              \
      Q74 = PSTEP_I_SLL (Q_10, _X, 1, VMVI_VX (Q_11, (vlen)), (vlen));         \
      Q74 = PSTEP_I_SLL (Q_9, _X, 1, Q74, (vlen));                            \
      Q74 = PSTEP_I_SLL (Q_8, _X, 1, Q74, (vlen));                            \
      Q74 = PSTEP_I_SLL (Q_7, _X, 1, Q74, (vlen));                            \
      Q74 = PSTEP_I_SLL (Q_6, _X, 1, Q74, (vlen));                            \
      Q74 = PSTEP_I_SLL (Q_5, _X, 1, Q74, (vlen));                            \
      Q74 = PSTEP_I_SLL (Q_4, _X, 1, Q74, (vlen));                            \
                                                                              \
      VFLOAT _q_right_hi, _q_right_lo;                                        \
      FIX2FLT (Q74, RVV_EXPINT1_RAT_SCALE_Q74_D, _q_right_hi, _q_right_lo, (vlen));             \
      VFLOAT _q_tmp1_hi, _q_tmp1_lo;                                          \
      PROD_X2Y2 (_xsq_hi, _xsq_lo, _q_right_hi, _q_right_lo, _q_tmp1_hi,       \
                 _q_tmp1_lo, (vlen));                                         \
                                                                              \
      Q81 = PSTEP_I_SLL (Q_2, _X, 1, VMVI_VX (Q_3, (vlen)), (vlen));           \
      VFLOAT _q_mid_hi, _q_mid_lo;                                            \
      FIX2FLT (Q81, RVV_EXPINT1_RAT_SCALE_Q81_D, _q_mid_hi, _q_mid_lo, (vlen));                 \
      VFLOAT _q_tmp2_hi, _q_tmp2_lo;                                          \
      POS2SUM (_q_tmp1_hi, _q_mid_hi, _q_tmp2_hi, _q_tmp2_lo, (vlen));        \
      _q_tmp2_lo = __riscv_vfadd (_q_tmp2_lo, _q_tmp1_lo, (vlen));            \
      _q_tmp2_lo = __riscv_vfadd (_q_tmp2_lo, _q_mid_lo, (vlen));             \
                                                                              \
      PROD_X2Y2 (_xsq_hi, _xsq_lo, _q_tmp2_hi, _q_tmp2_lo, _q_tmp1_hi,         \
                 _q_tmp1_lo, (vlen));                                         \
      VFLOAT _q_left_hi, _q_left_lo;                                          \
      Q89 = PSTEP_I_SLL (Q_0, _X, 1, VMVI_VX (Q_1, (vlen)), (vlen));           \
      FIX2FLT (Q89, RVV_EXPINT1_RAT_SCALE_Q89_D, _q_left_hi, _q_left_lo, (vlen));               \
      POS2SUM (_q_left_hi, _q_tmp1_hi, _q_tmp2_hi, _q_tmp2_lo, (vlen));       \
      _q_tmp2_lo = __riscv_vfadd (_q_tmp2_lo, _q_left_lo, (vlen));            \
      _q_tmp2_lo = __riscv_vfadd (_q_tmp2_lo, _q_tmp1_lo, (vlen));            \
      _q_tmp2_lo = __riscv_vfadd (_q_tmp2_lo, delta_Q_0, (vlen));             \
      VFLOAT _q_hi, _q_lo;                                                    \
      FAST2SUM (_q_tmp2_hi, _q_tmp2_lo, _q_hi, _q_lo, (vlen));                \
                                                                              \
      DIV2_N2D2 (_p_hi, _p_lo, _q_hi, _q_lo, (rat_hi), (rat_lo), (vlen));      \
    }                                                                         \
  while (0)

#define EXP_NEGX(vx, n, y_hi, y_lo, vlen)                                     \
  do                                                                          \
    {                                                                         \
      VFLOAT _n_flt = __riscv_vfmul ((vx), NEG_LOG2_INV, (vlen));             \
      (n) = __riscv_vfcvt_x (_n_flt, (vlen));                                 \
      _n_flt = __riscv_vfcvt_f ((n), (vlen));                                 \
      VFLOAT _r_hi = __riscv_vfnmadd (_n_flt, LOG2_HI, (vx), (vlen));          \
      VFLOAT _r_lo = __riscv_vfmul (_n_flt, LOG2_LO, (vlen));                 \
      VFLOAT _r = __riscv_vfsub (_r_hi, _r_lo, (vlen));                       \
      _r_lo = __riscv_vfsgnjx (_r_lo, fp_negOne, (vlen));                     \
      VFLOAT _p_right = PSTEP (                                               \
          P_exp_9, _r,                                                        \
          PSTEP (P_exp_10, _r, PSTEP (P_exp_11, P_exp_12, _r, (vlen)),         \
                 (vlen)),                                                    \
          (vlen));                                                           \
      VFLOAT _p_mid =                                                        \
          PSTEP (P_exp_5, _r,                                                 \
                 PSTEP (P_exp_6, _r, PSTEP (P_exp_7, P_exp_8, _r, (vlen)),     \
                        (vlen)),                                             \
                 (vlen));                                                    \
      VFLOAT _rsq = __riscv_vfmul (_r, _r, (vlen));                           \
      VFLOAT _r4 = __riscv_vfmul (_rsq, _rsq, (vlen));                        \
      _p_mid = __riscv_vfmacc (_p_mid, _r4, _p_right, (vlen));                \
      VFLOAT _p_left = PSTEP (P_exp_3, P_exp_4, _r, (vlen));                  \
      _p_left = __riscv_vfmacc (_p_left, _rsq, _p_mid, (vlen));               \
      VINT _P = __riscv_vfcvt_x (_p_left, (vlen));                            \
      VINT _R;                                                               \
      FLT2FIX (_r_hi, _r_lo, RVV_EXPINT1_NEAR0_SCALE_D, _R, (vlen));                           \
      _P = PSTEP_I (P_one, _R,                                                \
                    PSTEP_I (P_one, _R, PSTEP_I (P_half, _R, _P, (vlen)),      \
                             (vlen)),                                        \
                    (vlen));                                                 \
      FIX2FLT (_P, RVV_EXPINT1_NEAR0_INV_SCALE_D, (y_hi), (y_lo), (vlen));                        \
    }                                                                         \
  while (0)


#endif

# else
# endif
#endif
#undef RVV_MATH_IS_EXPINT1
#undef RVV_MATH_IS_EXPINT1_INNER
#undef RVV_MATH_IS_EXPINT1_expint1
/* The double-precision cbrt importer keeps its source EH envelope in this
   shared helper so the skeleton contains only the operation graph and its
   data reference.  */
#ifndef RVV_MATH_IS_CBRT_cbrt
# define RVV_MATH_IS_CBRT_cbrt 1
# define RVV_MATH_IS_CBRT_INNER(name) RVV_MATH_IS_CBRT_##name
# define RVV_MATH_IS_CBRT(name) RVV_MATH_IS_CBRT_INNER (name)
#endif

#if RVV_MATH_IS_CBRT (V_MATH_FUNC)
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_exception_handling, LMUL)
  (VFLOAT *vx, VBOOL *special_args, VFLOAT *vy_special, VINT *n_adjust,
   size_t vlen)
{
  VUINT vclass = __riscv_vfclass (*vx, vlen);
  IDENTIFY (vclass, 0x399, *special_args, vlen);
  VBOOL denorm;
  IDENTIFY (vclass, 0x24, denorm, vlen);
  VBOOL special_n_denorm = __riscv_vmor (*special_args, denorm, vlen);
  *n_adjust = __riscv_vxor (*n_adjust, *n_adjust, vlen);
  if (__riscv_vcpop (special_n_denorm, vlen) > 0)
    {
      VFLOAT vx_normal = __riscv_vfmul
        (denorm, *vx, VFMV_VF (0x1.0p60, vlen), vlen);
      *vx = __riscv_vmerge (*vx, vx_normal, denorm, vlen);
      *n_adjust = __riscv_vmerge (*n_adjust, -20, denorm, vlen);
      *vy_special = __riscv_vfadd (*special_args, *vx, *vx, vlen);
      *vx = __riscv_vfmerge (*vx, fp_posOne, *special_args, vlen);
    }
}

/* Keep the cbrt polynomial operation graph unchanged for active finite lanes
   while preventing sanitized exceptional lanes from executing the FMA.  The
   initial scalar polynomial constants are materialized by the caller so this
   helper uses the same vector-vector fused operation as the common PSTEP
   spelling.  */
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_cbrt_pstep_masked, LMUL)
  (VBOOL active, double coeff, VFLOAT x, VFLOAT poly, size_t vlen)
{
  return VFMADD_VV_M (active, poly, x, VFMV_VF (coeff, vlen), vlen);
}
#endif

/* The double-precision erfc importer keeps its source helper operation
   graphs here so the skeleton contains only the main computation and its
   data reference.  */
#define RVV_MATH_IS_ERFC_erfc 1
#define RVV_MATH_IS_ERFC_INNER(name) RVV_MATH_IS_ERFC_##name
#define RVV_MATH_IS_ERFC(name) RVV_MATH_IS_ERFC_INNER (name)

#if (BIT_WIDTH == 64) && RVV_MATH_IS_ERFC (V_MATH_FUNC)
# define RVV_ERFC_PSTEP_I __PASTE2 (rvv_erfc_pstep_i_m, LMUL)
# define RVV_ERFC_PSTEP_I_SCALAR \
  __PASTE2 (rvv_erfc_pstep_i_scalar_m, LMUL)
# define RVV_ERFC_EXCEPTION_HANDLING \
  __PASTE2 (rvv_erfc_exception_handling_m, LMUL)
# define RVV_ERFC_EXP_NEG_AB __PASTE2 (rvv_erfc_exp_neg_ab_m, LMUL)
# define RVV_ERFC_X_TRANSFORM __PASTE2 (rvv_erfc_x_transform_m, LMUL)
# define RVV_ERFC_RECIP_SCALE __PASTE2 (rvv_erfc_recip_scale_m, LMUL)

static inline VINT
RVV_ERFC_PSTEP_I (int64_t coeff, VINT x, VINT poly, size_t vlen)
{
  return __riscv_vsadd (__riscv_vsmul (poly, x, 1, vlen), coeff, vlen);
}

static inline VINT
RVV_ERFC_PSTEP_I_SCALAR (int64_t coeff, int64_t x, VINT poly, size_t vlen)
{
  return __riscv_vsadd (__riscv_vsmul (poly, x, 1, vlen), coeff, vlen);
}

static inline void
RVV_ERFC_EXCEPTION_HANDLING (VFLOAT *vx,
                             VBOOL *special_args, VFLOAT *vy_special,
                             size_t vlen)
{
  VUINT vclass = __riscv_vfclass (*vx, vlen);
  IDENTIFY (vclass, class_NaN | class_Inf, *special_args, vlen);
  if (__riscv_vcpop (*special_args, vlen) > 0)
    {
      VBOOL id_mask;
      IDENTIFY (vclass, class_NaN, id_mask, vlen);
      VFLOAT nan_val = __riscv_vfadd (id_mask, *vx, *vx, vlen);
      *vy_special = __riscv_vmerge (*vy_special, nan_val, id_mask, vlen);
      IDENTIFY (vclass, class_posInf, id_mask, vlen);
      *vy_special = __riscv_vfmerge (*vy_special, fp_posZero, id_mask, vlen);
      IDENTIFY (vclass, class_negInf, id_mask, vlen);
      *vy_special = __riscv_vfmerge (*vy_special, RVV_TWO_D, id_mask, vlen);
      *vx = __riscv_vfmerge (*vx, fp_posZero, *special_args, vlen);
    }
}

static inline void
RVV_ERFC_EXP_NEG_AB (VFLOAT va, VFLOAT vb, VINT *n, VINT *z, size_t vlen)
{
  VFLOAT r = __riscv_vfmul (va, vb, vlen);
  VFLOAT delta_r = __riscv_vfmsub (va, vb, r, vlen);
  VFLOAT n_flt = __riscv_vfmul (r, -RVV_LN2_INV_D, vlen);
  *n = __riscv_vfcvt_x (n_flt, vlen);
  n_flt = __riscv_vfcvt_f (*n, vlen);
  r = __riscv_vfnmacc (r, RVV_LN2_HI_D, n_flt, vlen);
  delta_r = __riscv_vfnmacc (delta_r, RVV_LN2_LO_D, n_flt, vlen);

  VINT r_fixed = __riscv_vfcvt_x
    (__riscv_vfmul (r, RVV_Q63_SCALE_D, vlen), vlen);
  VINT delta_r_fixed = __riscv_vfcvt_x
    (__riscv_vfmul (delta_r, RVV_Q63_SCALE_D, vlen), vlen);
  r_fixed = __riscv_vadd (r_fixed, delta_r_fixed, vlen);

  VINT p_right = RVV_ERFC_PSTEP_I
    (RVV_ERFC_EXP_P_RIGHT0_I, r_fixed,
     RVV_ERFC_PSTEP_I
       (RVV_ERFC_EXP_P_RIGHT1_I, r_fixed,
        RVV_ERFC_PSTEP_I
          (RVV_ERFC_EXP_P_RIGHT2_I, r_fixed,
              RVV_ERFC_PSTEP_I
               (RVV_ERFC_EXP_P_RIGHT3_I, r_fixed,
               RVV_ERFC_PSTEP_I_SCALAR (RVV_ERFC_EXP_P_RIGHT4_I,
                                       RVV_ERFC_EXP_P_RIGHT_BASE_I, r_fixed,
                                       vlen),
              vlen),
           vlen),
        vlen),
     vlen);
  VINT r_sq = __riscv_vsmul (r_fixed, r_fixed, 1, vlen);
  VINT r_6 = __riscv_vsmul (r_sq, r_sq, 1, vlen);
  r_6 = __riscv_vsmul (r_6, r_sq, 1, vlen);
  VINT p_left = RVV_ERFC_PSTEP_I
    (RVV_ERFC_EXP_P_LEFT0_I, r_fixed,
     RVV_ERFC_PSTEP_I
       (RVV_ERFC_EXP_P_LEFT1_I, r_fixed,
        RVV_ERFC_PSTEP_I
          (RVV_ERFC_EXP_P_LEFT2_I, r_fixed,
             RVV_ERFC_PSTEP_I
               (RVV_ERFC_EXP_P_LEFT3_I, r_fixed,
               RVV_ERFC_PSTEP_I_SCALAR (RVV_ERFC_EXP_P_LEFT4_I,
                                       RVV_ERFC_EXP_P_LEFT_BASE_I, r_fixed,
                                       vlen),
              vlen),
           vlen),
        vlen),
     vlen);
  p_right = __riscv_vsmul (p_right, r_6, 1, vlen);
  p_right = __riscv_vadd (p_left, p_right, vlen);
  *z = p_right;
}

static inline VINT
RVV_ERFC_X_TRANSFORM (VFLOAT vx, size_t vlen)
{
  VFLOAT denom = __riscv_vfadd (vx, RVV_ERFC_B_FOR_TRANS_D, vlen);
  VFLOAT delta_denom = __riscv_vfrsub (denom, RVV_ERFC_B_FOR_TRANS_D, vlen);
  delta_denom = __riscv_vfadd (delta_denom, vx, vlen);
  VFLOAT one = VFMV_VF (fp_posOne, vlen);
  VFLOAT recip, delta_recip;
  DIV_N1D2 (one, denom, delta_denom, recip, delta_recip, vlen);
  VFLOAT numer = __riscv_vfmul (vx, RVV_Q63_SCALE_D, vlen);
  numer = __riscv_vfadd (numer, RVV_ERFC_NEG_A_SCALED_D, vlen);
  VFLOAT delta_numer = __riscv_vfrsub (numer, RVV_ERFC_NEG_A_SCALED_D, vlen);
  delta_numer = __riscv_vfmacc (delta_numer, RVV_Q63_SCALE_D, vx, vlen);
  VFLOAT q = __riscv_vfmul (numer, recip, vlen);
  VFLOAT delta_q = __riscv_vfmsub (numer, recip, q, vlen);
  delta_q = __riscv_vfmacc (delta_q, numer, delta_recip, vlen);
  delta_q = __riscv_vfmacc (delta_q, delta_numer, recip, vlen);
  VINT r = __riscv_vfcvt_x (q, vlen);
  return __riscv_vadd (r, __riscv_vfcvt_x (delta_q, vlen), vlen);
}

static inline void
RVV_ERFC_RECIP_SCALE (VFLOAT vx, VINT *b, VINT *m, size_t vlen)
{
  VFLOAT one = VFMV_VF (fp_posOne, vlen);
  VFLOAT denom = __riscv_vfmadd (vx, RVV_TWO_D, one, vlen);
  VFLOAT delta_denom = __riscv_vfsub (one, denom, vlen);
  delta_denom = __riscv_vfmacc (delta_denom, RVV_TWO_D, vx, vlen);
  VFLOAT recip, delta_recip;
  DIV_N1D2 (one, denom, delta_denom, recip, delta_recip, vlen);
  *m = __riscv_vsra (F_AS_I (vx), MAN_LEN, vlen);
  *m = __riscv_vmax (*m, EXP_BIAS - 1, vlen);
  *m = __riscv_vadd (*m, 63, vlen);
  VFLOAT scale = I_AS_F (__riscv_vsll (*m, MAN_LEN, vlen));
  *m = __riscv_vsub (*m, EXP_BIAS, vlen);
  *b = __riscv_vfcvt_x (__riscv_vfmul (recip, scale, vlen), vlen);
  delta_recip = __riscv_vfmul (delta_recip, scale, vlen);
  *b = __riscv_vadd (*b, __riscv_vfcvt_x (delta_recip, vlen), vlen);
}
#endif

#undef RVV_MATH_IS_ERFC
#undef RVV_MATH_IS_ERFC_INNER
#undef RVV_MATH_IS_ERFC_erfc

/* The double-precision erfcinv importer keeps its source helper operation
   graphs here so the skeleton contains only the main computation and its
   data reference.  */
#define RVV_MATH_IS_ERFCINV_erfcinv 1
#define RVV_MATH_IS_ERFCINV_INNER(name) RVV_MATH_IS_ERFCINV_##name
#define RVV_MATH_IS_ERFCINV(name) RVV_MATH_IS_ERFCINV_INNER (name)

#if (BIT_WIDTH == 64) && RVV_MATH_IS_ERFCINV (V_MATH_FUNC)
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_exception_handling, LMUL)
  (VFLOAT *vx, VBOOL *special_args, VFLOAT *vy_special, size_t vlen)
{
  VUINT vclass = __riscv_vfclass (*vx, vlen);
  IDENTIFY (vclass, 0x39F, *special_args, vlen);
  VBOOL x_ge_2 = __riscv_vmfge (*vx, RVV_TWO_D, vlen);
  *special_args = __riscv_vmor (*special_args, x_ge_2, vlen);

  if (__riscv_vcpop (*special_args, vlen) > 0)
    {
      VBOOL x_gt_2 = __riscv_vmfgt (*vx, RVV_TWO_D, vlen);
      VBOOL x_lt_0 = __riscv_vmflt (*vx, fp_posZero, vlen);
      *vx = VFMERGE_VFM (*vx, fp_sNaN, x_gt_2, vlen);
      *vx = VFMERGE_VFM (*vx, fp_sNaN, x_lt_0, vlen);

      VFLOAT tmp = VFMV_VF (fp_posZero, vlen);
      VFLOAT x_tmp = __riscv_vfrsub (*vx, fp_posOne, vlen);
      tmp = __riscv_vfsgnj (tmp, x_tmp, vlen);

      VBOOL x_eq_2 = __riscv_vmfeq (*vx, RVV_TWO_D, vlen);
      VBOOL x_eq_0 = __riscv_vmfeq (*vx, fp_posZero, vlen);
      VBOOL pm_inf = __riscv_vmor (x_eq_2, x_eq_0, vlen);
      tmp = VFREC7_M (pm_inf, tmp, vlen);

      /* Keep the clean definedness/operation-scope adaptation: the source
         spells this exceptional base as a masked vfsub, while the clean
         importer computes it unmasked and selects the exceptional lanes at
         the final merge.  */
      *vy_special = __riscv_vfsub (*vx, *vx, vlen);
      *vy_special = VMERGE_VVM_F (*vy_special, tmp, pm_inf, vlen);

      *vx = VFMERGE_VFM (*vx, fp_posOne, *special_args, vlen);
    }
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_neg_logx, LMUL)
  (VFLOAT vx, int n_adjust,
   VFLOAT *y_hi, VFLOAT *y_lo, size_t vlen)
{
  VFLOAT vx_in = vx;
  VINT n
    = __riscv_vadd (__riscv_vsra (F_AS_I (vx_in), MAN_LEN - 8, vlen), 0x96,
                    vlen);
  n = __riscv_vsub (__riscv_vsra (n, 8, vlen), EXP_BIAS, vlen);
  VFLOAT scale
    = I_AS_F (__riscv_vsll (__riscv_vrsub (n, EXP_BIAS, vlen), MAN_LEN,
                            vlen));
  vx_in = __riscv_vfmul (vx_in, scale, vlen);
  n = __riscv_vsub (n, n_adjust, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT numer = __riscv_vfrsub (vx_in, fp_posOne, vlen);
  numer = __riscv_vfadd (numer, numer, vlen);
  VFLOAT denom = __riscv_vfadd (vx_in, fp_posOne, vlen);
  VFLOAT delta_denom
    = __riscv_vfadd (__riscv_vfrsub (denom, fp_posOne, vlen), vx_in, vlen);
  VFLOAT r_hi, r_lo, r;
  DIV_N1D2 (numer, denom, delta_denom, r_hi, r_lo, vlen);
  r = __riscv_vfadd (r_hi, r_lo, vlen);
  VFLOAT a = __riscv_vfmadd (n_flt, -RVV_LOG2_HI_D, r_hi, vlen);
  VFLOAT b = __riscv_vfmsub (n_flt, -RVV_LOG2_HI_D, a, vlen);
  b = __riscv_vfadd (b, r_hi, vlen);
  b = __riscv_vfadd (r_lo, b, vlen);
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r6 = __riscv_vfmul (rcube, rcube, vlen);
  VFLOAT poly_right
    = PSTEP (RVV_ERFCINV_LOG_POLY_RIGHT0_D, rsq,
             PSTEP (RVV_ERFCINV_LOG_POLY_RIGHT1_D,
                    RVV_ERFCINV_LOG_POLY_RIGHT2_D, rsq,
                    vlen),
             vlen);
  VFLOAT poly_left
    = PSTEP (RVV_ERFCINV_LOG_POLY_LEFT0_D, rsq,
             PSTEP (RVV_ERFCINV_LOG_POLY_LEFT1_D,
                    RVV_ERFCINV_LOG_POLY_LEFT2_D, rsq,
                    vlen),
             vlen);
  poly_left = __riscv_vfmacc (poly_left, r6, poly_right, vlen);
  poly_left = PSTEP (RVV_ERFCINV_LOG_POLY_LEFT3_D, rsq, poly_left, vlen);
  b = __riscv_vfmacc (b, -RVV_LOG2_LO_D, n_flt, vlen);
  b = __riscv_vfmacc (b, rcube, poly_left, vlen);
  FAST2SUM (a, b, *y_hi, *y_lo, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_sqrtx, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT *w_hi, VFLOAT *w_lo, VINT *t,
   double t_scale, double t_scale_inv, size_t vlen)
{
  *w_hi = __riscv_vfsqrt (y_hi, vlen);
  *w_lo = __riscv_vfnmsub (*w_hi, *w_hi, y_hi, vlen);
  *w_lo = __riscv_vfadd (*w_lo, y_lo, vlen);
  VFLOAT recip = __riscv_vfadd (y_hi, y_hi, vlen);
  recip = __riscv_vfrec7 (recip, vlen);
  recip = __riscv_vfmul (recip, *w_hi, vlen);
  *w_lo = __riscv_vfmul (*w_lo, recip, vlen);
  VFLOAT t_lo = VFMV_VF (t_scale, vlen);
  VFLOAT t_hi = __riscv_vfdiv (t_lo, *w_hi, vlen);
  *t = __riscv_vfcvt_x (t_hi, vlen);
  t_lo = __riscv_vfnmsac (t_lo, *w_hi, t_hi, vlen);
  t_lo = __riscv_vfnmsac (t_lo, *w_lo, t_hi, vlen);
  t_lo = __riscv_vfmul (t_lo, t_hi, vlen);
  t_lo = __riscv_vfmul (t_lo, t_scale_inv, vlen);
  *t = __riscv_vadd (*t, __riscv_vfcvt_x (t_lo, vlen), vlen);
}

static inline VINT
RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
  (int64_t coeff, VINT t, VINT p_hi, VINT p_lo, size_t vlen)
{
  p_lo = PSTEP_I_HI (coeff, t, p_lo, vlen);
  VINT p_tmp = __riscv_vmul (t, p_hi, vlen);
  p_tmp = VSRL_I_AS_U (p_tmp, 30, vlen);
  return __riscv_vadd (p_lo, p_tmp, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_pq_hilo_tiny, LMUL)
  (VINT t, VFLOAT *p_hi_tiny,
   VFLOAT *p_lo_tiny, VFLOAT *q_hi_tiny, VFLOAT *q_lo_tiny, size_t vlen)
{
  VINT p_hi = VMVI_VX (RVV_ERFCINV_P_TINY_HI9_I, vlen);
  VINT p_lo = VMVI_VX (RVV_ERFCINV_P_TINY_LO9_I, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO8_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI8_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO7_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI7_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO6_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI6_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO5_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI5_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO4_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI4_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO3_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI3_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO2_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI2_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO1_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI1_I, t, p_hi, vlen);
  p_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_P_TINY_LO0_I, t, p_hi, p_lo, vlen);
  p_hi = PSTEP_I_HI (RVV_ERFCINV_P_TINY_HI0_I, t, p_hi, vlen);

  VINT q_hi = VMVI_VX (RVV_ERFCINV_Q_TINY_HI9_I, vlen);
  VINT q_lo = VMVI_VX (RVV_ERFCINV_Q_TINY_LO9_I, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO8_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI8_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO7_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI7_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO6_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI6_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO5_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI5_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO4_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI4_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO3_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI3_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO2_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI2_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO1_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI1_I, t, q_hi, vlen);
  q_lo = RVV_LIBMVEC_IMPL_NAME (rvv_erfcinv_update_p_lo, LMUL)
    (RVV_ERFCINV_Q_TINY_LO0_I, t, q_hi, q_lo, vlen);
  q_hi = PSTEP_I_HI (RVV_ERFCINV_Q_TINY_HI0_I, t, q_hi, vlen);

  VFLOAT a = __riscv_vfcvt_f (p_hi, vlen);
  *p_lo_tiny = __riscv_vfcvt_f (p_lo, vlen);
  *p_hi_tiny = __riscv_vfmadd
    (*p_lo_tiny, RVV_ERFCINV_TINY_DF_SCALE_D, a, vlen);
  *p_lo_tiny = __riscv_vfmadd
    (*p_lo_tiny, RVV_ERFCINV_TINY_DF_SCALE_D,
     __riscv_vfsub (a, *p_hi_tiny, vlen), vlen);

  VFLOAT b = __riscv_vfcvt_f (q_hi, vlen);
  *q_lo_tiny = __riscv_vfcvt_f (q_lo, vlen);
  *q_hi_tiny = __riscv_vfmadd
    (*q_lo_tiny, RVV_ERFCINV_TINY_DF_SCALE_D, b, vlen);
  *q_lo_tiny = __riscv_vfmadd
    (*q_lo_tiny, RVV_ERFCINV_TINY_DF_SCALE_D,
     __riscv_vfsub (b, *q_hi_tiny, vlen), vlen);
}
#endif

#undef RVV_MATH_IS_ERFCINV
#undef RVV_MATH_IS_ERFCINV_INNER
#undef RVV_MATH_IS_ERFCINV_erfcinv

/* The double-precision erfinv importer keeps its source helper operation
   graphs here so the skeleton contains only the main computation and its
   data reference.  */
#define RVV_MATH_IS_ERFINV_erfinv 1
#define RVV_MATH_IS_ERFINV_INNER(name) RVV_MATH_IS_ERFINV_##name
#define RVV_MATH_IS_ERFINV(name) RVV_MATH_IS_ERFINV_INNER (name)

#if (BIT_WIDTH == 64) && RVV_MATH_IS_ERFINV (V_MATH_FUNC)
static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_exception_handling, LMUL)
  (VFLOAT *vx, VBOOL *special_args, VFLOAT *vy_special, size_t vlen)
{
  VUINT expo_x
    = __riscv_vand (__riscv_vsrl (F_AS_U (*vx), MAN_LEN, vlen), 0x7FF, vlen);
  VBOOL x_large = __riscv_vmsgeu (expo_x, EXP_BIAS, vlen);
  VBOOL x_small
    = __riscv_vmsltu (expo_x, EXP_BIAS - RVV_ERFINV_SMALL_EXPONENT, vlen);
  *special_args = __riscv_vmor (x_large, x_small, vlen);

  if (__riscv_vcpop (*special_args, vlen) > 0)
    {
      VFLOAT abs_x = __riscv_vfsgnj (*vx, fp_posOne, vlen);
      VBOOL x_gt_1 = __riscv_vmfgt (abs_x, fp_posOne, vlen);
      VBOOL x_eq_1 = __riscv_vmfeq (abs_x, fp_posOne, vlen);
      *vx = VFMERGE_VFM (*vx, fp_sNaN, x_gt_1, vlen);

      VFLOAT tmp = VFMV_VF (fp_posZero, vlen);
      tmp = __riscv_vfsgnj (tmp, *vx, vlen);
      tmp = __riscv_vfrec7 (x_eq_1, tmp, vlen);

      *vy_special = VFADD_VV_M (*special_args, *vx, *vx, vlen);
      *vy_special = VMERGE_VVM_F (*vy_special, tmp, x_eq_1, vlen);

      tmp = __riscv_vfmul (x_small, *vx, RVV_ERFINV_RTPI_BY2_LO_D, vlen);
      tmp = __riscv_vfmacc (x_small, tmp, RVV_ERFINV_RTPI_BY2_HI_D, *vx,
                            vlen);
      *vy_special = VMERGE_VVM_F (*vy_special, tmp, x_small, vlen);

      *vx = VFMERGE_VFM (*vx, fp_posZero, *special_args, vlen);
    }
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_neg_logx, LMUL)
  (VFLOAT vx, VINT n_adjust,
   VFLOAT *y_hi, VFLOAT *y_lo, size_t vlen)
{
  VFLOAT vx_in = vx;
  VINT n
    = __riscv_vadd (__riscv_vsra (F_AS_I (vx_in), MAN_LEN - 8, vlen), 0x96,
                    vlen);
  n = __riscv_vsub (__riscv_vsra (n, 8, vlen), EXP_BIAS, vlen);
  VFLOAT scale
    = I_AS_F (__riscv_vsll (__riscv_vrsub (n, EXP_BIAS, vlen), MAN_LEN,
                            vlen));
  vx_in = __riscv_vfmul (vx_in, scale, vlen);
  n = __riscv_vsub (n, n_adjust, vlen);
  VFLOAT n_flt = __riscv_vfcvt_f (n, vlen);
  VFLOAT numer = __riscv_vfrsub (vx_in, fp_posOne, vlen);
  numer = __riscv_vfadd (numer, numer, vlen);
  VFLOAT denom = __riscv_vfadd (vx_in, fp_posOne, vlen);
  VFLOAT delta_denom
    = __riscv_vfadd (__riscv_vfrsub (denom, fp_posOne, vlen), vx_in, vlen);
  VFLOAT r_hi, r_lo, r;
  DIV_N1D2 (numer, denom, delta_denom, r_hi, r_lo, vlen);
  r = __riscv_vfadd (r_hi, r_lo, vlen);
  VFLOAT a = __riscv_vfmadd (n_flt, -RVV_LOG2_HI_D, r_hi, vlen);
  VFLOAT b = __riscv_vfmsub (n_flt, -RVV_LOG2_HI_D, a, vlen);
  b = __riscv_vfadd (b, r_hi, vlen);
  b = __riscv_vfadd (r_lo, b, vlen);
  VFLOAT rsq = __riscv_vfmul (r, r, vlen);
  VFLOAT rcube = __riscv_vfmul (rsq, r, vlen);
  VFLOAT r6 = __riscv_vfmul (rcube, rcube, vlen);
  VFLOAT poly_right
    = PSTEP (RVV_ERFINV_LOG_POLY_RIGHT0_D, rsq,
             PSTEP (RVV_ERFINV_LOG_POLY_RIGHT1_D,
                    RVV_ERFINV_LOG_POLY_RIGHT2_D, rsq,
                    vlen),
             vlen);
  VFLOAT poly_left
    = PSTEP (RVV_ERFINV_LOG_POLY_LEFT0_D, rsq,
             PSTEP (RVV_ERFINV_LOG_POLY_LEFT1_D,
                    RVV_ERFINV_LOG_POLY_LEFT2_D, rsq,
                    vlen),
             vlen);
  poly_left = __riscv_vfmacc (poly_left, r6, poly_right, vlen);
  poly_left = PSTEP (RVV_ERFINV_LOG_POLY_LEFT3_D, rsq, poly_left, vlen);
  b = __riscv_vfmacc (b, -RVV_LOG2_LO_D, n_flt, vlen);
  b = __riscv_vfmacc (b, rcube, poly_left, vlen);
  FAST2SUM (a, b, *y_hi, *y_lo, vlen);
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_erfinv_sqrtx, LMUL)
  (VFLOAT y_hi, VFLOAT y_lo, VFLOAT *w_hi, VFLOAT *w_lo, VINT *t,
   double t_scale, double t_scale_inv, size_t vlen)
{
  *w_hi = __riscv_vfsqrt (y_hi, vlen);
  *w_lo = __riscv_vfnmsub (*w_hi, *w_hi, y_hi, vlen);
  *w_lo = __riscv_vfadd (*w_lo, y_lo, vlen);
  VFLOAT recip = __riscv_vfadd (y_hi, y_hi, vlen);
  recip = __riscv_vfrec7 (recip, vlen);
  recip = __riscv_vfmul (recip, *w_hi, vlen);
  *w_lo = __riscv_vfmul (*w_lo, recip, vlen);
  VFLOAT t_lo = VFMV_VF (t_scale, vlen);
  VFLOAT t_hi = __riscv_vfdiv (t_lo, *w_hi, vlen);
  *t = __riscv_vfcvt_x (t_hi, vlen);
  t_lo = __riscv_vfnmsac (t_lo, *w_hi, t_hi, vlen);
  t_lo = __riscv_vfnmsac (t_lo, *w_lo, t_hi, vlen);
  t_lo = __riscv_vfmul (t_lo, t_hi, vlen);
  t_lo = __riscv_vfmul (t_lo, t_scale_inv, vlen);
  *t = __riscv_vadd (*t, __riscv_vfcvt_x (t_lo, vlen), vlen);
}
#endif

#undef RVV_MATH_IS_ERFINV
#undef RVV_MATH_IS_ERFINV_INNER
#undef RVV_MATH_IS_ERFINV_erfinv

/* Shared single-precision sinf helper selection.  The public clean and
   source-original skeletons keep their EH/final-result policy separate, but
   their per-LMUL source-shaped helpers belong in this common header.  */
#ifndef RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_sinf
# define RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_sinf 1
# define RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_sinf_clean 1
# define RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_rvv_libmvec_fastsinf_clean_sinf 1
# define RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_INNER(name) \
  RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_##name
# define RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG(name) \
  RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG_INNER (name)
#endif

#ifndef RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_sinf_original
# define RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_sinf_original 1
# define RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_rvv_sleef_orig_fastsinf_sinf 1
# define RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_INNER(name) \
  RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_##name
# define RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG(name) \
  RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG_INNER (name)
#endif

/* Clean sinf finite helpers.  Keep the helper graph source-shaped; the
   public sinf skeleton owns only classify/sanitize/merge and the finite call.  */
#if (BIT_WIDTH == 32) \
  && RVV_LIBMVEC_SINF_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_LIBMVEC_SINF_HELPERS_M1
#   define RVV_LIBMVEC_SINF_HELPERS_M1
#   define RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_LIBMVEC_SINF_HELPERS_M2
#   define RVV_LIBMVEC_SINF_HELPERS_M2
#   define RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_LIBMVEC_SINF_HELPERS_M4
#   define RVV_LIBMVEC_SINF_HELPERS_M4
#   define RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_LIBMVEC_SINF_HELPERS_M8
#   define RVV_LIBMVEC_SINF_HELPERS_M8
#   define RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_vrint, LMUL) (VFLOAT value, size_t vlen)
{
  VINT rounded = RVV_VFCVT_X_F_RNE (value, vlen);
  return __riscv_vfcvt_f (rounded, vlen);
}

static inline VINT
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_vilogb2k, LMUL) (VFLOAT a, size_t vlen)
{
  VUINT bits = F_AS_U_F (a);
  VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
  exponent = __riscv_vand (exponent, 0xff, vlen);
  return __riscv_vsub (U_AS_I_F (exponent), EXP_BIAS_F, vlen);
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_vldexp3, LMUL)
  (VFLOAT a, VINT q, size_t vlen)
{
  return I_AS_F_F (__riscv_vadd
                   (F_AS_I_F (a), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_vorsign, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  VUINT sign = __riscv_vand (F_AS_U_F (y), 0x80000000u, vlen);
  return U_AS_F_F (__riscv_vor (F_AS_U_F (x), sign, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_rempisubf, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_vorsign, LMUL) (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_vorsign, LMUL) (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x,
                vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x,
                vlen), vlen);
  return rem;
}

static inline void
RVV_LIBMVEC_IMPL_NAME (rvv_sinf_rempif, LMUL)
  (VFLOAT a, VFLOAT *out_hi, VFLOAT *out_lo, VINT *out_q, size_t vlen)
{
  const struct v_sleef_trig_data *trig_data =
    ptr_barrier (&__v_sleef_trig_data);

  VINT ex = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_vilogb2k, LMUL) (a, vlen);
  ex = __riscv_vsub (ex, RVV_SLEEF_SINF_REMPI_EXPONENT_BIAS, vlen);

  VBOOL scale_mask = __riscv_vmsgt
    (ex, RVV_SLEEF_SINF_REMPI_SCALE_LIMIT, vlen);
  VINT scale = VMERGE_VXM_I
    (VMVI_VX_F (0, vlen), RVV_SLEEF_SINF_REMPI_SCALE, scale_mask, vlen);
  a = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_vldexp3, LMUL) (a, scale, vlen);

  VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
  ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
  ex = __riscv_vsll (ex, RVV_SLEEF_SINF_REMPI_INDEX_SHIFT, vlen);
  VUINT ex_u = I_AS_U_F (ex);

  VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
  VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
  VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
  VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab0, x_hi, x_lo, vlen);
  VINT q;
  x_hi = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_rempisubf, LMUL) (x_hi, &q, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab1, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  VINT q_part;
  x_hi = RVV_LIBMVEC_IMPL_NAME
    (rvv_sinf_rempisubf, LMUL) (x_hi, &q_part, vlen);
  q = __riscv_vadd (q, q_part, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  y_hi = tab2;
  y_lo = tab3;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF (y_hi, y_lo, a, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo,
     VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO_PI_HI, vlen),
     VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO_PI_LO, vlen),
     x_hi, x_lo, vlen);

  VBOOL a_small = __riscv_vmflt
    (__riscv_vfabs (a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
  x_hi = VMERGE_VVM_F (x_hi, a, a_small, vlen);
  x_lo = VFMERGE_VFM (x_lo, fp_posZero_f, a_small, vlen);

  *out_hi = x_hi;
  *out_lo = x_lo;
  *out_q = q;
}

#  undef RVV_LIBMVEC_SINF_HELPERS_THIS_LMUL
# endif
#endif

/* Source-original SLEEF sinf finite helpers.  This graph intentionally keeps
   the source comparator's names and late special-value behavior independent
   from the clean public entry.  */
#if (BIT_WIDTH == 32) \
  && RVV_SLEEF_ORIG_SINF_HELPER_FUNC_FLAG (V_MATH_FUNC)
# if (LMUL == 1)
#  ifndef RVV_SLEEF_ORIG_SINF_HELPERS_M1
#   define RVV_SLEEF_ORIG_SINF_HELPERS_M1
#   define RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_SLEEF_ORIG_SINF_HELPERS_M2
#   define RVV_SLEEF_ORIG_SINF_HELPERS_M2
#   define RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_SLEEF_ORIG_SINF_HELPERS_M4
#   define RVV_SLEEF_ORIG_SINF_HELPERS_M4
#   define RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_SLEEF_ORIG_SINF_HELPERS_M8
#   define RVV_SLEEF_ORIG_SINF_HELPERS_M8
#   define RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_vrint, LMUL)
  (VFLOAT value, size_t vlen)
{
  VINT rounded = RVV_VFCVT_X_F_RNE (value, vlen);
  return __riscv_vfcvt_f (rounded, vlen);
}

static inline VINT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_vilogb2k, LMUL)
  (VFLOAT a, size_t vlen)
{
  VUINT bits = F_AS_U (a);
  VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
  exponent = __riscv_vand (exponent, 0xff, vlen);
  return __riscv_vsub (U_AS_I (exponent), EXP_BIAS_F, vlen);
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_vldexp3, LMUL)
  (VFLOAT a, VINT q, size_t vlen)
{
  return I_AS_F (__riscv_vadd
                 (F_AS_I (a), __riscv_vsll (q, MAN_LEN_F, vlen), vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_vorsign, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  VUINT sign = __riscv_vand (F_AS_U (y), 0x80000000u, vlen);
  return U_AS_F (__riscv_vor (F_AS_U (x), sign, vlen));
}

static inline VFLOAT
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_rempisubf, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_vorsign, LMUL) (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_vorsign, LMUL) (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x,
              vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x,
                vlen), vlen);
  return rem;
}

static inline void
RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_sinf_rempif, LMUL)
  (VFLOAT a, VFLOAT *out_hi, VFLOAT *out_lo, VINT *out_q, size_t vlen)
{
  const struct v_sleef_trig_data *trig_data =
    ptr_barrier (&__v_sleef_trig_data);

  VINT ex = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_vilogb2k, LMUL) (a, vlen);
  ex = __riscv_vsub (ex, RVV_SLEEF_SINF_REMPI_EXPONENT_BIAS, vlen);

  VBOOL scale_mask = __riscv_vmsgt
    (ex, RVV_SLEEF_SINF_REMPI_SCALE_LIMIT, vlen);
  VINT scale = VMERGE_VXM_I
    (VMVI_VX (0, vlen), RVV_SLEEF_SINF_REMPI_SCALE, scale_mask, vlen);
  a = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_vldexp3, LMUL) (a, scale, vlen);

  VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
  ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
  ex = __riscv_vsll (ex, RVV_SLEEF_SINF_REMPI_INDEX_SHIFT, vlen);
  VUINT ex_u = I_AS_U (ex);

  VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
  VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
  VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
  VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab0, x_hi, x_lo, vlen);
  VINT q;
  x_hi = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_rempisubf, LMUL) (x_hi, &q, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab1, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  VINT q_part;
  x_hi = RVV_SLEEF_ORIG_MANGLE
    (rvv_sleef_orig_sinf_rempisubf, LMUL) (x_hi, &q_part, vlen);
  q = __riscv_vadd (q, q_part, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  y_hi = tab2;
  y_lo = tab3;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (y_hi, y_lo, a, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo,
     VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO_PI_HI, vlen),
     VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO_PI_LO, vlen),
     x_hi, x_lo, vlen);

  VBOOL a_small = __riscv_vmflt
    (__riscv_vfabs (a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
  x_hi = VMERGE_VVM_F (x_hi, a, a_small, vlen);
  x_lo = VFMERGE_VFM (x_lo, fp_posZero_f, a_small, vlen);

  *out_hi = x_hi;
  *out_lo = x_lo;
  *out_q = q;
}

#  undef RVV_SLEEF_ORIG_SINF_HELPERS_THIS_LMUL
# endif
#endif
#ifndef RVV_MATH_IS_LGAMMA_lgamma
# define RVV_MATH_IS_LGAMMA_lgamma 1
# define RVV_MATH_IS_LGAMMA_INNER(name) RVV_MATH_IS_LGAMMA_##name
# define RVV_MATH_IS_LGAMMA(name) RVV_MATH_IS_LGAMMA_INNER (name)
#endif
#if RVV_MATH_IS_LGAMMA (V_MATH_FUNC)
#ifndef V_D_LGAMMA_SKELETON_ONCE
#define V_D_LGAMMA_SKELETON_ONCE

#define RVV_LGAMMA_P_LEFT0_I ((int64_t) 0x717ad571ef5dc61)
#define RVV_LGAMMA_P_LEFT1_I ((int64_t) 0x13b685202b6fd6a4)
#define RVV_LGAMMA_P_LEFT2_I ((int64_t) 0x16296a1296488970)
#define RVV_LGAMMA_P_LEFT3_I ((int64_t) 0x19ca4fa8bc759bde)
#define RVV_LGAMMA_P_LEFT4_I ((int64_t) 0x1079dccb79c3089c)
#define RVV_LGAMMA_P_LEFT5_I ((int64_t) 0x16847f57936dc8fb)
#define RVV_LGAMMA_P_LEFT6_I ((int64_t) 0x1d546ec89ba5d14e)
#define RVV_LGAMMA_P_LEFT7_I ((int64_t) 0x1ad0cdf5663cfacf)
#define RVV_LGAMMA_P_LEFT8_I ((int64_t) 0x1f1f571e9999b6c7)
#define RVV_LGAMMA_Q_LEFT0_I ((int64_t) 0x3877618277fdb07)
#define RVV_LGAMMA_Q_LEFT1_I ((int64_t) 0x15e97ae7797617d9)
#define RVV_LGAMMA_Q_LEFT2_I ((int64_t) 0x1c1a630c4311e499)
#define RVV_LGAMMA_Q_LEFT3_I ((int64_t) 0x133731ff844280b6)
#define RVV_LGAMMA_Q_LEFT4_I ((int64_t) 0x1e248b8ebb59488b)
#define RVV_LGAMMA_Q_LEFT5_I ((int64_t) 0x1b1f236d04a448b5)
#define RVV_LGAMMA_Q_LEFT6_I ((int64_t) 0x1a6052955f6252ac)
#define RVV_LGAMMA_Q_LEFT7_I ((int64_t) 0x17e477b664d95b52)
#define RVV_LGAMMA_Q_LEFT8_I ((int64_t) 0x1b97bd09b9b48410)
#define RVV_LGAMMA_P_LOG0_I ((int64_t) 0x5555555555555090)
#define RVV_LGAMMA_P_LOG1_I ((int64_t) 0x666666666686863a)
#define RVV_LGAMMA_P_LOG2_I ((int64_t) 0x49249248fc99ba4b)
#define RVV_LGAMMA_P_LOG3_I ((int64_t) 0x71c71ca402e164fa)
#define RVV_LGAMMA_P_LOG4_I ((int64_t) 0x5d1733e3ae94dde0)
#define RVV_LGAMMA_P_LOG5_I ((int64_t) 0x4ec8b69784234032)
#define RVV_LGAMMA_P_LOG6_I ((int64_t) 0x43cc44a056dc3c93)
#define RVV_LGAMMA_P_LOG7_I ((int64_t) 0x4432439bb76e7d74)
#define RVV_LGAMMA_P_LS_CORR0_I ((int64_t) 0x13eb19ce38760e4)
#define RVV_LGAMMA_P_LS_CORR1_I ((int64_t) 0x54ebdd91a33a236)
#define RVV_LGAMMA_P_LS_CORR2_I ((int64_t) 0xf5302c2f3171924)
#define RVV_LGAMMA_P_LS_CORR3_I ((int64_t) 0x17e6ca6c67d42c45)
#define RVV_LGAMMA_P_LS_CORR4_I ((int64_t) 0x18e683b7eb793968)
#define RVV_LGAMMA_P_LS_CORR5_I ((int64_t) 0xe6a7d68df697b37)
#define RVV_LGAMMA_P_LS_CORR6_I ((int64_t) 0x48f07444527e095)
#define RVV_LGAMMA_P_LS_CORR7_I ((int64_t) 0x5ac3ca10d36d7d)
#define RVV_LGAMMA_P_LS_CORR8_I ((int64_t) -0x115718edb07d53)
#define RVV_LGAMMA_Q_LS_CORR0_I ((int64_t) 0x2f8b7a297052f62)
#define RVV_LGAMMA_Q_LS_CORR1_I ((int64_t) 0xc13fa37f8190cf5)
#define RVV_LGAMMA_Q_LS_CORR2_I ((int64_t) 0x222d203fd991122d)
#define RVV_LGAMMA_Q_LS_CORR3_I ((int64_t) 0x32462b6d38e0bfd3)
#define RVV_LGAMMA_Q_LS_CORR4_I ((int64_t) 0x31844bff55d6561a)
#define RVV_LGAMMA_Q_LS_CORR5_I ((int64_t) 0x18c83406788ab40e)
#define RVV_LGAMMA_Q_LS_CORR6_I ((int64_t) 0x643329f595fac69)
#define RVV_LGAMMA_Q_LS_CORR7_I ((int64_t) -0x21b0b1bff373cd)
#define RVV_LGAMMA_Q_LS_CORR8_I ((int64_t) -0xc9c05b696db07)
#define RVV_LGAMMA_P_LOGSIN0_I ((int64_t) -0x34a34cc4a60fa863)
#define RVV_LGAMMA_P_LOGSIN1_I ((int64_t) -0x1151322ac7d51d2e)
#define RVV_LGAMMA_P_LOGSIN2_I ((int64_t) -0xada0658820c4e34)
#define RVV_LGAMMA_P_LOGSIN3_I ((int64_t) -0x80859b50a7b1918)
#define RVV_LGAMMA_P_LOGSIN4_I ((int64_t) -0x66807a019daf246)
#define RVV_LGAMMA_P_LOGSIN5_I ((int64_t) -0x555a97e7d8482c8)
#define RVV_LGAMMA_P_LOGSIN6_I ((int64_t) -0x4927ceefdc18f62)
#define RVV_LGAMMA_P_LOGSIN7_I ((int64_t) -0x3fe5862d4e702a2)
#define RVV_LGAMMA_P_LOGSIN8_I ((int64_t) -0x39da522c5099734)
#define RVV_LGAMMA_P_LOGSIN9_I ((int64_t) -0x2cbb6825e3efaad)
#define RVV_LGAMMA_P_LOGSIN10_I ((int64_t) -0x4df815d2f21e674)
#define RVV_LGAMMA_P_LOGSIN11_I ((int64_t) 0x41cf7e791cb446c)
#define RVV_LGAMMA_P_LOGSIN12_I ((int64_t) -0x126ea0159b1a7052)
#define RVV_LGAMMA_P_LOGSIN13_I ((int64_t) 0x155103f2634da2c6)
#define RVV_LGAMMA_P_LOGSIN14_I ((int64_t) -0x13e497482ec6dff4)
#define RVV_LGAMMA_P_EXP3_I ((int64_t) 0x1111111110ec10d)
#define RVV_LGAMMA_P_EXP4_I ((int64_t) 0x2d82d82d87a9b5)
#define RVV_LGAMMA_P_EXP5_I ((int64_t) 0x6806806ce6d6f)
#define RVV_LGAMMA_P_EXP6_I ((int64_t) 0xd00d00841fcf)
#define RVV_LGAMMA_P_EXP7_I ((int64_t) 0x171ddefda54b)
#define RVV_LGAMMA_P_EXP8_I ((int64_t) 0x24fcc01d627)
#define RVV_LGAMMA_P_EXP9_I ((int64_t) 0x35ed8bbd24)
#define RVV_LGAMMA_P_EXP10_I ((int64_t) 0x477745b6c)

#define RVV_LGAMMA_NEG_BND_TO_2_X0_0_D 0x1.3a7fc9600f86cp+1
#define RVV_LGAMMA_NEG_BND_TO_2_X0_1_D 0x1.55f64f98af8dp-55
#define RVV_LGAMMA_NEG_BND_TO_2_X0_2_D 0x1.c4b0cd201366ap-110
#define RVV_LGAMMA_NEG_3_TO_BND_X0_0_D 0x1.5fb410a1bd901p+1
#define RVV_LGAMMA_NEG_3_TO_BND_X0_1_D -0x1.a19a96d2e6f85p-54
#define RVV_LGAMMA_NEG_3_TO_BND_X0_2_D -0x1.140b4ff4b7d6p-108

/* Gamma(x) ~ (x-1)(x-2) * P(t)/Q(t), t = x - (3/2 - 1/8) x in [1-1/4, 2+1/4]
   Coefficients P in ascending order */
#define P_left_0 RVV_LGAMMA_P_LEFT0_I
#define P_left_1 RVV_LGAMMA_P_LEFT1_I
#define P_left_2 RVV_LGAMMA_P_LEFT2_I
#define P_left_3 RVV_LGAMMA_P_LEFT3_I
#define P_left_4 RVV_LGAMMA_P_LEFT4_I
#define P_left_5 RVV_LGAMMA_P_LEFT5_I
#define P_left_6 RVV_LGAMMA_P_LEFT6_I
#define P_left_7 RVV_LGAMMA_P_LEFT7_I
#define P_left_8 RVV_LGAMMA_P_LEFT8_I

#define Q_left_0 RVV_LGAMMA_Q_LEFT0_I
#define Q_left_1 RVV_LGAMMA_Q_LEFT1_I
#define Q_left_2 RVV_LGAMMA_Q_LEFT2_I
#define Q_left_3 RVV_LGAMMA_Q_LEFT3_I
#define Q_left_4 RVV_LGAMMA_Q_LEFT4_I
#define Q_left_5 RVV_LGAMMA_Q_LEFT5_I
#define Q_left_6 RVV_LGAMMA_Q_LEFT6_I
#define Q_left_7 RVV_LGAMMA_Q_LEFT7_I
#define Q_left_8 RVV_LGAMMA_Q_LEFT8_I

/* Approximate log(x) by w + w^3 poly(w^2)
   w = 2(x-1)/(x+1), x roughly in [1/rt(2), rt(2)] */
#define P_log_0 RVV_LGAMMA_P_LOG0_I
#define P_log_1 RVV_LGAMMA_P_LOG1_I
#define P_log_2 RVV_LGAMMA_P_LOG2_I
#define P_log_3 RVV_LGAMMA_P_LOG3_I
#define P_log_4 RVV_LGAMMA_P_LOG4_I
#define P_log_5 RVV_LGAMMA_P_LOG5_I
#define P_log_6 RVV_LGAMMA_P_LOG6_I
#define P_log_7 RVV_LGAMMA_P_LOG7_I

#define LOG2_HI RVV_LOG2_HI_D
#define LOG2_LO RVV_LOG2_LO_D

/* lgamma(x) overflows for |x| >= ~2.56e305 (loggamma > DBL_MAX). */
#define LGAMMA_OVERFLOW_X RVV_LGAMMA_OVERFLOW_X_D

/* Correction to log_stirling formula
   logGamma(x) - log_sterling(x) ~ P(t)/Q(t), t = 1/x
   x in [2.25 x_max], t in (0, 1/2.25) */
#define P_LS_corr_0 RVV_LGAMMA_P_LS_CORR0_I
#define P_LS_corr_1 RVV_LGAMMA_P_LS_CORR1_I
#define P_LS_corr_2 RVV_LGAMMA_P_LS_CORR2_I
#define P_LS_corr_3 RVV_LGAMMA_P_LS_CORR3_I
#define P_LS_corr_4 RVV_LGAMMA_P_LS_CORR4_I
#define P_LS_corr_5 RVV_LGAMMA_P_LS_CORR5_I
#define P_LS_corr_6 RVV_LGAMMA_P_LS_CORR6_I
#define P_LS_corr_7 RVV_LGAMMA_P_LS_CORR7_I
#define P_LS_corr_8 RVV_LGAMMA_P_LS_CORR8_I

#define Q_LS_corr_0 RVV_LGAMMA_Q_LS_CORR0_I
#define Q_LS_corr_1 RVV_LGAMMA_Q_LS_CORR1_I
#define Q_LS_corr_2 RVV_LGAMMA_Q_LS_CORR2_I
#define Q_LS_corr_3 RVV_LGAMMA_Q_LS_CORR3_I
#define Q_LS_corr_4 RVV_LGAMMA_Q_LS_CORR4_I
#define Q_LS_corr_5 RVV_LGAMMA_Q_LS_CORR5_I
#define Q_LS_corr_6 RVV_LGAMMA_Q_LS_CORR6_I
#define Q_LS_corr_7 RVV_LGAMMA_Q_LS_CORR7_I
#define Q_LS_corr_8 RVV_LGAMMA_Q_LS_CORR8_I

/* Approximate log(sin(pi x)/(pi x)) as x^2 poly(x^2) */
#define P_logsin_0 RVV_LGAMMA_P_LOGSIN0_I
#define P_logsin_1 RVV_LGAMMA_P_LOGSIN1_I
#define P_logsin_2 RVV_LGAMMA_P_LOGSIN2_I
#define P_logsin_3 RVV_LGAMMA_P_LOGSIN3_I
#define P_logsin_4 RVV_LGAMMA_P_LOGSIN4_I
#define P_logsin_5 RVV_LGAMMA_P_LOGSIN5_I
#define P_logsin_6 RVV_LGAMMA_P_LOGSIN6_I
#define P_logsin_7 RVV_LGAMMA_P_LOGSIN7_I
#define P_logsin_8 RVV_LGAMMA_P_LOGSIN8_I
#define P_logsin_9 RVV_LGAMMA_P_LOGSIN9_I
#define P_logsin_10 RVV_LGAMMA_P_LOGSIN10_I
#define P_logsin_11 RVV_LGAMMA_P_LOGSIN11_I
#define P_logsin_12 RVV_LGAMMA_P_LOGSIN12_I
#define P_logsin_13 RVV_LGAMMA_P_LOGSIN13_I
#define P_logsin_14 RVV_LGAMMA_P_LOGSIN14_I
#define P_exp_3 RVV_LGAMMA_P_EXP3_I
#define P_exp_4 RVV_LGAMMA_P_EXP4_I
#define P_exp_5 RVV_LGAMMA_P_EXP5_I
#define P_exp_6 RVV_LGAMMA_P_EXP6_I
#define P_exp_7 RVV_LGAMMA_P_EXP7_I
#define P_exp_8 RVV_LGAMMA_P_EXP8_I
#define P_exp_9 RVV_LGAMMA_P_EXP9_I
#define P_exp_10 RVV_LGAMMA_P_EXP10_I

#define LGAMMA_NEG_BND_TO_2_SC RVV_LGAMMA_NEG_BND_TO_2_SC_D
#define LGAMMA_NEG_3_TO_BND_SC RVV_LGAMMA_NEG_3_TO_BND_SC_D

/* The finite lgamma kernel is evaluated only for lanes selected by the
   importer-owned `lgamma_mask`.  Masked RVV arithmetic has agnostic inactive
   destinations on this toolchain, so every masked result is immediately
   sanitized to +0 before it can feed an unmasked integer conversion or the
   next reduction step.  Active-lane order and operation selection remain the
   source order.  */
#define LGAMMA_MASKED_RESULT(mask, value, vlen)                              \
  VMERGE_VVM_F (VFMV_VF (fp_posZero, (vlen)), (value), (mask), (vlen))
#define LGAMMA_VFADD_VV_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFADD_VV_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFADD_VF_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFADD_VF_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFSUB_VV_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFSUB_VV_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFSUB_VF_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFSUB_VF_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFMUL_VV_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFMUL_VV_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFMUL_VF_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFMUL_VF_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFMSUB_VV_M(mask, x, y, z, vlen)                               \
  LGAMMA_MASKED_RESULT ((mask), VFMSUB_VV_M ((mask), (x), (y), (z), (vlen)),  \
                        (vlen))
#define LGAMMA_VFMSUB_VF_M(mask, x, y, z, vlen)                               \
  LGAMMA_MASKED_RESULT ((mask), VFMSUB_VF_M ((mask), (x), (y), (z), (vlen)),  \
                        (vlen))
#define LGAMMA_VFMACC_VV_M(mask, x, y, z, vlen)                               \
  LGAMMA_MASKED_RESULT ((mask), VFMACC_VV_M ((mask), (x), (y), (z), (vlen)),  \
                        (vlen))
#define LGAMMA_VFMACC_VF_M(mask, x, y, z, vlen)                               \
  LGAMMA_MASKED_RESULT ((mask), VFMACC_VF_M ((mask), (x), (y), (z), (vlen)),  \
                        (vlen))
#define LGAMMA_VFDIV_VV_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFDIV_VV_M ((mask), (x), (y), (vlen)),        \
                        (vlen))
#define LGAMMA_VFNMSUB_VV_M(mask, x, y, z, vlen)                              \
  LGAMMA_MASKED_RESULT ((mask), VFNMSUB_VV_M ((mask), (x), (y), (z), (vlen)), \
                        (vlen))
#define LGAMMA_VFNMSAC_VV_M(mask, x, y, z, vlen)                              \
  LGAMMA_MASKED_RESULT ((mask), VFNMSAC_VV_M ((mask), (x), (y), (z), (vlen)), \
                        (vlen))
#define LGAMMA_VFRDIV_VF_M(mask, x, y, vlen)                                  \
  LGAMMA_MASKED_RESULT ((mask), VFRDIV_VF_M ((mask), (x), (y), (vlen)),       \
                        (vlen))
#define LGAMMA_VFREC7_M(mask, x, vlen)                                       \
  LGAMMA_MASKED_RESULT ((mask), VFREC7_M ((mask), (x), (vlen)), (vlen))

#define LGAMMA_VFADD(x, y, vlen)                                             \
  LGAMMA_VFADD_VV_M (lgamma_mask, (x), (y), (vlen))
#define LGAMMA_VFSUB(x, y, vlen)                                             \
  LGAMMA_VFSUB_VV_M (lgamma_mask, (x), (y), (vlen))
#define LGAMMA_VFMUL(x, y, vlen)                                             \
  LGAMMA_VFMUL_VV_M (lgamma_mask, (x), (y), (vlen))
#define LGAMMA_VFMUL_VF(x, y, vlen)                                          \
  LGAMMA_VFMUL_VF_M (lgamma_mask, (x), (y), (vlen))
#define LGAMMA_VFMSUB(x, y, z, vlen)                                         \
  LGAMMA_VFMSUB_VV_M (lgamma_mask, (x), (y), (z), (vlen))
#define LGAMMA_VFNMSUB(x, y, z, vlen)                                        \
  LGAMMA_VFNMSUB_VV_M (lgamma_mask, (x), (y), (z), (vlen))

#define LGAMMA_PROD_X1Y1_M(mask, x, y, prod_hi, prod_lo, vlen)               \
  do                                                                          \
    {                                                                         \
      (prod_hi) = LGAMMA_VFMUL_VV_M ((mask), (x), (y), (vlen));              \
      (prod_lo) = LGAMMA_VFMSUB_VV_M ((mask), (x), (y), (prod_hi), (vlen));   \
    }                                                                         \
  while (0)

#define LGAMMA_PROD_X1Y2_M(mask, x, y_hi, y_lo, prod_hi, prod_lo, vlen)       \
  do                                                                          \
    {                                                                         \
      (prod_hi) = LGAMMA_VFMUL_VV_M ((mask), (x), (y_hi), (vlen));            \
      (prod_lo) = LGAMMA_VFMSUB_VV_M ((mask), (x), (y_hi), (prod_hi),        \
                                      (vlen));                                \
      (prod_lo) = LGAMMA_VFMACC_VV_M ((mask), (prod_lo), (x), (y_lo),        \
                                      (vlen));                                \
    }                                                                         \
  while (0)

#define LGAMMA_PROD_X2Y2_M(mask, x_hi, x_lo, y_hi, y_lo, prod_hi, prod_lo,   \
                           vlen)                                              \
  do                                                                          \
    {                                                                         \
      (prod_hi) = LGAMMA_VFMUL_VV_M ((mask), (x_hi), (y_hi), (vlen));         \
      (prod_lo) = LGAMMA_VFMSUB_VV_M ((mask), (x_hi), (y_hi), (prod_hi),     \
                                      (vlen));                                \
      (prod_lo) = LGAMMA_VFMACC_VV_M ((mask), (prod_lo), (x_hi), (y_lo),     \
                                      (vlen));                                \
      (prod_lo) = LGAMMA_VFMACC_VV_M ((mask), (prod_lo), (x_lo), (y_hi),     \
                                      (vlen));                                \
    }                                                                         \
  while (0)

#define LGAMMA_KNUTH2SUM_M(mask, x, y, sum_hi, sum_lo, vlen)                 \
  do                                                                          \
    {                                                                         \
      (sum_hi) = LGAMMA_VFADD_VV_M ((mask), (x), (y), (vlen));                \
      VFLOAT _x_hat = LGAMMA_VFSUB_VV_M ((mask), (sum_hi), (y), (vlen));      \
      VFLOAT _y_hat = LGAMMA_VFSUB_VV_M ((mask), (sum_hi), _x_hat, (vlen));   \
      (sum_lo) = LGAMMA_VFADD_VV_M                                          \
        ((mask), LGAMMA_VFSUB_VV_M ((mask), (x), _x_hat, (vlen)),             \
         LGAMMA_VFSUB_VV_M ((mask), (y), _y_hat, (vlen)), (vlen));              \
    }                                                                         \
  while (0)

#define LGAMMA_FAST2SUM_M(mask, x, y, sum_hi, sum_lo, vlen)                  \
  do                                                                          \
    {                                                                         \
      (sum_hi) = LGAMMA_VFADD_VV_M ((mask), (x), (y), (vlen));                \
      (sum_lo) = LGAMMA_VFADD_VV_M                                         \
        ((mask), LGAMMA_VFSUB_VV_M ((mask), (x), (sum_hi), (vlen)), (y),     \
         (vlen));                                                             \
    }                                                                         \
  while (0)

#define LGAMMA_DIV2_N2D2_M(mask, numer, delta_n, denom, delta_d, q_hi,       \
                           q_lo, vlen)                                       \
  do                                                                          \
    {                                                                         \
      VFLOAT _q;                                                             \
      (q_hi) = LGAMMA_VFDIV_VV_M ((mask), (numer), (denom), (vlen));         \
      _q = LGAMMA_VFNMSUB_VV_M ((mask), (q_hi), (denom), (numer), (vlen));    \
      _q = LGAMMA_VFNMSAC_VV_M ((mask), _q, (q_hi), (delta_d), (vlen));       \
      _q = LGAMMA_VFADD_VV_M ((mask), _q, (delta_n), (vlen));                 \
      (q_lo) = LGAMMA_VFMUL_VV_M                                            \
        ((mask), _q, LGAMMA_VFREC7_M ((mask), (denom), (vlen)), (vlen));      \
    }                                                                         \
  while (0)

#define LGAMMA_ACC_DIV2_N2D2_M(mask, numer, delta_n, denom, delta_d, q_hi,   \
                               q_lo, vlen)                                    \
  do                                                                          \
    {                                                                         \
      VFLOAT _recip;                                                         \
      VFLOAT _q;                                                             \
      _recip = LGAMMA_VFRDIV_VF_M ((mask), (denom), 0x1.0p0, (vlen));        \
      (q_hi) = LGAMMA_VFMUL_VV_M ((mask), (numer), _recip, (vlen));           \
      _q = LGAMMA_VFNMSUB_VV_M ((mask), (q_hi), (denom), (numer), (vlen));    \
      _q = LGAMMA_VFNMSAC_VV_M ((mask), _q, (q_hi), (delta_d), (vlen));       \
      _q = LGAMMA_VFADD_VV_M ((mask), _q, (delta_n), (vlen));                 \
      (q_lo) = LGAMMA_VFMUL_VV_M ((mask), _q, _recip, (vlen));                \
    }                                                                         \
  while (0)

/* This macro computes loggamma(x) for 0 < x <= 2.25
   It uses a rational approximation for x in [0.75, 2.25]
   For x < 0.75, it uses the relation gamma(x) = gamma(x+1)/x */
#define LGAMMA_LE_225(x, logx_hi, logx_lo, y_hi, y_lo, vlen)                   \
  do                                                                           \
    {                                                                          \
      VBOOL x_lt_75 = __riscv_vmflt ((x), 0x1.8p-1, (vlen));                   \
      VFLOAT c = VFMV_VF (0x1.6p63, (vlen));                                   \
      c = __riscv_vfmerge (c, 0x1.8p61, x_lt_75, (vlen));                      \
      VFLOAT rt1 = VFMV_VF (0x1.0p0, (vlen));                                  \
      VFLOAT rt2 = VFMV_VF (0x1.0p1, (vlen));                                  \
      rt1 = __riscv_vfmerge (rt1, fp_posZero, x_lt_75, (vlen));                \
      rt2 = __riscv_vfmerge (rt2, fp_posOne, x_lt_75, (vlen));                 \
      VFLOAT t_hi = LGAMMA_VFMSUB_VF_M                                      \
        (lgamma_mask, (x), 0x1.0p63, c, (vlen));                              \
      VFLOAT t_lo = LGAMMA_VFADD (t_hi, c, (vlen));                           \
      t_lo = LGAMMA_VFMSUB_VF_M (lgamma_mask, (x), 0x1.0p63, t_lo,           \
                                 (vlen));                                     \
      VFLOAT fact1 = LGAMMA_VFSUB ((x), rt1, (vlen));                         \
      VINT T = __riscv_vfcvt_x (t_hi, (vlen));                                 \
      T = __riscv_vadd (T, __riscv_vfcvt_x (t_lo, (vlen)), (vlen));            \
      VFLOAT fact2_hi = LGAMMA_VFSUB ((x), rt2, (vlen));                      \
      VFLOAT fact2_lo = LGAMMA_VFADD (fact2_hi, rt2, (vlen));                 \
      fact2_lo = LGAMMA_VFSUB ((x), fact2_lo, (vlen));                        \
      VFLOAT fact_hi, fact_lo;                                                 \
      LGAMMA_PROD_X1Y2_M (lgamma_mask, fact1, fact2_hi, fact2_lo, fact_hi,    \
                          fact_lo, (vlen));                                   \
      VINT P = PSTEP_I_SRA (P_left_7, P_left_8, 9, T, (vlen));                 \
      P = PSTEP_I_SRA (P_left_6, T, 5, P, (vlen));                             \
      P = PSTEP_I_SRA (P_left_5, T, 4, P, (vlen));                             \
      P = PSTEP_I_SRA (P_left_4, T, 3, P, (vlen));                             \
      P = PSTEP_I_SRA (P_left_3, T, 1, P, (vlen));                             \
      P = PSTEP_I_SRA (P_left_2, T, 1, P, (vlen));                             \
      P = PSTEP_I (P_left_1, T, P, (vlen));                                    \
      P = PSTEP_I (P_left_0, T, P, (vlen));                                    \
      VINT Q = PSTEP_I_SRA (Q_left_7, Q_left_8, 6, T, (vlen));                 \
      Q = PSTEP_I_SRA (Q_left_6, T, 4, Q, (vlen));                             \
      Q = PSTEP_I_SRA (Q_left_5, T, 3, Q, (vlen));                             \
      Q = PSTEP_I_SRA (Q_left_4, T, 2, Q, (vlen));                             \
      Q = PSTEP_I_SRA (Q_left_3, T, 2, Q, (vlen));                             \
      Q = PSTEP_I (Q_left_1, T, PSTEP_I (Q_left_2, T, Q, (vlen)), (vlen));     \
      Q = PSTEP_I_SRA (Q_left_0, T, 1, Q, (vlen));                             \
      /* P is in Q69 and Q is in Q67 */                                        \
      P = VMERGE_VVM_I (VMVI_VX (0, (vlen)), P, lgamma_mask, (vlen));           \
      Q = VMERGE_VVM_I (VMVI_VX (0, (vlen)), Q, lgamma_mask, (vlen));           \
      VFLOAT p_hi = __riscv_vfcvt_f (P, (vlen));                               \
      VFLOAT p_lo = __riscv_vfcvt_f (                                          \
          __riscv_vsub (P, __riscv_vfcvt_x (p_hi, (vlen)), (vlen)), (vlen));   \
      VFLOAT q_hi = __riscv_vfcvt_f (Q, (vlen));                               \
      VFLOAT q_lo = __riscv_vfcvt_f (                                          \
          __riscv_vsub (Q, __riscv_vfcvt_x (q_hi, (vlen)), (vlen)), (vlen));   \
      VFLOAT z_hi, z_lo;                                                       \
      LGAMMA_DIV2_N2D2_M (lgamma_mask, p_hi, p_lo, q_hi, q_lo, z_hi, z_lo,   \
                          (vlen));                                           \
      z_hi = LGAMMA_VFMUL_VF (z_hi, 0x1.0p-2, (vlen));                        \
      z_lo = LGAMMA_VFMUL_VF (z_lo, 0x1.0p-2, (vlen));                        \
      LGAMMA_PROD_X2Y2_M (lgamma_mask, z_hi, z_lo, fact_hi, fact_lo,           \
                          (y_hi), (y_lo), (vlen));                            \
      /* if original input is in (0, 3/4), need to add -log(x) */             \
      VFLOAT A = VFMV_VF (fp_posZero, (vlen));                                \
      VFLOAT a = VFMV_VF (fp_posZero, (vlen));                                \
      A = __riscv_vmerge (A, (logx_hi), x_lt_75, (vlen));                      \
      a = __riscv_vmerge (a, (logx_lo), x_lt_75, (vlen));                      \
      /* y_hi + y_lo - (A + a), A is either 0 or dominates */                 \
      z_hi = LGAMMA_VFSUB ((y_hi), A, (vlen));                                 \
      z_lo = LGAMMA_VFADD (z_hi, A, (vlen));                                   \
      z_lo = LGAMMA_VFSUB ((y_hi), z_lo, (vlen));                              \
      (y_lo) = LGAMMA_VFADD ((y_lo), z_lo, (vlen));                            \
      (y_lo) = LGAMMA_VFSUB ((y_lo), a, (vlen));                               \
      (y_hi) = z_hi;                                                           \
    }                                                                          \
  while (0)

/* Compute log(x/e) or log(x) to 2^(-65) absolute accuracy.
   log(x) - c, c is 1 or 0; x > 0. */
#define LGAMMA_LOG(x_hi, x_lo, c, y_hi, y_lo, vlen)                            \
  do                                                                           \
    {                                                                          \
      /* need x_hi, x_lo as input */                                           \
      VFLOAT x_in_hi = (x_hi);                                                 \
      VFLOAT x_in_lo = (x_lo);                                                 \
      VINT n_adjust = VMVI_VX (0, (vlen));                                    \
      VBOOL x_tiny = __riscv_vmflt (x_in_hi, 0x1.0p-1020, (vlen));             \
      if (__riscv_vcpop (x_tiny, (vlen)) > 0)                                  \
        {                                                                      \
          VBOOL x_tiny_active = __riscv_vmand                         \
            (lgamma_mask, x_tiny, (vlen));                                     \
          VFLOAT x_adjust = LGAMMA_VFMUL_VF_M                                  \
            (x_tiny_active, x_in_hi, 0x1.0p60, (vlen));                         \
          x_in_hi = __riscv_vmerge (x_in_hi, x_adjust, x_tiny, (vlen));        \
          x_adjust = LGAMMA_VFMUL_VF_M                                         \
            (x_tiny_active, x_in_lo, 0x1.0p60, (vlen));                        \
          x_in_lo = __riscv_vmerge (x_in_lo, x_adjust, x_tiny, (vlen));        \
          n_adjust = __riscv_vmerge (n_adjust, 60, x_tiny, (vlen));            \
        }                                                                      \
      VINT n = __riscv_vadd (                                                  \
          __riscv_vsra (F_AS_I (x_in_hi), MAN_LEN - 8, (vlen)), 0x96, vlen);   \
      n = __riscv_vsub (__riscv_vsra (n, 8, vlen), EXP_BIAS, vlen);            \
      VFLOAT scale = I_AS_F (                                                  \
          __riscv_vsll (__riscv_vrsub (n, EXP_BIAS, (vlen)), MAN_LEN, (vlen)));\
      x_in_hi = LGAMMA_VFMUL (x_in_hi, scale, (vlen));                         \
      x_in_lo = LGAMMA_VFMUL (x_in_lo, scale, (vlen));                         \
      n = __riscv_vsub (n, n_adjust, (vlen));                                  \
      /* x is scaled, and log(x) is 2 atanh(w/2); w = 2(x-1)/(x+1) */         \
                                                                             \
      VFLOAT numer, denom, denom_delta;                                       \
      numer = LGAMMA_VFSUB_VF_M (lgamma_mask, x_in_hi, fp_posOne, (vlen));     \
      denom = LGAMMA_VFADD_VF_M (lgamma_mask, x_in_hi, fp_posOne, (vlen));     \
      denom_delta                                                              \
          = LGAMMA_VFADD_VV_M (lgamma_mask,                                   \
              __riscv_vfrsub (denom, fp_posOne, (vlen)), x_in_hi, (vlen));     \
      denom_delta = LGAMMA_VFADD (denom_delta, x_in_lo, (vlen));               \
      VFLOAT w_hi, w_lo;                                                      \
      LGAMMA_ACC_DIV2_N2D2_M (lgamma_mask, numer, x_in_lo, denom,             \
                              denom_delta, w_hi, w_lo, vlen);                  \
      /* w_hi + w_lo is at this point (x-1)/(x+1) */                          \
      /* Next get 2(x-1)/(x+1) in Q64 fixed point */                          \
      VINT W = __riscv_vfcvt_x                                         \
        (LGAMMA_VFMUL_VF (w_hi, 0x1.0p65, (vlen)), (vlen));                    \
      W = __riscv_vadd (                                                       \
          W, __riscv_vfcvt_x (LGAMMA_VFMUL_VF (w_lo, 0x1.0p65, (vlen)),       \
                              (vlen)),                                       \
          (vlen));                                                            \
      /* W is in Q64 because W is 2(x-1)/(x+1) */                             \
                                                                             \
      n = VMERGE_VVM_I (VMVI_VX (0, (vlen)), n, lgamma_mask, (vlen));           \
      VFLOAT n_flt = __riscv_vfcvt_f (n, (vlen));                              \
      VINT W2 = __riscv_vsmul (W, W, 1, (vlen)); /* Q65 */                     \
                                                                             \
      VINT P_right, P_left, W8;                                               \
      P_right = PSTEP_I_SRA (P_log_6, P_log_7, 4, W2, (vlen));                 \
      P_right = PSTEP_I_SRA (P_log_5, W2, 4, P_right, (vlen));                 \
      P_right = PSTEP_I_SRA (P_log_4, W2, 4, P_right, (vlen));                 \
      /* P_right in Q76 */                                                    \
      P_left = PSTEP_I_SRA (P_log_2, P_log_3, 5, W2, (vlen));                  \
      P_left = PSTEP_I_SRA (P_log_1, W2, 4, P_left, (vlen));                   \
      P_left = PSTEP_I_SRA (P_log_0, W2, 5, P_left, (vlen));                   \
      /* P_left in Q66 */                                                     \
      W8 = __riscv_vsmul (W2, W2, 1, (vlen));           /* Q67 */             \
      W8 = __riscv_vsmul (W8, W8, 1, (vlen));           /* Q71 */             \
      P_right = __riscv_vsmul (P_right, W8, 1, (vlen)); /* Q84 */             \
      P_right = __riscv_vsra (P_right, 18, (vlen));     /* Q66 */             \
      P_left = __riscv_vadd (P_left, P_right, (vlen));  /* Q66 */             \
                                                                             \
      VINT W3 = __riscv_vsmul (W2, W, 1, (vlen));     /* Q66 */                \
      P_left = __riscv_vsmul (P_left, W3, 1, (vlen)); /* Q69 */                \
      P_left = VMERGE_VVM_I (VMVI_VX (0, (vlen)), P_left, lgamma_mask,         \
                             (vlen));                                         \
      VFLOAT poly_hi = __riscv_vfcvt_f (P_left, (vlen));                       \
      P_left = __riscv_vsub (P_left, __riscv_vfcvt_x (poly_hi, (vlen)), (vlen));\
      VFLOAT poly_lo = __riscv_vfcvt_f (P_left, (vlen));                       \
      poly_hi = LGAMMA_VFMUL_VF (poly_hi, 0x1.0p-69, (vlen));                  \
      poly_lo = LGAMMA_VFMUL_VF (poly_lo, 0x1.0p-69, (vlen));                  \
                                                                             \
      /* n*log(2) - c + w + poly is the desired result */                     \
      VFLOAT A, B;                                                            \
      A = LGAMMA_VFMUL_VF (n_flt, LOG2_HI, (vlen)); /* exact */               \
      A = LGAMMA_VFSUB (A, (c), (vlen)); /* exact due to A's range */         \
      w_hi = LGAMMA_VFADD (w_hi, w_hi, (vlen));                               \
      w_lo = LGAMMA_VFADD (w_lo, w_lo, (vlen));                               \
      LGAMMA_FAST2SUM_M (lgamma_mask, A, w_hi, B, (y_lo), (vlen));             \
      w_lo = LGAMMA_VFADD ((y_lo), w_lo, (vlen));                             \
      w_lo = LGAMMA_VFMACC_VF_M (lgamma_mask, w_lo, LOG2_LO, n_flt,           \
                                 (vlen));                                     \
      poly_lo = LGAMMA_VFADD (poly_lo, w_lo, (vlen));                         \
      LGAMMA_FAST2SUM_M (lgamma_mask, B, poly_hi, (y_hi), (y_lo), (vlen));    \
      (y_lo) = LGAMMA_VFADD ((y_lo), poly_lo, (vlen));                        \
    }                                                                          \
  while (0)

/* Use Stirling approximation with correction when x >= 9/4.
   On input logx_hi, logx_lo is log(x)-1. Result is returned in y_hi, y_lo. */
#define LGAMMA_LOG_STIRLING(x, logx_hi, logx_lo, y_hi, y_lo, expo_adj, vlen)   \
  do                                                                           \
    {                                                                          \
      VFLOAT x_in = (x);                                                       \
      VBOOL adjust_x = __riscv_vmfge (x_in, 0x1.0p+200, (vlen));               \
      (expo_adj) = __riscv_vmerge ((expo_adj), 100, adjust_x, (vlen));         \
      VINT m = __riscv_vrsub ((expo_adj), EXP_BIAS, (vlen));                   \
      VFLOAT scale = I_AS_F (__riscv_vsll (m, MAN_LEN, (vlen)));               \
      x_in = LGAMMA_VFMUL (x_in, scale, (vlen));                               \
      VFLOAT w_hi, w_lo;                                                       \
      w_hi = LGAMMA_VFSUB_VF_M (lgamma_mask, x_in, 0x1.0p-1, (vlen));           \
      w_lo = LGAMMA_VFSUB (x_in, w_hi, (vlen));                                \
      w_lo = LGAMMA_VFSUB_VF_M (lgamma_mask, w_lo, 0x1.0p-1, (vlen));           \
      LGAMMA_PROD_X2Y2_M (lgamma_mask, w_hi, w_lo, (logx_hi), (logx_lo),       \
                          (y_hi), (y_lo), (vlen));                             \
    }                                                                          \
  while (0)

/* Compute log(x*r) + log(|sin(pi r)/(pi r)|) where x = N + r, |r| <= 1/2.
   Used for negative arguments: log(|x sin(pi x)/pi|) term. */
#define LGAMMA_LOGSIN(x, y_hi, y_lo, vy_special, special_args, vlen)           \
  do                                                                           \
    {                                                                          \
      VFLOAT x_in = __riscv_vfmin ((x), 0x1.0p+52, (vlen));                    \
      VFLOAT n_flt;                                                            \
      VINT n = __riscv_vfcvt_x (x_in, (vlen));                                 \
      n_flt = __riscv_vfcvt_f (n, (vlen));                                     \
      VFLOAT r = LGAMMA_VFSUB (x_in, n_flt, (vlen));                            \
      VBOOL pole = __riscv_vmand (lgamma_mask,                              \
                                  __riscv_vmfeq (r, fp_posZero, (vlen)),       \
                                  (vlen));                                     \
      if (__riscv_vcpop (pole, (vlen)) > 0)                                    \
        {                                                                      \
          r = __riscv_vfmerge (r, 0x1.0p-1, pole, (vlen));                      \
          (special_args) = __riscv_vmor ((special_args), pole, (vlen));        \
          (vy_special) = __riscv_vfmerge ((vy_special), fp_posInf, pole,       \
                                          (vlen));                             \
        }                                                                      \
      VFLOAT rsq = LGAMMA_VFMUL (r, r, (vlen));                                \
      VFLOAT rsq_lo = LGAMMA_VFMSUB (r, r, rsq, (vlen));                       \
      VINT Rsq = __riscv_vfcvt_x (LGAMMA_VFMUL_VF (rsq, 0x1.0p63, (vlen)),    \
                                  (vlen));                                      \
      Rsq = __riscv_vadd (                                                     \
          Rsq, __riscv_vfcvt_x (LGAMMA_VFMUL_VF (rsq_lo, 0x1.0p63, (vlen)),    \
                                (vlen)),                                       \
          (vlen));                                                             \
      VINT P_right                                                            \
          = PSTEP_I (P_logsin_8, Rsq,                                          \
                     PSTEP_I (P_logsin_9, Rsq,                                 \
                              PSTEP_I (P_logsin_10, Rsq,                       \
                                       PSTEP_I (P_logsin_11, Rsq,              \
                                                PSTEP_I (P_logsin_12, Rsq,     \
                                                         PSTEP_I (P_logsin_13,\
                                                                  P_logsin_14, \
                                                                  Rsq, (vlen)),\
                                                         (vlen)),             \
                                                (vlen)),                      \
                                       (vlen)),                               \
                              (vlen)),                                        \
                     (vlen));                                                 \
                                                                             \
      VINT P_left = PSTEP_I (                                                 \
          P_logsin_0, Rsq,                                                    \
          PSTEP_I (                                                           \
              P_logsin_1, Rsq,                                                \
              PSTEP_I (P_logsin_2, Rsq,                                       \
                       PSTEP_I (P_logsin_3, Rsq,                              \
                                PSTEP_I (P_logsin_4, Rsq,                     \
                                         PSTEP_I (P_logsin_5, Rsq,            \
                                                  PSTEP_I (P_logsin_6,        \
                                                           P_logsin_7, Rsq,   \
                                                           (vlen)),           \
                                                  (vlen)),                    \
                                         (vlen)),                             \
                                (vlen)),                                      \
                       (vlen)),                                               \
              (vlen)),                                                        \
          (vlen));                                                            \
      VINT R16 = __riscv_vsmul (Rsq, Rsq, 1, (vlen));                          \
      R16 = __riscv_vsmul (R16, R16, 1, (vlen));                               \
      R16 = __riscv_vsmul (R16, R16, 1, (vlen));                               \
      P_right = __riscv_vsmul (R16, P_right, 1, (vlen));                       \
      P_left = __riscv_vadd (P_left, P_right, (vlen));                         \
      P_left = VMERGE_VVM_I (VMVI_VX (0, (vlen)), P_left, lgamma_mask,         \
                             (vlen));                                         \
      VFLOAT z_hi = __riscv_vfcvt_f (P_left, (vlen));                          \
      P_right = __riscv_vfcvt_x (z_hi, (vlen));                                \
      VFLOAT z_lo = __riscv_vfcvt_f (                                         \
          __riscv_vsub (P_left, P_right, (vlen)), (vlen));                     \
      z_hi = LGAMMA_VFMUL_VF (z_hi, 0x1.0p-61, (vlen));                        \
      z_lo = LGAMMA_VFMUL_VF (z_lo, 0x1.0p-61, (vlen));                        \
      VFLOAT ls_hi, ls_lo;                                                    \
      LGAMMA_PROD_X2Y2_M (lgamma_mask, z_hi, z_lo, rsq, rsq_lo, ls_hi,        \
                          ls_lo, (vlen));                                     \
      /* At this point we have log|sin(pi r)/(pi r)| */                       \
                                                                             \
      /* we now compute log(|x r|); 2^(-60) <= x <= 2^52 by design */         \
      VFLOAT xr_hi, xr_lo;                                                    \
      r = __riscv_vfsgnj (r, fp_posOne, (vlen));                               \
      LGAMMA_PROD_X1Y1_M (lgamma_mask, r, x_in, xr_hi, xr_lo, (vlen));         \
      VFLOAT logx_hi, logx_lo;                                                 \
      VFLOAT c = VFMV_VF (fp_posZero, (vlen));                                 \
      LGAMMA_LOG (xr_hi, xr_lo, c, logx_hi, logx_lo, (vlen));                  \
      VFLOAT S_hi, S_lo;                                                      \
      LGAMMA_KNUTH2SUM_M (lgamma_mask, ls_hi, logx_hi, S_hi, S_lo, (vlen));    \
      logx_lo = LGAMMA_VFADD (logx_lo, ls_lo, (vlen));                         \
      (y_lo) = LGAMMA_VFADD (S_lo, logx_lo, (vlen));                          \
      (y_hi) = S_hi;                                                          \
    }                                                                         \
  while (0)

/* LogGamma based on Stirling formula:
   logGamma(x) ~ (x-1/2)*(log(x)-1) + poly(1/x) */
#define LOG_STIRLING_CORRECTION(x, y_hi, y_lo, vlen)                           \
  do                                                                           \
    {                                                                          \
      VFLOAT x_in = __riscv_vfmin ((x), 0x1.0p80, (vlen));                     \
      VFLOAT z_hi = LGAMMA_VFRDIV_VF_M (lgamma_mask, x_in, 0x1.0p0, (vlen));  \
      VFLOAT z_lo = VFMV_VF (fp_posOne, (vlen));                               \
      z_lo = LGAMMA_VFNMSUB (x_in, z_hi, z_lo, (vlen));                        \
      z_lo = LGAMMA_VFMUL (z_hi, z_lo, (vlen));                                \
      z_hi = LGAMMA_VFMUL_VF (z_hi, 0x1.0p63, (vlen));                         \
      z_lo = LGAMMA_VFMUL_VF (z_lo, 0x1.0p63, (vlen));                         \
      VINT R = __riscv_vfcvt_x (z_hi, (vlen));                                 \
      R = __riscv_vadd (R, __riscv_vfcvt_x (z_lo, (vlen)), (vlen));            \
      VINT P_SC, Q_SC;                                                        \
      /* R is 1/x in Q63 */                                                   \
      P_SC = PSTEP_I (                                                        \
          P_LS_corr_4, R,                                                     \
          PSTEP_I (P_LS_corr_5, R,                                            \
                   PSTEP_I (P_LS_corr_6, R,                                   \
                            PSTEP_I (P_LS_corr_7, P_LS_corr_8, R, (vlen)),     \
                            (vlen)),                                          \
                   (vlen)),                                                   \
          (vlen));                                                            \
      P_SC = PSTEP_I (                                                        \
          P_LS_corr_0, R,                                                     \
          PSTEP_I (P_LS_corr_1, R,                                            \
                   PSTEP_I (P_LS_corr_2, R,                                   \
                            PSTEP_I (P_LS_corr_3, R, P_SC, (vlen)), (vlen)),  \
                   (vlen)),                                                   \
          (vlen));                                                            \
      Q_SC = PSTEP_I (                                                        \
          Q_LS_corr_4, R,                                                     \
          PSTEP_I (Q_LS_corr_5, R,                                            \
                   PSTEP_I (Q_LS_corr_6, R,                                   \
                            PSTEP_I (Q_LS_corr_7, Q_LS_corr_8, R, (vlen)),     \
                            (vlen)),                                          \
                   (vlen)),                                                   \
          (vlen));                                                            \
      Q_SC = PSTEP_I (                                                        \
          Q_LS_corr_0, R,                                                     \
          PSTEP_I (Q_LS_corr_1, R,                                            \
                   PSTEP_I (Q_LS_corr_2, R,                                   \
                            PSTEP_I (Q_LS_corr_3, R, Q_SC, (vlen)), (vlen)),  \
                   (vlen)),                                                   \
          (vlen));                                                            \
      VFLOAT p_hi, p_lo, q_hi, q_lo;                                          \
      VINT P_tmp, Q_tmp;                                                      \
      P_SC = VMERGE_VVM_I (VMVI_VX (0, (vlen)), P_SC, lgamma_mask, (vlen));    \
      Q_SC = VMERGE_VVM_I (VMVI_VX (0, (vlen)), Q_SC, lgamma_mask, (vlen));    \
      p_hi = __riscv_vfcvt_f (P_SC, (vlen));                                   \
      P_tmp = __riscv_vfcvt_x (p_hi, (vlen));                                  \
      p_lo = __riscv_vfcvt_f (                                                 \
          __riscv_vsub (P_SC, P_tmp, (vlen)), (vlen));                         \
      q_hi = __riscv_vfcvt_f (Q_SC, (vlen));                                   \
      Q_tmp = __riscv_vfcvt_x (q_hi, (vlen));                                  \
      q_lo = __riscv_vfcvt_f (                                                 \
          __riscv_vsub (Q_SC, Q_tmp, (vlen)), (vlen));                         \
      LGAMMA_ACC_DIV2_N2D2_M (lgamma_mask, p_hi, p_lo, q_hi, q_lo,             \
                              (y_hi), (y_lo), (vlen));                         \
    }                                                                          \
  while (0)

#endif

#if (LMUL == 1)
# ifndef RVV_LIBMVEC_LGAMMA_HELPERS_M1
#  define RVV_LIBMVEC_LGAMMA_HELPERS_M1
#  define RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
# endif
#elif (LMUL == 2)
# ifndef RVV_LIBMVEC_LGAMMA_HELPERS_M2
#  define RVV_LIBMVEC_LGAMMA_HELPERS_M2
#  define RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
# endif
#elif (LMUL == 4)
# ifndef RVV_LIBMVEC_LGAMMA_HELPERS_M4
#  define RVV_LIBMVEC_LGAMMA_HELPERS_M4
#  define RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
# endif
#elif (LMUL == 8)
# ifndef RVV_LIBMVEC_LGAMMA_HELPERS_M8
#  define RVV_LIBMVEC_LGAMMA_HELPERS_M8
#  define RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
# endif
#endif
#ifdef RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
/* This helper is instantiated once per importer LMUL.  Keeping the vector
   types in the generated name avoids reusing the m1 signature for m2/m4/m8.  */
static inline void
RVV_LIBMVEC_IMPL_NAME (lgamma_zero_apply, LMUL)
  (VFLOAT sx, VBOOL mask, double root_hi, double root_lo1, double root_lo2,
   const double c[][2], int c_len, double sc,
   VFLOAT *out_hi, VFLOAT *out_lo, size_t vlen)
{
  if (__riscv_vcpop (mask, vlen) == 0)
    return;
  VFLOAT z_hi = VFMV_VF (fp_posZero, vlen);
  VFLOAT z_lo = VFMV_VF (fp_posZero, vlen);
  VFLOAT x0_hi = VFMV_VF (root_hi, vlen);
  VFLOAT x0_lo1 = VFMV_VF (root_lo1, vlen);
  VFLOAT x0_lo2 = VFMV_VF (root_lo2, vlen);
  z_hi = VFADD_VV_M (mask, sx, x0_hi, vlen);
  VFLOAT x_hat = VFSUB_VV_M (mask, z_hi, x0_hi, vlen);
  VFLOAT y_hat = VFSUB_VV_M (mask, z_hi, x_hat, vlen);
  VFLOAT x_delta = VFSUB_VV_M (mask, sx, x_hat, vlen);
  VFLOAT y_delta = VFSUB_VV_M (mask, x0_hi, y_hat, vlen);
  z_lo = VFADD_VV_M (mask, x_delta, y_delta, vlen);
  z_lo = VFADD_VV_M (mask, z_lo, x0_lo1, vlen);
  z_lo = VFADD_VV_M (mask, z_lo, x0_lo2, vlen);

  VFLOAT scv = VFMV_VF (sc, vlen);
  VFLOAT sh = VFMUL_VV_M (mask, z_hi, scv, vlen);
  VFLOAT sl = VFMUL_VV_M (mask, z_lo, scv, vlen);

  VFLOAT ph = VFMV_VF (c[c_len - 1][0], vlen);
  VFLOAT pl = VFMV_VF (c[c_len - 1][1], vlen);
  for (int i = c_len - 2; i >= 0; --i)
    {
      VFLOAT prod_hi, prod_lo;
      prod_hi = VFMUL_VV_M (mask, sh, ph, vlen);
      prod_lo = VFMSUB_VV_M (mask, sh, ph, prod_hi, vlen);
      prod_lo = VFMACC_VV_M (mask, prod_lo, sh, pl, vlen);
      prod_lo = VFMACC_VV_M (mask, prod_lo, sl, ph, vlen);
      VFLOAT ch = VFMV_VF (c[i][0], vlen);
      VFLOAT cl = VFMV_VF (c[i][1], vlen);
      VFLOAT sum_hi = VFADD_VV_M (mask, prod_hi, ch, vlen);
      x_hat = VFSUB_VV_M (mask, sum_hi, ch, vlen);
      y_hat = VFSUB_VV_M (mask, sum_hi, x_hat, vlen);
      x_delta = VFSUB_VV_M (mask, prod_hi, x_hat, vlen);
      y_delta = VFSUB_VV_M (mask, ch, y_hat, vlen);
      VFLOAT sum_lo = VFADD_VV_M (mask, x_delta, y_delta, vlen);
      sum_lo = VFADD_VV_M (mask, sum_lo, prod_lo, vlen);
      sum_lo = VFADD_VV_M (mask, sum_lo, cl, vlen);
      ph = sum_hi;
      pl = sum_lo;
    }

  VFLOAT lgamma_zero_hi = VFMUL_VV_M (mask, z_hi, ph, vlen);
  VFLOAT lgamma_zero_lo = VFMSUB_VV_M (mask, z_hi, ph,
                                       lgamma_zero_hi, vlen);
  lgamma_zero_lo = VFMACC_VV_M (mask, lgamma_zero_lo, z_hi, pl, vlen);
  lgamma_zero_lo = VFMACC_VV_M (mask, lgamma_zero_lo, z_lo, ph, vlen);
  *out_hi = __riscv_vmerge (*out_hi, lgamma_zero_hi, mask, vlen);
  *out_lo = __riscv_vmerge (*out_lo, lgamma_zero_lo, mask, vlen);
}
# undef RVV_LIBMVEC_LGAMMA_HELPERS_THIS_LMUL
#endif
#endif

/* Exact small binary32 operations shared by clean and source kernels.  */
#if BIT_WIDTH == 32
# if (LMUL == 1)
#  ifndef RVV_F32_COMMON_HELPERS_M1
#   define RVV_F32_COMMON_HELPERS_M1
#   define RVV_F32_COMMON_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_F32_COMMON_HELPERS_M2
#   define RVV_F32_COMMON_HELPERS_M2
#   define RVV_F32_COMMON_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_F32_COMMON_HELPERS_M4
#   define RVV_F32_COMMON_HELPERS_M4
#   define RVV_F32_COMMON_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_F32_COMMON_HELPERS_M8
#   define RVV_F32_COMMON_HELPERS_M8
#   define RVV_F32_COMMON_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_F32_COMMON_HELPERS_THIS_LMUL
static inline VINT
RVV_LIBMVEC_IMPL_NAME (rvv_select_i_lt_f32, LMUL)
  (VFLOAT f0, VFLOAT f1, VINT x, VINT y, size_t vlen)
{
  VBOOL m = __riscv_vmflt (f0, f1, vlen);
  return __riscv_vmerge (y, x, m, vlen);
}
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (rvv_vrint_f32, LMUL)
  (VFLOAT value, size_t vlen)
{
  VINT rounded = RVV_VFCVT_X_F_RNE (value, vlen);
  return __riscv_vfcvt_f (rounded, vlen);
}

#  undef RVV_F32_COMMON_HELPERS_THIS_LMUL
# endif
#endif

/* Quarter-period reduction shared by clean sinf, fastsinf and tanf.  */
#if (BIT_WIDTH == 32) \
  && defined (RVV_NEED_CLEAN_F32_TRIG_REDUCTION)
# if (LMUL == 1)
#  ifndef RVV_LIBMVEC_TRIG_HELPERS_M1
#   define RVV_LIBMVEC_TRIG_HELPERS_M1
#   define RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_LIBMVEC_TRIG_HELPERS_M2
#   define RVV_LIBMVEC_TRIG_HELPERS_M2
#   define RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_LIBMVEC_TRIG_HELPERS_M4
#   define RVV_LIBMVEC_TRIG_HELPERS_M4
#   define RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_LIBMVEC_TRIG_HELPERS_M8
#   define RVV_LIBMVEC_TRIG_HELPERS_M8
#   define RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sinf_or_sign, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  VUINT sign = __riscv_vand (F_AS_U_F (y), 0x80000000u, vlen);
  return U_AS_F_F (__riscv_vor (F_AS_U_F (x), sign, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sinf_reduce_quarter, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_LIBMVEC_IMPL_NAME (sinf_or_sign, LMUL) (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_LIBMVEC_IMPL_NAME (sinf_or_sign, LMUL) (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x,
     vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF_F (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x, vlen), vlen);
  return rem;
}
#  undef RVV_LIBMVEC_TRIG_HELPERS_THIS_LMUL
# endif
#endif

/* Shared source-original SLEEF single-precision trig reduction helpers.
   The source comparator's conversion points, bitwise sign insertion, table
   reduction, and late special-value behavior remain unchanged.  */
#if (BIT_WIDTH == 32) \
  && defined (RVV_NEED_SLEEF_F32_TRIG_REDUCTION)
# if (LMUL == 1)
#  ifndef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M1
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M1
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL
#  endif
# elif (LMUL == 2)
#  ifndef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M2
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M2
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL
#  endif
# elif (LMUL == 4)
#  ifndef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M4
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M4
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL
#  endif
# elif (LMUL == 8)
#  ifndef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M8
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_M8
#   define RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL
#  endif
# endif
# ifdef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sinf_original_or_sign, LMUL)
  (VFLOAT x, VFLOAT y, size_t vlen)
{
  VUINT sign = __riscv_vand (F_AS_U (y), 0x80000000u, vlen);
  return U_AS_F (__riscv_vor (F_AS_U (x), sign, vlen));
}

static inline VFLOAT
RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_quarter, LMUL)
  (VFLOAT x, VINT *q_out, size_t vlen)
{
  VFLOAT c = __riscv_vfsgnjx
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_TWO23, vlen), x, vlen);
  VFLOAT x4 = VFMUL_VF (x, RVV_SLEEF_TRIG_REMPI_FOUR, vlen);
  VFLOAT abs_x4 = __riscv_vfabs (x4, vlen);
  VFLOAT abs_x = __riscv_vfabs (x, vlen);
  VBOOL x4_large = __riscv_vmfgt
    (abs_x4, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);
  VBOOL x_large = __riscv_vmfgt
    (abs_x, RVV_SLEEF_TRIG_REMPI_TWO23, vlen);

  VFLOAT rint4x = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_FOUR, vlen), x, c, vlen);
  rint4x = VFSUB_VV (rint4x, c, vlen);
  rint4x = RVV_LIBMVEC_IMPL_NAME (sinf_original_or_sign, LMUL)
    (rint4x, x, vlen);
  rint4x = VMERGE_VVM_F (rint4x, x4, x4_large, vlen);

  VFLOAT rintx = VFADD_VV (x, c, vlen);
  rintx = VFSUB_VV (rintx, c, vlen);
  rintx = RVV_LIBMVEC_IMPL_NAME (sinf_original_or_sign, LMUL)
    (rintx, x, vlen);
  rintx = VMERGE_VVM_F (rintx, x, x_large, vlen);

  VFLOAT rem = VFMADD_VV
    (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_QUARTER, vlen), rint4x, x,
              vlen);
  *q_out = __riscv_vfcvt_rtz_x
    (VFMADD_VV (VFMV_VF (RVV_SLEEF_TRIG_REMPI_NEG_FOUR, vlen), rintx,
                rint4x,
                vlen), vlen);
  return rem;
}

static inline void
RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_large, LMUL)
  (VFLOAT a, VFLOAT *out_hi, VFLOAT *out_lo, VINT *out_q, size_t vlen)
{
  const struct v_sleef_trig_data *trig_data =
    ptr_barrier (&__v_sleef_trig_data);

  /* Exact exponent-field adapter for the source reduction.  */
  VINT ex;
  {
    VUINT bits = F_AS_U (a);
    VUINT exponent = __riscv_vsrl (bits, MAN_LEN_F, vlen);
    exponent = __riscv_vand (exponent, 0xff, vlen);
    ex = __riscv_vsub (U_AS_I (exponent), EXP_BIAS_F, vlen);
  }
  ex = __riscv_vsub (ex, RVV_SLEEF_SINF_REMPI_EXPONENT_BIAS, vlen);

  VBOOL scale_mask = __riscv_vmsgt
    (ex, RVV_SLEEF_SINF_REMPI_SCALE_LIMIT, vlen);
  VINT scale = VMERGE_VXM_I
    (VMVI_VX (0, vlen), RVV_SLEEF_SINF_REMPI_SCALE, scale_mask, vlen);
  /* Apply the exact exponent-field scaling.  */
  a = I_AS_F (__riscv_vadd
    (F_AS_I (a), __riscv_vsll (scale, MAN_LEN_F, vlen), vlen));

  VBOOL ex_negative = __riscv_vmslt (ex, 0, vlen);
  ex = VMERGE_VXM_I (ex, 0, ex_negative, vlen);
  ex = __riscv_vsll (ex, RVV_SLEEF_SINF_REMPI_INDEX_SHIFT, vlen);
  VUINT ex_u = I_AS_U (ex);

  VFLOAT tab0 = RVV_VLUXEI32_F (trig_data->rempitabsp, ex_u, vlen);
  VFLOAT tab1 = RVV_VLUXEI32_F (trig_data->rempitabsp + 1, ex_u, vlen);
  VFLOAT tab2 = RVV_VLUXEI32_F (trig_data->rempitabsp + 2, ex_u, vlen);
  VFLOAT tab3 = RVV_VLUXEI32_F (trig_data->rempitabsp + 3, ex_u, vlen);

  VFLOAT x_hi, x_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab0, x_hi, x_lo, vlen);
  VINT q;
  x_hi = RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_quarter, LMUL)
    (x_hi, &q, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  VFLOAT y_hi, y_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF (a, tab1, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  VINT q_part;
  x_hi = RVV_LIBMVEC_IMPL_NAME (sinf_original_reduce_quarter, LMUL)
    (x_hi, &q_part, vlen);
  q = __riscv_vadd (q, q_part, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  y_hi = tab2;
  y_lo = tab3;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (y_hi, y_lo, a, y_hi, y_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (x_hi, x_lo, y_hi, y_lo, x_hi, x_lo, vlen);
  RVV_SLEEF_ORIG_DFNORMALIZE_VF2 (x_hi, x_lo, x_hi, x_lo, vlen);

  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (x_hi, x_lo,
     VFMV_VF ((RVV_PI_F * 2.0f), vlen),
     VFMV_VF ((RVV_PI_LO_F * 2.0f), vlen),
     x_hi, x_lo, vlen);

  VBOOL a_small = __riscv_vmflt
    (__riscv_vfabs (a, vlen), RVV_SLEEF_TRIG_REMPI_SMALL, vlen);
  x_hi = VMERGE_VVM_F (x_hi, a, a_small, vlen);
  x_lo = VFMERGE_VFM (x_lo, fp_posZero_f, a_small, vlen);

  *out_hi = x_hi;
  *out_lo = x_lo;
  *out_q = q;
}

#  undef RVV_SLEEF_ORIG_SINF_REDUCTION_ADAPTERS_THIS_LMUL
# endif
#endif
