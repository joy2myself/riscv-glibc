/* Single-precision source-original SLEEF atan2-family skeleton.

   The finite and post-core special-value path follows SLEEF xatan2f_u1.
   atan2pif_original is the explicit target-side composition of that radians
   source entry with 1/pi; it is not presented as a native SLEEF entry.  */

#include <vecmath_config.h>
#include <v_math.h>

#ifndef RVV_LIBMVEC_IMPL_NAME
# error "RVV_LIBMVEC_IMPL_NAME must be defined for RVV libmvec skeleton"
#endif
#ifndef V_MATH_FUNC
# error "V_MATH_FUNC must be defined for RVV libmvec skeleton"
#endif

static VFLOAT
RVV_LIBMVEC_IMPL_NAME (V_MATH_FUNC, LMUL) (VFLOAT vy, VFLOAT vx,
                                           size_t vlen)
{
#if defined(COMPILE_FOR_ATAN2PI_ORIGINAL)
  const float angle_pi = RVV_ATANF_NORM_PI_F;
  const float angle_pi2 = RVV_ATANF_NORM_PI2_F;
  const float angle_pi4 = RVV_ATANF_NORM_PI4_F;
#else
  const float angle_pi = RVV_ATANF_PI_F;
  const float angle_pi2 = RVV_ATANF_PI2_F;
  const float angle_pi4 = RVV_ATANF_PI4_F;
#endif
  VBOOL tiny = __riscv_vmflt (__riscv_vfabs (vx, vlen),
                              RVV_ATANF_TINY_F, vlen);
  VFLOAT scale = VFMV_VF (RVV_ATANF_SCALE_F, vlen);
  VFLOAT scaled_x = VFMUL_VV (vx, scale, vlen);
  VFLOAT scaled_y = VFMUL_VV (vy, scale, vlen);
  vx = __riscv_vmerge (vx, scaled_x, tiny, vlen);
  vy = __riscv_vmerge (vy, scaled_y, tiny, vlen);

  VFLOAT r_hi, r_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_atan2kf_u1_m, LMUL)
    (__riscv_vfabs (vy, vlen), VFMV_VF (0.0f, vlen), vx,
     VFMV_VF (0.0f, vlen), &r_hi, &r_lo, vlen);
#if defined(COMPILE_FOR_ATAN2PI_ORIGINAL)
  /* atan2pif is a target-side normalized-pi composition.  Scale the
     source double-float angle before collapsing it; scaling the collapsed
     binary32 result loses two ULPs for a small ordinary finite subset. */
  VFLOAT scaled_hi, scaled_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF2
    (r_hi, r_lo,
     VFMV_VF (RVV_ATANF_INV_PI_HI_F, vlen),
     VFMV_VF (RVV_ATANF_INV_PI_LO_F, vlen),
     scaled_hi, scaled_lo, vlen);
  VFLOAT r = VFADD_VV (scaled_hi, scaled_lo, vlen);
#else
  VFLOAT r = VFADD_VV (r_hi, r_lo, vlen);
#endif
  r = __riscv_vfsgnjx (r, vx, vlen);

  VUINT ux = F_AS_U (vx);
  VUINT uy = F_AS_U (vy);
  VUINT ax = __riscv_vand (ux, 0x7fffffff, vlen);
  VUINT ay = __riscv_vand (uy, 0x7fffffff, vlen);
  VBOOL x_inf = __riscv_vmseq (ax, 0x7f800000, vlen);
  VBOOL y_inf = __riscv_vmseq (ay, 0x7f800000, vlen);
  VBOOL x_zero = __riscv_vmseq (ax, 0, vlen);
  VBOOL y_zero = __riscv_vmseq (ay, 0, vlen);

  /* SLEEF vmulsign(pi/2, x): preserve sign-XOR, not sign-copy.  */
  VFLOAT x_pi2 = __riscv_vfsgnjx
    (VFMV_VF (angle_pi2, vlen), vx, vlen);
  VFLOAT x_inf_part
    = VFSUB_VV (VFMV_VF (angle_pi2, vlen),
                RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_visinf2_m, LMUL)
                  (vx, x_pi2, vlen), vlen);
  VBOOL x_inf_or_zero = __riscv_vmor (x_inf, x_zero, vlen);
  r = __riscv_vmerge (r, x_inf_part, x_inf_or_zero, vlen);

  /* SLEEF vmulsign(pi/4, x): preserve sign-XOR, not sign-copy.  */
  VFLOAT x_pi4 = __riscv_vfsgnjx
    (VFMV_VF (angle_pi4, vlen), vx, vlen);
  VFLOAT y_inf_part
    = VFSUB_VV (VFMV_VF (angle_pi2, vlen),
                RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_visinf2_m, LMUL)
                  (vx, x_pi4, vlen), vlen);
  r = __riscv_vmerge (r, y_inf_part, y_inf, vlen);

  VBOOL x_sign = __riscv_vmsne (__riscv_vand (ux, 0x80000000, vlen), 0,
                                vlen);
  VFLOAT x_negative_pi
    = __riscv_vmerge (VFMV_VF (0.0f, vlen),
                      VFMV_VF (angle_pi, vlen), x_sign,
                      vlen);
  r = __riscv_vmerge (r, x_negative_pi, y_zero, vlen);
  r = __riscv_vfsgnjx (r, vy, vlen);

  VBOOL x_nan = __riscv_vmsgtu (ax, 0x7f800000, vlen);
  VBOOL y_nan = __riscv_vmsgtu (ay, 0x7f800000, vlen);
  VBOOL any_nan = __riscv_vmor (x_nan, y_nan, vlen);
  VUINT nan_bits = __riscv_vmerge (VMVU_VX (0, vlen),
                                    VMVU_VX (~(uint32_t) 0, vlen), any_nan,
                                    vlen);
  r = U_AS_F (__riscv_vor (nan_bits, F_AS_U (r), vlen));

  return r;
}

#undef LMUL
#undef SIMDLEN
