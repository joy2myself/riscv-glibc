/* Clean-side numerical repairs for the SLEEF-shaped lgammaf kernel.

   The source-original comparator deliberately does not include this file.
   These helpers only repair two finite-path limitations exposed by the clean
   accuracy contract:

     - gammafk's large-x expression can overflow in the intermediate
       (x - 1/2) * log(x), even when the final lgammaf value is finite;
     - the source double-float result loses relative accuracy when reflection
       terms cancel near a zero of lgammaf.

   The first repair is an algebraically equivalent source o2 expression with
   the source Stirling correction retained.  The second is the pre-existing
   widened-double reflection evaluation used by the clean implementation.
   Neither helper changes the source-original skeleton.  */

/* The subnormal branch of lgammaf is dominated by -log(|x|).  Evaluate the
   logarithm after an exact power-of-two normalization so that the source
   logk2f helper never sees a subnormal operand.  The correction terms in
   log(Gamma(x)) are below the binary32 ulp throughout the subnormal range.  */
static inline VFLOAT
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_subnormal, LMUL)
  (VFLOAT x, size_t vlen)
{
  VFLOAT ax = __riscv_vfabs (x, vlen);
  VFLOAT normalized = VFMUL_VF (ax, 0x1.0p+64f, vlen);
  VFLOAT zero = VFMV_VF (0.0f, vlen);
  VFLOAT log_hi, log_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m, LMUL)
    (normalized, zero, &log_hi, &log_lo, vlen);

  VFLOAT ln2_hi = VFMV_VF (RVV_TGAMMAF_LOG_L2U_F, vlen);
  VFLOAT ln2_lo = VFMV_VF (RVV_TGAMMAF_LOG_L2L_F, vlen);
  VFLOAT sixty_four = VFMV_VF (64.0f, vlen);
  VFLOAT scale_hi, scale_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (ln2_hi, ln2_lo, sixty_four, scale_hi, scale_lo, vlen);
  RVV_SLEEF_ORIG_DFNEG_VF2
    (scale_hi, scale_lo, scale_hi, scale_lo, vlen);

  VFLOAT neg_log_hi, neg_log_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (log_hi, log_lo, scale_hi, scale_lo, neg_log_hi, neg_log_lo, vlen);
  RVV_SLEEF_ORIG_DFNEG_VF2
    (neg_log_hi, neg_log_lo, neg_log_hi, neg_log_lo, vlen);
  return VFADD_VV (neg_log_hi, neg_log_lo, vlen);
}

/* Evaluate the source o2 lgammaf expression after the overflow-prone product
   has been rearranged as

       x * (log(x) - 1) - log(x) / 2 + log(sqrt(2*pi)).

   The `u*t` correction is still formed as gammafk would form clld, so the
   returned value remains the same source approximation away from overflow.  */
static inline VFLOAT
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_large, LMUL)
  (VFLOAT x, size_t vlen)
{
  VFLOAT log_arg_scale = VFMV_VF (0.0625f, vlen);
  VBOOL scale_log = __riscv_vmfge (x, 0x1.2ced32p+126f, vlen);
  VFLOAT scaled_x = VFMUL_VV (x, log_arg_scale, vlen);
  VFLOAT log_arg = VMERGE_VVM_F (x, scaled_x, scale_log, vlen);

  VFLOAT zero = VFMV_VF (0.0f, vlen);
  VFLOAT log_hi, log_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m, LMUL)
    (log_arg, zero, &log_hi, &log_lo, vlen);

  VFLOAT ln2_hi = VFMV_VF (RVV_TGAMMAF_LOG_L2U_F, vlen);
  VFLOAT ln2_lo = VFMV_VF (RVV_TGAMMAF_LOG_L2L_F, vlen);
  VFLOAT four = VFMV_VF (4.0f, vlen);
  VFLOAT four_ln2_hi, four_ln2_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (ln2_hi, ln2_lo, four, four_ln2_hi, four_ln2_lo, vlen);
  VFLOAT scaled_log_hi, scaled_log_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (log_hi, log_lo, four_ln2_hi, four_ln2_lo,
     scaled_log_hi, scaled_log_lo, vlen);
  log_hi = VMERGE_VVM_F (log_hi, scaled_log_hi, scale_log, vlen);
  log_lo = VMERGE_VVM_F (log_lo, scaled_log_lo, scale_log, vlen);

  VFLOAT log_minus_one_hi, log_minus_one_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (log_hi, log_lo, -1.0f, log_minus_one_hi, log_minus_one_lo, vlen);
  VFLOAT xterm_hi, xterm_lo;
  RVV_SLEEF_ORIG_DFMUL_VF2_VF
    (log_minus_one_hi, log_minus_one_lo, x, xterm_hi, xterm_lo, vlen);

  VFLOAT half_log_hi, half_log_lo;
  RVV_SLEEF_ORIG_DFSCALE_VF2_VF
    (log_hi, log_lo, VFMV_VF (0.5f, vlen),
     half_log_hi, half_log_lo, vlen);
  RVV_SLEEF_ORIG_DFNEG_VF2
    (half_log_hi, half_log_lo, half_log_hi, half_log_lo, vlen);

  VFLOAT sum_hi, sum_lo;
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (xterm_hi, xterm_lo, half_log_hi, half_log_lo,
     sum_hi, sum_lo, vlen);

  VFLOAT c_hi = VFMV_VF ((float) RVV_TGAMMAF_LOG_SQRT_2PI_D, vlen);
  VFLOAT c_lo = VFMV_VF
    ((float) (RVV_TGAMMAF_LOG_SQRT_2PI_D
              - (double) (float) RVV_TGAMMAF_LOG_SQRT_2PI_D), vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (sum_hi, sum_lo, c_hi, c_lo, sum_hi, sum_lo, vlen);

  /* Retain gammafk's o2 correction: clln/clld = 1 + u*t.  */
  VFLOAT t = VFDIV_VV (VFMV_VF (1.0f, vlen), x, vlen);
  VFLOAT u = VFMV_VF (+0.000839498720672087279971000786f, vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (-5.17179090826059219329394422e-05f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (-0.000592166437353693882857342347f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (+6.97281375836585777403743539e-05f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (+0.000784039221720066627493314301f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (-0.000229472093621399176949318732f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (-0.002681327160493827160473958490f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (+0.003472222222222222222175164840f, vlen), vlen);
  u = VFMADD_VV
    (u, t, VFMV_VF (+0.083333333333333333335592087900f, vlen), vlen);

  VFLOAT ut_hi, ut_lo, denominator_hi, denominator_lo;
  RVV_SLEEF_ORIG_DFMUL_VF_VF
    (u, t, ut_hi, ut_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF_S
    (ut_hi, ut_lo, 1.0f, denominator_hi, denominator_lo, vlen);
  VFLOAT log_correction_hi, log_correction_lo;
  RVV_SLEEF_ORIG_MANGLE (rvv_sleef_orig_tgammaf_logk2f_m, LMUL)
    (denominator_hi, denominator_lo,
     &log_correction_hi, &log_correction_lo, vlen);
  RVV_SLEEF_ORIG_DFADD2_VF2_VF2
    (sum_hi, sum_lo, log_correction_hi, log_correction_lo,
     sum_hi, sum_lo, vlen);
  return VFADD_VV (sum_hi, sum_lo, vlen);
}

#if LMUL == 1 || LMUL == 2 || LMUL == 4
/* A widened-double reflection repair.  It keeps the old clean association:
   log(pi) - log(abs(sin(pi*x))) - lgamma(1-x), with a double Lanczos and
   double range-reduced sine.  The caller sanitizes lanes outside its mask.  */
static inline RVV_LGAMMAF_REPAIR_D
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_logd, LMUL)
  (RVV_LGAMMAF_REPAIR_D x, size_t vlen)
{
  RVV_LGAMMAF_REPAIR_U ix = RVV_LGAMMAF_REPAIR_F2U (x);
  RVV_LGAMMAF_REPAIR_U exponent = __riscv_vsrl (ix, 52, vlen);
  exponent = __riscv_vand (exponent, 0x7ff, vlen);
  RVV_LGAMMAF_REPAIR_I e = __riscv_vsub
    (RVV_LGAMMAF_REPAIR_U2I (exponent), 1023, vlen);
  RVV_LGAMMAF_REPAIR_U mantissa
    = __riscv_vand (ix, 0x000fffffffffffffull, vlen);
  RVV_LGAMMAF_REPAIR_U m_bits
    = __riscv_vor (mantissa, 0x3ff0000000000000ull, vlen);
  RVV_LGAMMAF_REPAIR_D m = RVV_LGAMMAF_REPAIR_U2F (m_bits);

  RVV_LGAMMAF_REPAIR_B gt_sqrt2 = __riscv_vmfgt
    (m, 0x1.6a09e667f3bcdp+0, vlen);
  m = RVV_LGAMMAF_REPAIR_MERGE
    (m, __riscv_vfmul (m, 0.5, vlen), gt_sqrt2, vlen);
#if LMUL == 1
  RVV_LGAMMAF_REPAIR_I e_inc = __riscv_vmv_v_x_i64m2 (0, vlen);
#elif LMUL == 2
  RVV_LGAMMAF_REPAIR_I e_inc = __riscv_vmv_v_x_i64m4 (0, vlen);
#elif LMUL == 4
  RVV_LGAMMAF_REPAIR_I e_inc = __riscv_vmv_v_x_i64m8 (0, vlen);
#endif
  e_inc = RVV_LGAMMAF_REPAIR_MERGE_X (e_inc, 1, gt_sqrt2, vlen);
  e = __riscv_vadd (e, e_inc, vlen);

  RVV_LGAMMAF_REPAIR_D num = __riscv_vfsub (m, 1.0, vlen);
  RVV_LGAMMAF_REPAIR_D den = __riscv_vfadd (m, 1.0, vlen);
  RVV_LGAMMAF_REPAIR_D y = __riscv_vfdiv (num, den, vlen);
  RVV_LGAMMAF_REPAIR_D y2 = __riscv_vfmul (y, y, vlen);
  RVV_LGAMMAF_REPAIR_D p = RVV_LGAMMAF_REPAIR_MV (1.0 / 17.0, vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 15.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 13.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 11.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 9.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 7.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 5.0, vlen), vlen);
  p = __riscv_vfmadd (p, y2,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 3.0, vlen), vlen);
  RVV_LGAMMAF_REPAIR_D y3 = __riscv_vfmul (y, y2, vlen);
  RVV_LGAMMAF_REPAIR_D logm = __riscv_vfmadd (y3, p, y, vlen);
  logm = __riscv_vfmul (logm, 2.0, vlen);
  RVV_LGAMMAF_REPAIR_D e_f = __riscv_vfcvt_f (e, vlen);
  return __riscv_vfmadd
    (e_f, 0x1.62e42fefa39efp-1, logm, vlen);
}

static inline RVV_LGAMMAF_REPAIR_D
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_sinpi, LMUL)
  (RVV_LGAMMAF_REPAIR_D x, size_t vlen)
{
  RVV_LGAMMAF_REPAIR_D ax = __riscv_vfabs (x, vlen);
  RVV_LGAMMAF_REPAIR_B x_is_int = __riscv_vmfge (ax, 0x1.0p23, vlen);
  RVV_LGAMMAF_REPAIR_D zero = RVV_LGAMMAF_REPAIR_MV (0.0, vlen);
  RVV_LGAMMAF_REPAIR_D x_small
    = RVV_LGAMMAF_REPAIR_MERGE (x, zero, x_is_int, vlen);
  RVV_LGAMMAF_REPAIR_I n = __riscv_vfcvt_x (x_small, vlen);
  RVV_LGAMMAF_REPAIR_D n_f = __riscv_vfcvt_f (n, vlen);
  RVV_LGAMMAF_REPAIR_D r = __riscv_vfsub (x_small, n_f, vlen);

  RVV_LGAMMAF_REPAIR_D y
    = __riscv_vfmul (r, 0x1.921fb54442d18p+1, vlen);
  RVV_LGAMMAF_REPAIR_D s = __riscv_vfmul (y, y, vlen);
  RVV_LGAMMAF_REPAIR_D p
    = RVV_LGAMMAF_REPAIR_MV (1.0 / 355687428096000.0, vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (-1.0 / 1307674368000.0,
                                               vlen), vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 6227020800.0, vlen),
                      vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (-1.0 / 39916800.0, vlen),
                      vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 362880.0, vlen), vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (-1.0 / 5040.0, vlen), vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (1.0 / 120.0, vlen), vlen);
  p = __riscv_vfmadd (p, s,
                      RVV_LGAMMAF_REPAIR_MV (-1.0 / 6.0, vlen), vlen);
  RVV_LGAMMAF_REPAIR_D sin_y = __riscv_vfmadd
    (__riscv_vfmul (y, s, vlen), p, y, vlen);

  RVV_LGAMMAF_REPAIR_I n_odd = __riscv_vand (n, 1, vlen);
  RVV_LGAMMAF_REPAIR_B odd = __riscv_vmsne (n_odd, 0, vlen);
  sin_y = RVV_LGAMMAF_REPAIR_MERGE
    (sin_y, __riscv_vfneg (sin_y, vlen), odd, vlen);
  return RVV_LGAMMAF_REPAIR_MERGE (sin_y, zero, x_is_int, vlen);
}

static inline RVV_LGAMMAF_REPAIR_D
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_lanczos, LMUL)
  (RVV_LGAMMAF_REPAIR_D u, size_t vlen)
{
  RVV_LGAMMAF_REPAIR_D z = __riscv_vfsub (u, 1.0, vlen);
  RVV_LGAMMAF_REPAIR_D a
    = RVV_LGAMMAF_REPAIR_MV (0.99999999999980993, vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (676.5203681218851, vlen),
                       __riscv_vfadd (z, 1.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (-1259.1392167224028, vlen),
                       __riscv_vfadd (z, 2.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (771.32342877765313, vlen),
                       __riscv_vfadd (z, 3.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (-176.61502916214059, vlen),
                       __riscv_vfadd (z, 4.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (12.507343278686905, vlen),
                       __riscv_vfadd (z, 5.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (-0.13857109526572012, vlen),
                       __riscv_vfadd (z, 6.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (9.9843695780195716e-6, vlen),
                       __riscv_vfadd (z, 7.0, vlen), vlen), vlen);
  a = __riscv_vfadd
    (a, __riscv_vfdiv (RVV_LGAMMAF_REPAIR_MV (1.5056327351493116e-7, vlen),
                       __riscv_vfadd (z, 8.0, vlen), vlen), vlen);
  RVV_LGAMMAF_REPAIR_D t = __riscv_vfadd (z, 7.5, vlen);
  RVV_LGAMMAF_REPAIR_D log_t
    = RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_logd, LMUL) (t, vlen);
  RVV_LGAMMAF_REPAIR_D log_a
    = RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_logd, LMUL) (a, vlen);
  RVV_LGAMMAF_REPAIR_D result
    = RVV_LGAMMAF_REPAIR_MV (0.9189385332046727417803297364, vlen);
  result = __riscv_vfmacc
    (result, __riscv_vfadd (z, 0.5, vlen), log_t, vlen);
  result = __riscv_vfsub (result, t, vlen);
  return __riscv_vfadd (result, log_a, vlen);
}

static inline VFLOAT
RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_reflection, LMUL)
  (VFLOAT x, size_t vlen)
{
  RVV_LGAMMAF_REPAIR_D xd
    = RVV_LGAMMAF_REPAIR_VL (x, vlen);
  RVV_LGAMMAF_REPAIR_D one
    = RVV_LGAMMAF_REPAIR_MV (1.0, vlen);
  RVV_LGAMMAF_REPAIR_D u = __riscv_vfsub (one, xd, vlen);
  RVV_LGAMMAF_REPAIR_D sin_abs = __riscv_vfabs
    (RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_sinpi, LMUL)
       (xd, vlen), vlen);
  RVV_LGAMMAF_REPAIR_D result = __riscv_vfsub
    (__riscv_vfsub
       (RVV_LGAMMAF_REPAIR_MV (RVV_TGAMMAF_LOG_PI_D, vlen),
        RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_logd, LMUL)
          (sin_abs, vlen), vlen),
     RVV_LGAMMAF_REPAIR_MANGLE (rvv_lgammaf_clean_lanczos, LMUL)
       (u, vlen), vlen);
  return RVV_LGAMMAF_REPAIR_NC (result, vlen);
}
#endif
