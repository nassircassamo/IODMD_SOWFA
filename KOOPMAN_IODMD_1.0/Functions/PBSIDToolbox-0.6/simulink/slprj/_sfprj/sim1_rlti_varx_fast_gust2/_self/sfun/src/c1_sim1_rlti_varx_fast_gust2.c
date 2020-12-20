/* Include files */

#include "blascompat32.h"
#include "sim1_rlti_varx_fast_gust2_sfun.h"
#include "c1_sim1_rlti_varx_fast_gust2.h"
#include "mwmathutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void initialize_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void initialize_params_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void enable_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void disable_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void set_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_st);
static void finalize_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void sf_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void c1_chartstep_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void initSimStructsc1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber);
static void c1_info_helper(c1_ResolvedFunctionInfo c1_info[102]);
static void c1_b_info_helper(c1_ResolvedFunctionInfo c1_info[102]);
static void c1_eml_error(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static real_T c1_mrdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_A, real_T c1_B);
static real_T c1_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y);
static void c1_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   real_T c1_I[100]);
static void c1_b_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[64]);
static void c1_logspace(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_d1, real_T c1_d2, real_T c1_y[300]);
static void c1_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static void c1_b_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y[2], real_T c1_z[2]);
static void c1_diag(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T c1_v[2], real_T c1_d[4]);
static void c1_b_diag(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, real_T c1_v[3], real_T c1_d[9]);
static void c1_c_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[9]);
static void c1_c_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x[157], real_T c1_y, real_T c1_z[157]);
static real_T c1_norm(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, real_T c1_x[3]);
static void c1_d_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x[936], real_T c1_y, real_T c1_z[936]);
static void c1_b_mrdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_A[153], real_T c1_B, real_T c1_y[153]);
static void c1_damp(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T c1_a[49], real_T c1_h, real_T c1_wn[7], real_T c1_z[7]);
static void c1_eig(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   real_T c1_A[49], creal_T c1_V[7]);
static void c1_eml_matlab_zlartg(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn,
  creal_T *c1_r);
static void c1_eml_matlab_zhgeqz(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_A[49], int32_T c1_ilo, int32_T c1_ihi, real_T
  *c1_info, creal_T c1_alpha1[7], creal_T c1_beta1[7]);
static real_T c1_eml_matlab_zlanhs(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, creal_T c1_A[49], int32_T c1_ilo, int32_T c1_ihi);
static creal_T c1_eml_div(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_x, real_T c1_y);
static void c1_b_eml_matlab_zlartg(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn);
static void c1_b_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static void c1_c_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static void c1_d_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[49]);
static void c1_mldivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_A[49], real_T c1_B[7], creal_T c1_Y[7]);
static creal_T c1_b_eml_div(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_x, creal_T c1_y);
static void c1_d_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static void c1_abs(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   creal_T c1_x[2], real_T c1_y[2]);
static void c1_b_eml_error(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);
static real_T c1_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_ERR, const char_T *c1_identifier);
static real_T c1_b_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_c_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_Wn, const char_T *c1_identifier, real_T
  c1_y[7]);
static void c1_d_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7]);
static void c1_e_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_ABK, const char_T *c1_identifier, real_T
  c1_y[70]);
static void c1_f_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[70]);
static void c1_g_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_CD, const char_T *c1_identifier, real_T
  c1_y[16]);
static void c1_h_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[16]);
static void c1_i_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Glk, const char_T *c1_identifier, real_T
  c1_y[459]);
static void c1_j_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[459]);
static void c1_k_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_LK, const char_T *c1_identifier, real_T
  c1_y[15000]);
static void c1_l_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[15000]);
static void c1_m_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Pabk, const char_T *c1_identifier, real_T
  c1_y[100]);
static void c1_n_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[100]);
static void c1_o_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Pcd, const char_T *c1_identifier, real_T
  c1_y[64]);
static void c1_p_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[64]);
static void c1_q_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Plk, const char_T *c1_identifier, real_T
  c1_y[1099]);
static void c1_r_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[1099]);
static real_T c1_s_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_U1, const char_T *c1_identifier);
static real_T c1_t_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_u_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_VARX, const char_T *c1_identifier, real_T
  c1_y[306]);
static void c1_v_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[306]);
static void c1_w_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_VARX1, const char_T *c1_identifier, real_T
  c1_y[300]);
static void c1_x_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[300]);
static void c1_y_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_X, const char_T *c1_identifier, real_T
  c1_y[7]);
static void c1_ab_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7]);
static void c1_bb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_X1, const char_T *c1_identifier, real_T
  c1_y[7]);
static void c1_cb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7]);
static void c1_db_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Y1, const char_T *c1_identifier, real_T
  c1_y[2]);
static void c1_eb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[2]);
static void c1_fb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Ylk, const char_T *c1_identifier, real_T
  c1_y[3]);
static void c1_gb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3]);
static void c1_hb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Z, const char_T *c1_identifier, real_T
  c1_y[156]);
static void c1_ib_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[156]);
static void c1_jb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_eta, const char_T *c1_identifier, real_T
  c1_y[468]);
static void c1_kb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[468]);
static real_T c1_lb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_k, const char_T *c1_identifier);
static real_T c1_mb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId);
static real_T c1_nb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_saw, const char_T *c1_identifier);
static real_T c1_ob_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId);
static real_T c1_pb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_start, const char_T *c1_identifier);
static real_T c1_qb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_rb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_w, const char_T *c1_identifier, real_T
  c1_y[300]);
static void c1_sb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[300]);
static uint8_T c1_tb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_is_active_c1_sim1_rlti_varx_fast_gust2, const char_T *c1_identifier);
static uint8_T c1_ub_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_sqrt(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T *c1_x);
static real_T c1_fastQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[156], real_T c1_u[3], real_T c1_y[2], real_T
  c1_theta[306], real_T c1_M[1099], real_T c1_lambda, real_T c1_tol, real_T
  c1_b_eta[468], real_T c1_ireg, real_T c1_reg, real_T c1_G[459], real_T c1_Y[3],
  real_T *c1_b_saw);
static void c1_inverseQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[8], real_T c1_y[2], real_T c1_theta[16], real_T
  c1_P[64], real_T c1_lambda);
static void c1_b_inverseQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[10], real_T c1_y[7], real_T c1_theta[70], real_T
  c1_P[100], real_T c1_lambda);
static void c1_b_sqrt(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, creal_T *c1_x);
static void c1_exp(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   creal_T *c1_x);
static void c1_log10(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T *c1_x);
static int32_T c1_div_s32_floor(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, int32_T c1_numerator, int32_T c1_denominator);
static void init_dsm_address_info(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c1_Plk_not_empty = FALSE;
  chartInstance->c1_Glk_not_empty = FALSE;
  chartInstance->c1_Ylk_not_empty = FALSE;
  chartInstance->c1_VARX_not_empty = FALSE;
  chartInstance->c1_VARX1_not_empty = FALSE;
  chartInstance->c1_LK_not_empty = FALSE;
  chartInstance->c1_Pcd_not_empty = FALSE;
  chartInstance->c1_Pabk_not_empty = FALSE;
  chartInstance->c1_ABK_not_empty = FALSE;
  chartInstance->c1_CD_not_empty = FALSE;
  chartInstance->c1_U1_not_empty = FALSE;
  chartInstance->c1_Y1_not_empty = FALSE;
  chartInstance->c1_X_not_empty = FALSE;
  chartInstance->c1_X1_not_empty = FALSE;
  chartInstance->c1_Z_not_empty = FALSE;
  chartInstance->c1_eta_not_empty = FALSE;
  chartInstance->c1_w_not_empty = FALSE;
  chartInstance->c1_k_not_empty = FALSE;
  chartInstance->c1_start_not_empty = FALSE;
  chartInstance->c1_saw_not_empty = FALSE;
  chartInstance->c1_is_active_c1_sim1_rlti_varx_fast_gust2 = 0U;
}

static void initialize_params_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  real_T c1_dv0[700];
  int32_T c1_i0;
  real_T c1_d0;
  real_T c1_dv1[4];
  sf_set_error_prefix_string(
    "Error evaluating data 'W' in the parent workspace.\n");
  sf_mex_import_named("W", sf_mex_get_sfun_param(chartInstance->S, 2, 0), c1_dv0,
                      0, 0, 0U, 1, 0U, 2, 7, 100);
  for (c1_i0 = 0; c1_i0 < 700; c1_i0++) {
    chartInstance->c1_W[c1_i0] = c1_dv0[c1_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Du' in the parent workspace.\n");
  sf_mex_import_named("Du", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c1_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c1_Du = c1_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Dy' in the parent workspace.\n");
  sf_mex_import_named("Dy", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      c1_dv1, 0, 0, 0U, 1, 0U, 2, 2, 2);
  for (c1_i0 = 0; c1_i0 < 4; c1_i0++) {
    chartInstance->c1_Dy[c1_i0] = c1_dv1[c1_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static const mxArray *get_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  real_T c1_u;
  const mxArray *c1_b_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_d_y = NULL;
  int32_T c1_i1;
  real_T c1_d_u[7];
  const mxArray *c1_e_y = NULL;
  real_T c1_e_u;
  const mxArray *c1_f_y = NULL;
  real_T c1_f_u[7];
  const mxArray *c1_g_y = NULL;
  real_T c1_g_u[70];
  const mxArray *c1_h_y = NULL;
  real_T c1_h_u[16];
  const mxArray *c1_i_y = NULL;
  real_T c1_i_u[459];
  const mxArray *c1_j_y = NULL;
  real_T c1_j_u[15000];
  const mxArray *c1_k_y = NULL;
  real_T c1_k_u[100];
  const mxArray *c1_l_y = NULL;
  real_T c1_l_u[64];
  const mxArray *c1_m_y = NULL;
  real_T c1_m_u[1099];
  const mxArray *c1_n_y = NULL;
  real_T c1_n_u;
  const mxArray *c1_o_y = NULL;
  real_T c1_o_u[306];
  const mxArray *c1_p_y = NULL;
  real_T c1_p_u[300];
  const mxArray *c1_q_y = NULL;
  real_T c1_q_u[7];
  const mxArray *c1_r_y = NULL;
  real_T c1_r_u[7];
  const mxArray *c1_s_y = NULL;
  real_T c1_s_u[2];
  const mxArray *c1_t_y = NULL;
  real_T c1_t_u[3];
  const mxArray *c1_u_y = NULL;
  real_T c1_u_u[156];
  const mxArray *c1_v_y = NULL;
  real_T c1_v_u[468];
  const mxArray *c1_w_y = NULL;
  real_T c1_w_u;
  const mxArray *c1_x_y = NULL;
  real_T c1_x_u;
  const mxArray *c1_y_y = NULL;
  real_T c1_y_u;
  const mxArray *c1_ab_y = NULL;
  real_T c1_ab_u[300];
  const mxArray *c1_bb_y = NULL;
  uint8_T c1_bb_u;
  const mxArray *c1_cb_y = NULL;
  real_T *c1_ERR;
  real_T *c1_Fs1;
  real_T *c1_Fs2;
  real_T *c1_Ws;
  real_T (*c1_Zn)[7];
  real_T (*c1_Wn)[7];
  c1_ERR = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c1_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c1_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c1_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellarray(27), FALSE);
  c1_u = *c1_ERR;
  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_b_u = *c1_Fs1;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  c1_c_u = *c1_Fs2;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 2, c1_d_y);
  for (c1_i1 = 0; c1_i1 < 7; c1_i1++) {
    c1_d_u[c1_i1] = (*c1_Wn)[c1_i1];
  }

  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_d_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c1_y, 3, c1_e_y);
  c1_e_u = *c1_Ws;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 4, c1_f_y);
  for (c1_i1 = 0; c1_i1 < 7; c1_i1++) {
    c1_f_u[c1_i1] = (*c1_Zn)[c1_i1];
  }

  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", c1_f_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c1_y, 5, c1_g_y);
  for (c1_i1 = 0; c1_i1 < 70; c1_i1++) {
    c1_g_u[c1_i1] = chartInstance->c1_ABK[c1_i1];
  }

  c1_h_y = NULL;
  if (!chartInstance->c1_ABK_not_empty) {
    sf_mex_assign(&c1_h_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_h_y, sf_mex_create("y", c1_g_u, 0, 0U, 1U, 0U, 2, 7, 10),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 6, c1_h_y);
  for (c1_i1 = 0; c1_i1 < 16; c1_i1++) {
    c1_h_u[c1_i1] = chartInstance->c1_CD[c1_i1];
  }

  c1_i_y = NULL;
  if (!chartInstance->c1_CD_not_empty) {
    sf_mex_assign(&c1_i_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_i_y, sf_mex_create("y", c1_h_u, 0, 0U, 1U, 0U, 2, 2, 8),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 7, c1_i_y);
  for (c1_i1 = 0; c1_i1 < 459; c1_i1++) {
    c1_i_u[c1_i1] = chartInstance->c1_Glk[c1_i1];
  }

  c1_j_y = NULL;
  if (!chartInstance->c1_Glk_not_empty) {
    sf_mex_assign(&c1_j_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_j_y, sf_mex_create("y", c1_i_u, 0, 0U, 1U, 0U, 2, 153, 3),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 8, c1_j_y);
  for (c1_i1 = 0; c1_i1 < 15000; c1_i1++) {
    c1_j_u[c1_i1] = chartInstance->c1_LK[c1_i1];
  }

  c1_k_y = NULL;
  if (!chartInstance->c1_LK_not_empty) {
    sf_mex_assign(&c1_k_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_k_y, sf_mex_create("y", c1_j_u, 0, 0U, 1U, 0U, 2, 100, 150),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 9, c1_k_y);
  for (c1_i1 = 0; c1_i1 < 100; c1_i1++) {
    c1_k_u[c1_i1] = chartInstance->c1_Pabk[c1_i1];
  }

  c1_l_y = NULL;
  if (!chartInstance->c1_Pabk_not_empty) {
    sf_mex_assign(&c1_l_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_l_y, sf_mex_create("y", c1_k_u, 0, 0U, 1U, 0U, 2, 10, 10),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 10, c1_l_y);
  for (c1_i1 = 0; c1_i1 < 64; c1_i1++) {
    c1_l_u[c1_i1] = chartInstance->c1_Pcd[c1_i1];
  }

  c1_m_y = NULL;
  if (!chartInstance->c1_Pcd_not_empty) {
    sf_mex_assign(&c1_m_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_m_y, sf_mex_create("y", c1_l_u, 0, 0U, 1U, 0U, 2, 8, 8),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 11, c1_m_y);
  for (c1_i1 = 0; c1_i1 < 1099; c1_i1++) {
    c1_m_u[c1_i1] = chartInstance->c1_Plk[c1_i1];
  }

  c1_n_y = NULL;
  if (!chartInstance->c1_Plk_not_empty) {
    sf_mex_assign(&c1_n_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_n_y, sf_mex_create("y", c1_m_u, 0, 0U, 1U, 0U, 2, 157, 7),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 12, c1_n_y);
  c1_n_u = chartInstance->c1_U1;
  c1_o_y = NULL;
  if (!chartInstance->c1_U1_not_empty) {
    sf_mex_assign(&c1_o_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_o_y, sf_mex_create("y", &c1_n_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c1_y, 13, c1_o_y);
  for (c1_i1 = 0; c1_i1 < 306; c1_i1++) {
    c1_o_u[c1_i1] = chartInstance->c1_VARX[c1_i1];
  }

  c1_p_y = NULL;
  if (!chartInstance->c1_VARX_not_empty) {
    sf_mex_assign(&c1_p_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_p_y, sf_mex_create("y", c1_o_u, 0, 0U, 1U, 0U, 2, 2, 153),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 14, c1_p_y);
  for (c1_i1 = 0; c1_i1 < 300; c1_i1++) {
    c1_p_u[c1_i1] = chartInstance->c1_VARX1[c1_i1];
  }

  c1_q_y = NULL;
  if (!chartInstance->c1_VARX1_not_empty) {
    sf_mex_assign(&c1_q_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_q_y, sf_mex_create("y", c1_p_u, 0, 0U, 1U, 0U, 2, 2, 150),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 15, c1_q_y);
  for (c1_i1 = 0; c1_i1 < 7; c1_i1++) {
    c1_q_u[c1_i1] = chartInstance->c1_X[c1_i1];
  }

  c1_r_y = NULL;
  if (!chartInstance->c1_X_not_empty) {
    sf_mex_assign(&c1_r_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_r_y, sf_mex_create("y", c1_q_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 16, c1_r_y);
  for (c1_i1 = 0; c1_i1 < 7; c1_i1++) {
    c1_r_u[c1_i1] = chartInstance->c1_X1[c1_i1];
  }

  c1_s_y = NULL;
  if (!chartInstance->c1_X1_not_empty) {
    sf_mex_assign(&c1_s_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_s_y, sf_mex_create("y", c1_r_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 17, c1_s_y);
  for (c1_i1 = 0; c1_i1 < 2; c1_i1++) {
    c1_s_u[c1_i1] = chartInstance->c1_Y1[c1_i1];
  }

  c1_t_y = NULL;
  if (!chartInstance->c1_Y1_not_empty) {
    sf_mex_assign(&c1_t_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_t_y, sf_mex_create("y", c1_s_u, 0, 0U, 1U, 0U, 1, 2),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 18, c1_t_y);
  for (c1_i1 = 0; c1_i1 < 3; c1_i1++) {
    c1_t_u[c1_i1] = chartInstance->c1_Ylk[c1_i1];
  }

  c1_u_y = NULL;
  if (!chartInstance->c1_Ylk_not_empty) {
    sf_mex_assign(&c1_u_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_u_y, sf_mex_create("y", c1_t_u, 0, 0U, 1U, 0U, 2, 1, 3),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 19, c1_u_y);
  for (c1_i1 = 0; c1_i1 < 156; c1_i1++) {
    c1_u_u[c1_i1] = chartInstance->c1_Z[c1_i1];
  }

  c1_v_y = NULL;
  if (!chartInstance->c1_Z_not_empty) {
    sf_mex_assign(&c1_v_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_v_y, sf_mex_create("y", c1_u_u, 0, 0U, 1U, 0U, 1, 156),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 20, c1_v_y);
  for (c1_i1 = 0; c1_i1 < 468; c1_i1++) {
    c1_v_u[c1_i1] = chartInstance->c1_eta[c1_i1];
  }

  c1_w_y = NULL;
  if (!chartInstance->c1_eta_not_empty) {
    sf_mex_assign(&c1_w_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_w_y, sf_mex_create("y", c1_v_u, 0, 0U, 1U, 0U, 2, 156, 3),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 21, c1_w_y);
  c1_w_u = chartInstance->c1_k;
  c1_x_y = NULL;
  if (!chartInstance->c1_k_not_empty) {
    sf_mex_assign(&c1_x_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_x_y, sf_mex_create("y", &c1_w_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c1_y, 22, c1_x_y);
  c1_x_u = chartInstance->c1_saw;
  c1_y_y = NULL;
  if (!chartInstance->c1_saw_not_empty) {
    sf_mex_assign(&c1_y_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_y_y, sf_mex_create("y", &c1_x_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c1_y, 23, c1_y_y);
  c1_y_u = chartInstance->c1_start;
  c1_ab_y = NULL;
  if (!chartInstance->c1_start_not_empty) {
    sf_mex_assign(&c1_ab_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_ab_y, sf_mex_create("y", &c1_y_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c1_y, 24, c1_ab_y);
  for (c1_i1 = 0; c1_i1 < 300; c1_i1++) {
    c1_ab_u[c1_i1] = chartInstance->c1_w[c1_i1];
  }

  c1_bb_y = NULL;
  if (!chartInstance->c1_w_not_empty) {
    sf_mex_assign(&c1_bb_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c1_bb_y, sf_mex_create("y", c1_ab_u, 0, 0U, 1U, 0U, 2, 1, 300),
                  FALSE);
  }

  sf_mex_setcell(c1_y, 25, c1_bb_y);
  c1_bb_u = chartInstance->c1_is_active_c1_sim1_rlti_varx_fast_gust2;
  c1_cb_y = NULL;
  sf_mex_assign(&c1_cb_y, sf_mex_create("y", &c1_bb_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 26, c1_cb_y);
  sf_mex_assign(&c1_st, c1_y, FALSE);
  return c1_st;
}

static void set_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv2[7];
  int32_T c1_i2;
  real_T c1_dv3[70];
  real_T c1_dv4[16];
  real_T c1_dv5[459];
  static real_T c1_dv6[15000];
  real_T c1_dv7[100];
  real_T c1_dv8[64];
  real_T c1_dv9[1099];
  real_T c1_dv10[306];
  real_T c1_dv11[300];
  real_T c1_dv12[2];
  real_T c1_dv13[3];
  real_T c1_dv14[156];
  real_T c1_dv15[468];
  real_T *c1_ERR;
  real_T *c1_Fs1;
  real_T *c1_Fs2;
  real_T *c1_Ws;
  real_T (*c1_Wn)[7];
  real_T (*c1_Zn)[7];
  c1_ERR = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c1_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c1_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c1_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_u = sf_mex_dup(c1_st);
  *c1_ERR = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)),
    "ERR");
  *c1_Fs1 = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 1)),
    "Fs1");
  *c1_Fs2 = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 2)),
    "Fs2");
  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 3)), "Wn",
                        c1_dv2);
  for (c1_i2 = 0; c1_i2 < 7; c1_i2++) {
    (*c1_Wn)[c1_i2] = c1_dv2[c1_i2];
  }

  *c1_Ws = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 4)),
    "Ws");
  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 5)), "Zn",
                        c1_dv2);
  for (c1_i2 = 0; c1_i2 < 7; c1_i2++) {
    (*c1_Zn)[c1_i2] = c1_dv2[c1_i2];
  }

  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 6)),
                        "ABK", c1_dv3);
  for (c1_i2 = 0; c1_i2 < 70; c1_i2++) {
    chartInstance->c1_ABK[c1_i2] = c1_dv3[c1_i2];
  }

  c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 7)), "CD",
                        c1_dv4);
  for (c1_i2 = 0; c1_i2 < 16; c1_i2++) {
    chartInstance->c1_CD[c1_i2] = c1_dv4[c1_i2];
  }

  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 8)),
                        "Glk", c1_dv5);
  for (c1_i2 = 0; c1_i2 < 459; c1_i2++) {
    chartInstance->c1_Glk[c1_i2] = c1_dv5[c1_i2];
  }

  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 9)), "LK",
                        c1_dv6);
  for (c1_i2 = 0; c1_i2 < 15000; c1_i2++) {
    chartInstance->c1_LK[c1_i2] = c1_dv6[c1_i2];
  }

  c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 10)),
                        "Pabk", c1_dv7);
  for (c1_i2 = 0; c1_i2 < 100; c1_i2++) {
    chartInstance->c1_Pabk[c1_i2] = c1_dv7[c1_i2];
  }

  c1_o_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 11)),
                        "Pcd", c1_dv8);
  for (c1_i2 = 0; c1_i2 < 64; c1_i2++) {
    chartInstance->c1_Pcd[c1_i2] = c1_dv8[c1_i2];
  }

  c1_q_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 12)),
                        "Plk", c1_dv9);
  for (c1_i2 = 0; c1_i2 < 1099; c1_i2++) {
    chartInstance->c1_Plk[c1_i2] = c1_dv9[c1_i2];
  }

  chartInstance->c1_U1 = c1_s_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_u, 13)), "U1");
  c1_u_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 14)),
                        "VARX", c1_dv10);
  for (c1_i2 = 0; c1_i2 < 306; c1_i2++) {
    chartInstance->c1_VARX[c1_i2] = c1_dv10[c1_i2];
  }

  c1_w_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 15)),
                        "VARX1", c1_dv11);
  for (c1_i2 = 0; c1_i2 < 300; c1_i2++) {
    chartInstance->c1_VARX1[c1_i2] = c1_dv11[c1_i2];
  }

  c1_y_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 16)), "X",
                        c1_dv2);
  for (c1_i2 = 0; c1_i2 < 7; c1_i2++) {
    chartInstance->c1_X[c1_i2] = c1_dv2[c1_i2];
  }

  c1_bb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 17)),
    "X1", c1_dv2);
  for (c1_i2 = 0; c1_i2 < 7; c1_i2++) {
    chartInstance->c1_X1[c1_i2] = c1_dv2[c1_i2];
  }

  c1_db_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 18)),
    "Y1", c1_dv12);
  for (c1_i2 = 0; c1_i2 < 2; c1_i2++) {
    chartInstance->c1_Y1[c1_i2] = c1_dv12[c1_i2];
  }

  c1_fb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 19)),
    "Ylk", c1_dv13);
  for (c1_i2 = 0; c1_i2 < 3; c1_i2++) {
    chartInstance->c1_Ylk[c1_i2] = c1_dv13[c1_i2];
  }

  c1_hb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 20)),
    "Z", c1_dv14);
  for (c1_i2 = 0; c1_i2 < 156; c1_i2++) {
    chartInstance->c1_Z[c1_i2] = c1_dv14[c1_i2];
  }

  c1_jb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 21)),
    "eta", c1_dv15);
  for (c1_i2 = 0; c1_i2 < 468; c1_i2++) {
    chartInstance->c1_eta[c1_i2] = c1_dv15[c1_i2];
  }

  chartInstance->c1_k = c1_lb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_u, 22)), "k");
  chartInstance->c1_saw = c1_nb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_u, 23)), "saw");
  chartInstance->c1_start = c1_pb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c1_u, 24)), "start");
  c1_rb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 25)),
    "w", c1_dv11);
  for (c1_i2 = 0; c1_i2 < 300; c1_i2++) {
    chartInstance->c1_w[c1_i2] = c1_dv11[c1_i2];
  }

  chartInstance->c1_is_active_c1_sim1_rlti_varx_fast_gust2 =
    c1_tb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 26)),
    "is_active_c1_sim1_rlti_varx_fast_gust2");
  sf_mex_destroy(&c1_u);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
}

static void sf_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  c1_chartstep_c1_sim1_rlti_varx_fast_gust2(chartInstance);
}

static void c1_chartstep_c1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
  real_T c1_ON;
  real_T c1_U;
  int32_T c1_i3;
  real_T c1_Fs[2];
  real_T c1_b_W[700];
  real_T c1_b_Du;
  real_T c1_b_Dy[4];
  real_T c1_beta1;
  int32_T c1_ldc;
  int32_T c1_i;
  real_T c1_dv16[300];
  real_T c1_c_Dy[2];
  real_T c1_dv17[2];
  real_T c1_a[4];
  real_T c1_b_VARX[306];
  real_T c1_b_Z[156];
  real_T c1_b_Plk[1099];
  real_T c1_b_eta[468];
  real_T c1_b_Glk[459];
  real_T c1_b_Ylk[3];
  real_T c1_dv18[3];
  real_T c1_ERR;
  int32_T c1_n;
  int32_T c1_b_k;
  int32_T c1_lda;
  int32_T c1_ldb;
  char_T c1_TRANSA;
  char_T c1_TRANSB;
  real_T c1_y[1050];
  real_T c1_b[150];
  real_T c1_b_CD[16];
  real_T c1_b_Pcd[64];
  real_T c1_dv19[8];
  real_T c1_b_ABK[70];
  real_T c1_b_Pabk[100];
  real_T c1_dv20[10];
  real_T c1_b_a[7];
  real_T c1_A[49];
  real_T c1_C[14];
  real_T c1_dv21[49];
  real_T c1_Zn[7];
  real_T c1_Wn[7];
  creal_T c1_b_y;
  creal_T c1_c_y[49];
  real_T c1_c_a[7];
  creal_T c1_b_b[7];
  creal_T c1_b_C[14];
  creal_T c1_c_C[2];
  real_T *c1_Ws;
  real_T *c1_Fs1;
  real_T *c1_b_ERR;
  real_T *c1_Fs2;
  real_T *c1_b_ON;
  real_T *c1_b_U;
  real_T *c1_RESET;
  real_T (*c1_b_Wn)[7];
  real_T (*c1_b_Zn)[7];
  real_T (*c1_Y)[2];
  c1_b_ERR = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c1_b_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c1_b_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c1_Y = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c1_b_U = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_RESET = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_ON = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c1_ON = *c1_b_ON;
  c1_U = *c1_b_U;
  for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
    c1_Fs[c1_i3] = (*c1_Y)[c1_i3];
  }

  for (c1_i3 = 0; c1_i3 < 700; c1_i3++) {
    c1_b_W[c1_i3] = chartInstance->c1_W[c1_i3];
  }

  c1_b_Du = chartInstance->c1_Du;
  for (c1_i3 = 0; c1_i3 < 4; c1_i3++) {
    c1_b_Dy[c1_i3] = chartInstance->c1_Dy[c1_i3];
  }

  if ((!chartInstance->c1_Plk_not_empty) || (*c1_RESET != 0.0)) {
    for (c1_i3 = 0; c1_i3 < 1099; c1_i3++) {
      chartInstance->c1_Plk[c1_i3] = 0.0;
    }

    chartInstance->c1_Plk_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 459; c1_i3++) {
      chartInstance->c1_Glk[c1_i3] = 0.0;
    }

    chartInstance->c1_Glk_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
      chartInstance->c1_Ylk[c1_i3] = 1.0;
    }

    chartInstance->c1_Ylk_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 306; c1_i3++) {
      chartInstance->c1_VARX[c1_i3] = 0.0;
    }

    chartInstance->c1_VARX_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 300; c1_i3++) {
      chartInstance->c1_VARX1[c1_i3] = 0.0;
    }

    chartInstance->c1_VARX1_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 156; c1_i3++) {
      chartInstance->c1_Z[c1_i3] = 0.0;
    }

    chartInstance->c1_Z_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 468; c1_i3++) {
      chartInstance->c1_eta[c1_i3] = 0.0;
    }

    chartInstance->c1_eta_not_empty = TRUE;
    chartInstance->c1_saw = 0.0;
    chartInstance->c1_saw_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 15000; c1_i3++) {
      chartInstance->c1_LK[c1_i3] = 0.0;
    }

    chartInstance->c1_LK_not_empty = TRUE;
    c1_eye(chartInstance, chartInstance->c1_Pabk);
    c1_beta1 = c1_mrdivide(chartInstance, 1.0, 0.001);
    for (c1_i3 = 0; c1_i3 < 100; c1_i3++) {
      chartInstance->c1_Pabk[c1_i3] *= c1_beta1;
    }

    chartInstance->c1_Pabk_not_empty = TRUE;
    c1_i3 = 0;
    for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
      for (c1_i = 0; c1_i < 7; c1_i++) {
        chartInstance->c1_ABK[c1_i + c1_i3] = 0.0;
      }

      c1_i3 += 7;
    }

    for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
      chartInstance->c1_ABK[c1_i3 + 49] = 0.0;
    }

    c1_i3 = 0;
    for (c1_ldc = 0; c1_ldc < 2; c1_ldc++) {
      for (c1_i = 0; c1_i < 7; c1_i++) {
        chartInstance->c1_ABK[(c1_i + c1_i3) + 56] = 0.0;
      }

      c1_i3 += 7;
    }

    chartInstance->c1_ABK_not_empty = TRUE;
    c1_b_eye(chartInstance, chartInstance->c1_Pcd);
    c1_beta1 = c1_mrdivide(chartInstance, 1.0, 0.001);
    for (c1_i3 = 0; c1_i3 < 64; c1_i3++) {
      chartInstance->c1_Pcd[c1_i3] *= c1_beta1;
    }

    chartInstance->c1_Pcd_not_empty = TRUE;
    c1_i3 = 0;
    for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
      for (c1_i = 0; c1_i < 2; c1_i++) {
        chartInstance->c1_CD[c1_i + c1_i3] = 0.0;
      }

      c1_i3 += 2;
    }

    chartInstance->c1_CD_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
      chartInstance->c1_CD[c1_i3 + 14] = 0.0;
      chartInstance->c1_Y1[c1_i3] = 0.0;
    }

    chartInstance->c1_Y1_not_empty = TRUE;
    chartInstance->c1_U1 = 0.0;
    chartInstance->c1_U1_not_empty = TRUE;
    chartInstance->c1_X1_not_empty = TRUE;
    for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
      chartInstance->c1_X1[c1_i3] = 0.0;
      chartInstance->c1_X[c1_i3] = 0.0;
    }

    chartInstance->c1_X_not_empty = TRUE;
    c1_logspace(chartInstance, -1.0, 1.0, c1_dv16);
    for (c1_i3 = 0; c1_i3 < 300; c1_i3++) {
      chartInstance->c1_w[c1_i3] = c1_dv16[c1_i3];
    }

    chartInstance->c1_w_not_empty = TRUE;
    chartInstance->c1_k = 1.0;
    chartInstance->c1_k_not_empty = TRUE;
    chartInstance->c1_start = 1.0;
    chartInstance->c1_start_not_empty = TRUE;
  }

  if (c1_ON > 0.5) {
    c1_U *= c1_rdivide(chartInstance, 1.0, c1_b_Du);
    c1_i3 = 0;
    for (c1_ldc = 0; c1_ldc < 2; c1_ldc++) {
      c1_c_Dy[c1_ldc] = c1_b_Dy[c1_i3];
      c1_i3 += 3;
    }

    c1_b_rdivide(chartInstance, 1.0, c1_c_Dy, c1_dv17);
    c1_diag(chartInstance, c1_dv17, c1_a);
    for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
      c1_c_Dy[c1_i3] = 0.0;
      c1_ldc = 0;
      for (c1_i = 0; c1_i < 2; c1_i++) {
        c1_c_Dy[c1_i3] += c1_a[c1_ldc + c1_i3] * c1_Fs[c1_i];
        c1_ldc += 2;
      }
    }

    for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
      c1_Fs[c1_i3] = c1_c_Dy[c1_i3];
    }

    if ((chartInstance->c1_start < 0.5) || (chartInstance->c1_k >= 50.0)) {
      for (c1_i3 = 0; c1_i3 < 306; c1_i3++) {
        c1_b_VARX[c1_i3] = chartInstance->c1_VARX[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 156; c1_i3++) {
        c1_b_Z[c1_i3] = chartInstance->c1_Z[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 1099; c1_i3++) {
        c1_b_Plk[c1_i3] = chartInstance->c1_Plk[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 468; c1_i3++) {
        c1_b_eta[c1_i3] = chartInstance->c1_eta[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 459; c1_i3++) {
        c1_b_Glk[c1_i3] = chartInstance->c1_Glk[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
        c1_b_Ylk[c1_i3] = chartInstance->c1_Ylk[c1_i3];
      }

      c1_ON = chartInstance->c1_saw;
      for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
        c1_dv18[c1_i3] = chartInstance->c1_Y1[c1_i3];
      }

      c1_dv18[2] = c1_U;
      c1_ERR = c1_fastQR(chartInstance, c1_b_Z, c1_dv18, c1_Fs, c1_b_VARX,
                         c1_b_Plk, 0.999, 0.0001, c1_b_eta, 100.0, 1.0, c1_b_Glk,
                         c1_b_Ylk, &c1_ON);
      for (c1_i3 = 0; c1_i3 < 306; c1_i3++) {
        chartInstance->c1_VARX[c1_i3] = c1_b_VARX[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 156; c1_i3++) {
        chartInstance->c1_Z[c1_i3] = c1_b_Z[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 1099; c1_i3++) {
        chartInstance->c1_Plk[c1_i3] = c1_b_Plk[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 468; c1_i3++) {
        chartInstance->c1_eta[c1_i3] = c1_b_eta[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 459; c1_i3++) {
        chartInstance->c1_Glk[c1_i3] = c1_b_Glk[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
        chartInstance->c1_Ylk[c1_i3] = c1_b_Ylk[c1_i3];
      }

      chartInstance->c1_saw = c1_ON;
      c1_i3 = 0;
      for (c1_ldc = 0; c1_ldc < 150; c1_ldc++) {
        for (c1_i = 0; c1_i < 2; c1_i++) {
          chartInstance->c1_VARX1[c1_i + c1_i3] = chartInstance->c1_VARX[(c1_i +
            c1_i3) + 4];
        }

        c1_i3 += 2;
      }
    } else {
      c1_ERR = 0.0;
    }

    if ((chartInstance->c1_start < 0.5) || (chartInstance->c1_k >= 150.0)) {
      for (c1_i3 = 0; c1_i3 < 15000; c1_i3++) {
        chartInstance->c1_LK[c1_i3] = 0.0;
      }

      for (c1_i = 0; c1_i < 50; c1_i++) {
        c1_n = 0;
        while (c1_n <= 49 - c1_i) {
          c1_b_k = (c1_i << 1) - 1;
          c1_lda = (c1_i + c1_n) * 3;
          c1_ldb = c1_n * 3 - 1;
          for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
            for (c1_ldc = 0; c1_ldc < 2; c1_ldc++) {
              chartInstance->c1_LK[((c1_ldc + c1_b_k) + 100 *
                                    (sf_mex_lw_bounds_check((c1_i3 + c1_lda) + 1,
                1, 150) - 1)) + 1] = chartInstance->c1_VARX1[c1_ldc + (((c1_i3 +
                c1_ldb) + 1) << 1)];
            }
          }

          c1_n++;
          sf_mex_listen_for_ctrl_c(chartInstance->S);
        }

        sf_mex_listen_for_ctrl_c(chartInstance->S);
      }

      c1_i = 7;
      c1_n = 150;
      c1_b_k = 100;
      c1_ON = 1.0;
      c1_lda = 7;
      c1_ldb = 100;
      c1_beta1 = 0.0;
      c1_ldc = 7;
      c1_TRANSA = 'N';
      c1_TRANSB = 'N';
      for (c1_i3 = 0; c1_i3 < 1050; c1_i3++) {
        c1_y[c1_i3] = 0.0;
      }

      dgemm32(&c1_TRANSA, &c1_TRANSB, &c1_i, &c1_n, &c1_b_k, &c1_ON, &c1_b_W[0],
              &c1_lda, &chartInstance->c1_LK[0], &c1_ldb, &c1_beta1, &c1_y[0],
              &c1_ldc);
      for (c1_i3 = 0; c1_i3 < 150; c1_i3++) {
        c1_b[c1_i3] = chartInstance->c1_Z[c1_i3 + 5];
      }

      c1_i = 7;
      c1_n = 1;
      c1_b_k = 150;
      c1_ON = 1.0;
      c1_lda = 7;
      c1_ldb = 150;
      c1_beta1 = 0.0;
      c1_ldc = 7;
      c1_TRANSA = 'N';
      c1_TRANSB = 'N';
      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        chartInstance->c1_X[c1_i3] = 0.0;
      }

      dgemm32(&c1_TRANSA, &c1_TRANSB, &c1_i, &c1_n, &c1_b_k, &c1_ON, &c1_y[0],
              &c1_lda, &c1_b[0], &c1_ldb, &c1_beta1, &chartInstance->c1_X[0],
              &c1_ldc);
    }

    if ((chartInstance->c1_start < 0.5) || (chartInstance->c1_k >= 151.0)) {
      for (c1_i3 = 0; c1_i3 < 16; c1_i3++) {
        c1_b_CD[c1_i3] = chartInstance->c1_CD[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 64; c1_i3++) {
        c1_b_Pcd[c1_i3] = chartInstance->c1_Pcd[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        c1_dv19[c1_i3] = chartInstance->c1_X1[c1_i3];
      }

      c1_dv19[7] = chartInstance->c1_U1;
      for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
        c1_c_Dy[c1_i3] = chartInstance->c1_Y1[c1_i3];
      }

      c1_inverseQR(chartInstance, c1_dv19, c1_c_Dy, c1_b_CD, c1_b_Pcd, 0.999);
      for (c1_i3 = 0; c1_i3 < 16; c1_i3++) {
        chartInstance->c1_CD[c1_i3] = c1_b_CD[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 64; c1_i3++) {
        chartInstance->c1_Pcd[c1_i3] = c1_b_Pcd[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 70; c1_i3++) {
        c1_b_ABK[c1_i3] = chartInstance->c1_ABK[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 100; c1_i3++) {
        c1_b_Pabk[c1_i3] = chartInstance->c1_Pabk[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        c1_dv19[c1_i3] = chartInstance->c1_X1[c1_i3];
      }

      c1_dv19[7] = chartInstance->c1_U1;
      for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
        c1_beta1 = 0.0;
        c1_ldc = 0;
        for (c1_i = 0; c1_i < 8; c1_i++) {
          c1_beta1 += chartInstance->c1_CD[c1_ldc + c1_i3] * c1_dv19[c1_i];
          c1_ldc += 2;
        }

        c1_c_Dy[c1_i3] = chartInstance->c1_Y1[c1_i3] - c1_beta1;
      }

      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        c1_dv20[c1_i3] = chartInstance->c1_X1[c1_i3];
      }

      c1_dv20[7] = chartInstance->c1_U1;
      for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
        c1_dv20[c1_i3 + 8] = c1_c_Dy[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        c1_b_a[c1_i3] = chartInstance->c1_X[c1_i3];
      }

      c1_b_inverseQR(chartInstance, c1_dv20, c1_b_a, c1_b_ABK, c1_b_Pabk, 0.999);
      for (c1_i3 = 0; c1_i3 < 70; c1_i3++) {
        chartInstance->c1_ABK[c1_i3] = c1_b_ABK[c1_i3];
      }

      for (c1_i3 = 0; c1_i3 < 100; c1_i3++) {
        chartInstance->c1_Pabk[c1_i3] = c1_b_Pabk[c1_i3];
      }
    }

    if ((chartInstance->c1_start < 0.5) || (chartInstance->c1_k >= 150.0)) {
      for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
        chartInstance->c1_X1[c1_i3] = chartInstance->c1_X[c1_i3];
      }
    }

    for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
      chartInstance->c1_Y1[c1_i3] = c1_Fs[c1_i3];
    }

    chartInstance->c1_U1 = c1_U;
  } else {
    c1_ERR = 0.0;
  }

  c1_ON = c1_rdivide(chartInstance, 1.0, c1_b_Du);
  c1_i3 = 0;
  for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_A[c1_i + c1_i3] = chartInstance->c1_ABK[c1_i + c1_i3];
    }

    c1_b_a[c1_ldc] = chartInstance->c1_ABK[c1_ldc + 49];
    c1_i3 += 7;
  }

  for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
    for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
      c1_C[c1_i3 + (c1_ldc << 1)] = 0.0;
      for (c1_i = 0; c1_i < 2; c1_i++) {
        c1_C[c1_i3 + (c1_ldc << 1)] += c1_b_Dy[c1_i3 + (c1_i << 1)] *
          chartInstance->c1_CD[c1_i + (c1_ldc << 1)];
      }
    }

    c1_Fs[c1_i3] = 0.0;
    for (c1_ldc = 0; c1_ldc < 2; c1_ldc++) {
      c1_Fs[c1_i3] += c1_b_Dy[c1_i3 + (c1_ldc << 1)] * chartInstance->c1_CD[14 +
        c1_ldc];
    }
  }

  c1_U = c1_rdivide(chartInstance, 1.0, c1_b_Du);
  c1_i3 = 0;
  for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_dv21[c1_i + c1_i3] = chartInstance->c1_ABK[c1_i + c1_i3];
    }

    c1_i3 += 7;
  }

  c1_damp(chartInstance, c1_dv21, 0.05, c1_Wn, c1_Zn);
  c1_b_Du = chartInstance->c1_w[sf_mex_lw_bounds_check((int32_T)
    chartInstance->c1_k, 1, 300) - 1];
  c1_b_y.re = 3.1415926535897931 * (2.0 * (chartInstance->c1_w[(int32_T)
    chartInstance->c1_k - 1] * 0.0));
  c1_b_y.im = 3.1415926535897931 * (2.0 * (chartInstance->c1_w[(int32_T)
    chartInstance->c1_k - 1] * 0.05));
  c1_exp(chartInstance, &c1_b_y);
  c1_d_eye(chartInstance, c1_dv21);
  for (c1_i3 = 0; c1_i3 < 49; c1_i3++) {
    c1_c_y[c1_i3].re = c1_b_y.re * c1_dv21[c1_i3] - c1_A[c1_i3];
    c1_c_y[c1_i3].im = c1_b_y.im * c1_dv21[c1_i3];
  }

  for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
    c1_c_a[c1_i3] = c1_b_a[c1_i3] * c1_ON;
  }

  c1_mldivide(chartInstance, c1_c_y, c1_c_a, c1_b_b);
  c1_i3 = 0;
  for (c1_ldc = 0; c1_ldc < 7; c1_ldc++) {
    for (c1_i = 0; c1_i < 2; c1_i++) {
      c1_b_C[c1_i + c1_i3].re = c1_C[c1_i + c1_i3];
      c1_b_C[c1_i + c1_i3].im = 0.0;
    }

    c1_i3 += 2;
  }

  for (c1_i3 = 0; c1_i3 < 2; c1_i3++) {
    c1_ON = 0.0;
    c1_beta1 = 0.0;
    c1_ldc = 0;
    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_ON += c1_b_C[c1_ldc + c1_i3].re * c1_b_b[c1_i].re - 0.0 * c1_b_b[c1_i].
        im;
      c1_beta1 += c1_b_C[c1_ldc + c1_i3].re * c1_b_b[c1_i].im + 0.0 *
        c1_b_b[c1_i].re;
      c1_ldc += 2;
    }

    c1_c_C[c1_i3].re = c1_ON + c1_Fs[c1_i3] * c1_U;
    c1_c_C[c1_i3].im = c1_beta1;
  }

  c1_abs(chartInstance, c1_c_C, c1_Fs);
  c1_log10(chartInstance, &c1_b_Du);
  if (c1_Fs[0] < 1.0E-5) {
    c1_beta1 = -100.0;
  } else {
    c1_ON = c1_Fs[0];
    c1_log10(chartInstance, &c1_ON);
    c1_beta1 = 20.0 * c1_ON;
  }

  if (c1_Fs[1] < 1.0E-5) {
    *c1_Fs2 = -100.0;
  } else {
    c1_ON = c1_Fs[1];
    c1_log10(chartInstance, &c1_ON);
    *c1_Fs2 = 20.0 * c1_ON;
  }

  if (chartInstance->c1_k >= 300.0) {
    chartInstance->c1_start = 0.0;
    chartInstance->c1_k = 1.0;
  } else {
    chartInstance->c1_k++;
  }

  *c1_Ws = c1_b_Du;
  *c1_Fs1 = c1_beta1;
  for (c1_i3 = 0; c1_i3 < 7; c1_i3++) {
    (*c1_b_Wn)[c1_i3] = c1_Wn[c1_i3];
    (*c1_b_Zn)[c1_i3] = c1_Zn[c1_i3];
  }

  *c1_b_ERR = c1_ERR;
}

static void initSimStructsc1_sim1_rlti_varx_fast_gust2
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber)
{
}

const mxArray *sf_c1_sim1_rlti_varx_fast_gust2_get_eml_resolved_functions_info
  (void)
{
  const mxArray *c1_nameCaptureInfo;
  c1_ResolvedFunctionInfo c1_info[102];
  const mxArray *c1_m0 = NULL;
  int32_T c1_i4;
  c1_ResolvedFunctionInfo *c1_r0;
  c1_nameCaptureInfo = NULL;
  c1_info_helper(c1_info);
  c1_b_info_helper(c1_info);
  sf_mex_assign(&c1_m0, sf_mex_createstruct("nameCaptureInfo", 1, 102), FALSE);
  for (c1_i4 = 0; c1_i4 < 102; c1_i4++) {
    c1_r0 = &c1_info[c1_i4];
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->context)), "context", "nameCaptureInfo",
                    c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c1_r0->name)), "name", "nameCaptureInfo", c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c1_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->resolved)), "resolved", "nameCaptureInfo",
                    c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c1_i4);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c1_i4);
  }

  sf_mex_assign(&c1_nameCaptureInfo, c1_m0, FALSE);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(c1_ResolvedFunctionInfo c1_info[102])
{
  c1_info[0].context = "";
  c1_info[0].name = "mtimes";
  c1_info[0].dominantType = "double";
  c1_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[0].fileTimeLo = 1289519692U;
  c1_info[0].fileTimeHi = 0U;
  c1_info[0].mFileTimeLo = 0U;
  c1_info[0].mFileTimeHi = 0U;
  c1_info[1].context = "";
  c1_info[1].name = "sqrt";
  c1_info[1].dominantType = "double";
  c1_info[1].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[1].fileTimeLo = 1286818752U;
  c1_info[1].fileTimeHi = 0U;
  c1_info[1].mFileTimeLo = 0U;
  c1_info[1].mFileTimeHi = 0U;
  c1_info[2].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[2].name = "eml_error";
  c1_info[2].dominantType = "char";
  c1_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c1_info[2].fileTimeLo = 1305318000U;
  c1_info[2].fileTimeHi = 0U;
  c1_info[2].mFileTimeLo = 0U;
  c1_info[2].mFileTimeHi = 0U;
  c1_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[3].name = "eml_scalar_sqrt";
  c1_info[3].dominantType = "double";
  c1_info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c1_info[3].fileTimeLo = 1286818738U;
  c1_info[3].fileTimeHi = 0U;
  c1_info[3].mFileTimeLo = 0U;
  c1_info[3].mFileTimeHi = 0U;
  c1_info[4].context = "";
  c1_info[4].name = "mrdivide";
  c1_info[4].dominantType = "double";
  c1_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c1_info[4].fileTimeLo = 1310137456U;
  c1_info[4].fileTimeHi = 0U;
  c1_info[4].mFileTimeLo = 1289519692U;
  c1_info[4].mFileTimeHi = 0U;
  c1_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c1_info[5].name = "rdivide";
  c1_info[5].dominantType = "double";
  c1_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c1_info[5].fileTimeLo = 1286818844U;
  c1_info[5].fileTimeHi = 0U;
  c1_info[5].mFileTimeLo = 0U;
  c1_info[5].mFileTimeHi = 0U;
  c1_info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c1_info[6].name = "eml_div";
  c1_info[6].dominantType = "double";
  c1_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c1_info[6].fileTimeLo = 1305318000U;
  c1_info[6].fileTimeHi = 0U;
  c1_info[6].mFileTimeLo = 0U;
  c1_info[6].mFileTimeHi = 0U;
  c1_info[7].context = "";
  c1_info[7].name = "eye";
  c1_info[7].dominantType = "double";
  c1_info[7].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m";
  c1_info[7].fileTimeLo = 1286818688U;
  c1_info[7].fileTimeHi = 0U;
  c1_info[7].mFileTimeLo = 0U;
  c1_info[7].mFileTimeHi = 0U;
  c1_info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c1_info[8].name = "eml_assert_valid_size_arg";
  c1_info[8].dominantType = "double";
  c1_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c1_info[8].fileTimeLo = 1286818694U;
  c1_info[8].fileTimeHi = 0U;
  c1_info[8].mFileTimeLo = 0U;
  c1_info[8].mFileTimeHi = 0U;
  c1_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral";
  c1_info[9].name = "isinf";
  c1_info[9].dominantType = "double";
  c1_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m";
  c1_info[9].fileTimeLo = 1286818760U;
  c1_info[9].fileTimeHi = 0U;
  c1_info[9].mFileTimeLo = 0U;
  c1_info[9].mFileTimeHi = 0U;
  c1_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c1_info[10].name = "eml_index_class";
  c1_info[10].dominantType = "";
  c1_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[10].fileTimeLo = 1286818778U;
  c1_info[10].fileTimeHi = 0U;
  c1_info[10].mFileTimeLo = 0U;
  c1_info[10].mFileTimeHi = 0U;
  c1_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c1_info[11].name = "intmax";
  c1_info[11].dominantType = "char";
  c1_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[11].fileTimeLo = 1286818756U;
  c1_info[11].fileTimeHi = 0U;
  c1_info[11].mFileTimeLo = 0U;
  c1_info[11].mFileTimeHi = 0U;
  c1_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c1_info[12].name = "eml_is_float_class";
  c1_info[12].dominantType = "char";
  c1_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c1_info[12].fileTimeLo = 1286818782U;
  c1_info[12].fileTimeHi = 0U;
  c1_info[12].mFileTimeLo = 0U;
  c1_info[12].mFileTimeHi = 0U;
  c1_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c1_info[13].name = "min";
  c1_info[13].dominantType = "double";
  c1_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c1_info[13].fileTimeLo = 1308747330U;
  c1_info[13].fileTimeHi = 0U;
  c1_info[13].mFileTimeLo = 0U;
  c1_info[13].mFileTimeHi = 0U;
  c1_info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c1_info[14].name = "eml_min_or_max";
  c1_info[14].dominantType = "char";
  c1_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  c1_info[14].fileTimeLo = 1303146212U;
  c1_info[14].fileTimeHi = 0U;
  c1_info[14].mFileTimeLo = 0U;
  c1_info[14].mFileTimeHi = 0U;
  c1_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[15].name = "eml_scalar_eg";
  c1_info[15].dominantType = "double";
  c1_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[15].fileTimeLo = 1286818796U;
  c1_info[15].fileTimeHi = 0U;
  c1_info[15].mFileTimeLo = 0U;
  c1_info[15].mFileTimeHi = 0U;
  c1_info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[16].name = "eml_scalexp_alloc";
  c1_info[16].dominantType = "double";
  c1_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c1_info[16].fileTimeLo = 1286818796U;
  c1_info[16].fileTimeHi = 0U;
  c1_info[16].mFileTimeLo = 0U;
  c1_info[16].mFileTimeHi = 0U;
  c1_info[17].context = "";
  c1_info[17].name = "logspace";
  c1_info[17].dominantType = "double";
  c1_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[17].fileTimeLo = 1305318000U;
  c1_info[17].fileTimeHi = 0U;
  c1_info[17].mFileTimeLo = 0U;
  c1_info[17].mFileTimeHi = 0U;
  c1_info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[18].name = "eml_scalar_floor";
  c1_info[18].dominantType = "double";
  c1_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c1_info[18].fileTimeLo = 1286818726U;
  c1_info[18].fileTimeHi = 0U;
  c1_info[18].mFileTimeLo = 0U;
  c1_info[18].mFileTimeHi = 0U;
  c1_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[19].name = "abs";
  c1_info[19].dominantType = "double";
  c1_info[19].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[19].fileTimeLo = 1286818694U;
  c1_info[19].fileTimeHi = 0U;
  c1_info[19].mFileTimeLo = 0U;
  c1_info[19].mFileTimeHi = 0U;
  c1_info[20].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[20].name = "eml_scalar_abs";
  c1_info[20].dominantType = "double";
  c1_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c1_info[20].fileTimeLo = 1286818712U;
  c1_info[20].fileTimeHi = 0U;
  c1_info[20].mFileTimeLo = 0U;
  c1_info[20].mFileTimeHi = 0U;
  c1_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[21].name = "eps";
  c1_info[21].dominantType = "char";
  c1_info[21].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c1_info[21].fileTimeLo = 1286818686U;
  c1_info[21].fileTimeHi = 0U;
  c1_info[21].mFileTimeLo = 0U;
  c1_info[21].mFileTimeHi = 0U;
  c1_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[22].name = "eml_warning";
  c1_info[22].dominantType = "char";
  c1_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c1_info[22].fileTimeLo = 1286818802U;
  c1_info[22].fileTimeHi = 0U;
  c1_info[22].mFileTimeLo = 0U;
  c1_info[22].mFileTimeHi = 0U;
  c1_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c1_info[23].name = "linspace";
  c1_info[23].dominantType = "double";
  c1_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c1_info[23].fileTimeLo = 1286818762U;
  c1_info[23].fileTimeHi = 0U;
  c1_info[23].mFileTimeLo = 0U;
  c1_info[23].mFileTimeHi = 0U;
  c1_info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c1_info[24].name = "realmax";
  c1_info[24].dominantType = "char";
  c1_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c1_info[24].fileTimeLo = 1286818766U;
  c1_info[24].fileTimeHi = 0U;
  c1_info[24].mFileTimeLo = 0U;
  c1_info[24].mFileTimeHi = 0U;
  c1_info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c1_info[25].name = "mpower";
  c1_info[25].dominantType = "double";
  c1_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c1_info[25].fileTimeLo = 1286818842U;
  c1_info[25].fileTimeHi = 0U;
  c1_info[25].mFileTimeLo = 0U;
  c1_info[25].mFileTimeHi = 0U;
  c1_info[26].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c1_info[26].name = "power";
  c1_info[26].dominantType = "double";
  c1_info[26].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  c1_info[26].fileTimeLo = 1294067944U;
  c1_info[26].fileTimeHi = 0U;
  c1_info[26].mFileTimeLo = 0U;
  c1_info[26].mFileTimeHi = 0U;
  c1_info[27].context = "";
  c1_info[27].name = "diag";
  c1_info[27].dominantType = "double";
  c1_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c1_info[27].fileTimeLo = 1286818686U;
  c1_info[27].fileTimeHi = 0U;
  c1_info[27].mFileTimeLo = 0U;
  c1_info[27].mFileTimeHi = 0U;
  c1_info[28].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c1_info[28].name = "eml_index_plus";
  c1_info[28].dominantType = "int32";
  c1_info[28].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[28].fileTimeLo = 1286818778U;
  c1_info[28].fileTimeHi = 0U;
  c1_info[28].mFileTimeLo = 0U;
  c1_info[28].mFileTimeHi = 0U;
  c1_info[29].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c1_info[29].name = "eml_index_times";
  c1_info[29].dominantType = "int32";
  c1_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c1_info[29].fileTimeLo = 1286818780U;
  c1_info[29].fileTimeHi = 0U;
  c1_info[29].mFileTimeLo = 0U;
  c1_info[29].mFileTimeHi = 0U;
  c1_info[30].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c1_info[30].name = "eml_index_minus";
  c1_info[30].dominantType = "int32";
  c1_info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[30].fileTimeLo = 1286818778U;
  c1_info[30].fileTimeHi = 0U;
  c1_info[30].mFileTimeLo = 0U;
  c1_info[30].mFileTimeHi = 0U;
  c1_info[31].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[31].name = "eml_xgemm";
  c1_info[31].dominantType = "int32";
  c1_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c1_info[31].fileTimeLo = 1299076772U;
  c1_info[31].fileTimeHi = 0U;
  c1_info[31].mFileTimeLo = 0U;
  c1_info[31].mFileTimeHi = 0U;
  c1_info[32].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c1_info[32].name = "eml_blas_inline";
  c1_info[32].dominantType = "";
  c1_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[32].fileTimeLo = 1299076768U;
  c1_info[32].fileTimeHi = 0U;
  c1_info[32].mFileTimeLo = 0U;
  c1_info[32].mFileTimeHi = 0U;
  c1_info[33].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c1_info[33].name = "eml_refblas_xgemm";
  c1_info[33].dominantType = "int32";
  c1_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c1_info[33].fileTimeLo = 1299076774U;
  c1_info[33].fileTimeHi = 0U;
  c1_info[33].mFileTimeLo = 0U;
  c1_info[33].mFileTimeHi = 0U;
  c1_info[34].context = "";
  c1_info[34].name = "length";
  c1_info[34].dominantType = "double";
  c1_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c1_info[34].fileTimeLo = 1303146206U;
  c1_info[34].fileTimeHi = 0U;
  c1_info[34].mFileTimeLo = 0U;
  c1_info[34].mFileTimeHi = 0U;
  c1_info[35].context = "";
  c1_info[35].name = "norm";
  c1_info[35].dominantType = "double";
  c1_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m";
  c1_info[35].fileTimeLo = 1286818826U;
  c1_info[35].fileTimeHi = 0U;
  c1_info[35].mFileTimeLo = 0U;
  c1_info[35].mFileTimeHi = 0U;
  c1_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm";
  c1_info[36].name = "eml_xnrm2";
  c1_info[36].dominantType = "int32";
  c1_info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m";
  c1_info[36].fileTimeLo = 1299076776U;
  c1_info[36].fileTimeHi = 0U;
  c1_info[36].mFileTimeLo = 0U;
  c1_info[36].mFileTimeHi = 0U;
  c1_info[37].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m";
  c1_info[37].name = "eml_refblas_xnrm2";
  c1_info[37].dominantType = "int32";
  c1_info[37].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c1_info[37].fileTimeLo = 1299076784U;
  c1_info[37].fileTimeHi = 0U;
  c1_info[37].mFileTimeLo = 0U;
  c1_info[37].mFileTimeHi = 0U;
  c1_info[38].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c1_info[38].name = "realmin";
  c1_info[38].dominantType = "char";
  c1_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c1_info[38].fileTimeLo = 1286818768U;
  c1_info[38].fileTimeHi = 0U;
  c1_info[38].mFileTimeLo = 0U;
  c1_info[38].mFileTimeHi = 0U;
  c1_info[39].context = "";
  c1_info[39].name = "colon";
  c1_info[39].dominantType = "double";
  c1_info[39].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c1_info[39].fileTimeLo = 1286818838U;
  c1_info[39].fileTimeHi = 0U;
  c1_info[39].mFileTimeLo = 0U;
  c1_info[39].mFileTimeHi = 0U;
  c1_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c1_info[40].name = "isfinite";
  c1_info[40].dominantType = "double";
  c1_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c1_info[40].fileTimeLo = 1286818758U;
  c1_info[40].fileTimeHi = 0U;
  c1_info[40].mFileTimeLo = 0U;
  c1_info[40].mFileTimeHi = 0U;
  c1_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c1_info[41].name = "isnan";
  c1_info[41].dominantType = "double";
  c1_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  c1_info[41].fileTimeLo = 1286818760U;
  c1_info[41].fileTimeHi = 0U;
  c1_info[41].mFileTimeLo = 0U;
  c1_info[41].mFileTimeHi = 0U;
  c1_info[42].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c1_info[42].name = "floor";
  c1_info[42].dominantType = "double";
  c1_info[42].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c1_info[42].fileTimeLo = 1286818742U;
  c1_info[42].fileTimeHi = 0U;
  c1_info[42].mFileTimeLo = 0U;
  c1_info[42].mFileTimeHi = 0U;
  c1_info[43].context = "";
  c1_info[43].name = "eig";
  c1_info[43].dominantType = "double";
  c1_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c1_info[43].fileTimeLo = 1305318000U;
  c1_info[43].fileTimeHi = 0U;
  c1_info[43].mFileTimeLo = 0U;
  c1_info[43].mFileTimeHi = 0U;
  c1_info[44].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c1_info[44].name = "eml_xgeev";
  c1_info[44].dominantType = "double";
  c1_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c1_info[44].fileTimeLo = 1286818804U;
  c1_info[44].fileTimeHi = 0U;
  c1_info[44].mFileTimeLo = 0U;
  c1_info[44].mFileTimeHi = 0U;
  c1_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c1_info[45].name = "eml_lapack_xgeev";
  c1_info[45].dominantType = "double";
  c1_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c1_info[45].fileTimeLo = 1301328468U;
  c1_info[45].fileTimeHi = 0U;
  c1_info[45].mFileTimeLo = 0U;
  c1_info[45].mFileTimeHi = 0U;
  c1_info[46].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c1_info[46].name = "eml_matlab_zggev";
  c1_info[46].dominantType = "double";
  c1_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[46].fileTimeLo = 1286818818U;
  c1_info[46].fileTimeHi = 0U;
  c1_info[46].mFileTimeLo = 0U;
  c1_info[46].mFileTimeHi = 0U;
  c1_info[47].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[47].name = "eml_matlab_zlangeM";
  c1_info[47].dominantType = "double";
  c1_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c1_info[47].fileTimeLo = 1286818820U;
  c1_info[47].fileTimeHi = 0U;
  c1_info[47].mFileTimeLo = 0U;
  c1_info[47].mFileTimeHi = 0U;
  c1_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c1_info[48].name = "eml_dlapy2";
  c1_info[48].dominantType = "double";
  c1_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m";
  c1_info[48].fileTimeLo = 1286818698U;
  c1_info[48].fileTimeHi = 0U;
  c1_info[48].mFileTimeLo = 0U;
  c1_info[48].mFileTimeHi = 0U;
  c1_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c1_info[49].name = "eml_guarded_nan";
  c1_info[49].dominantType = "char";
  c1_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m";
  c1_info[49].fileTimeLo = 1286818776U;
  c1_info[49].fileTimeHi = 0U;
  c1_info[49].mFileTimeLo = 0U;
  c1_info[49].mFileTimeHi = 0U;
  c1_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[50].name = "eml_matlab_zlascl";
  c1_info[50].dominantType = "double";
  c1_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m";
  c1_info[50].fileTimeLo = 1286818822U;
  c1_info[50].fileTimeHi = 0U;
  c1_info[50].mFileTimeLo = 0U;
  c1_info[50].mFileTimeHi = 0U;
  c1_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[51].name = "eml_matlab_zggbal";
  c1_info[51].dominantType = "double";
  c1_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m";
  c1_info[51].fileTimeLo = 1286818818U;
  c1_info[51].fileTimeHi = 0U;
  c1_info[51].mFileTimeLo = 0U;
  c1_info[51].mFileTimeHi = 0U;
  c1_info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[52].name = "eml_matlab_zgghrd";
  c1_info[52].dominantType = "int32";
  c1_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c1_info[52].fileTimeLo = 1286818820U;
  c1_info[52].fileTimeHi = 0U;
  c1_info[52].mFileTimeLo = 0U;
  c1_info[52].mFileTimeHi = 0U;
  c1_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c1_info[53].name = "eml_matlab_zlartg";
  c1_info[53].dominantType = "double";
  c1_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c1_info[53].fileTimeLo = 1286818822U;
  c1_info[53].fileTimeHi = 0U;
  c1_info[53].mFileTimeLo = 0U;
  c1_info[53].mFileTimeHi = 0U;
  c1_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c1_info[54].name = "fix";
  c1_info[54].dominantType = "double";
  c1_info[54].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c1_info[54].fileTimeLo = 1286818742U;
  c1_info[54].fileTimeHi = 0U;
  c1_info[54].mFileTimeLo = 0U;
  c1_info[54].mFileTimeHi = 0U;
  c1_info[55].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c1_info[55].name = "eml_scalar_fix";
  c1_info[55].dominantType = "double";
  c1_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m";
  c1_info[55].fileTimeLo = 1286818726U;
  c1_info[55].fileTimeHi = 0U;
  c1_info[55].mFileTimeLo = 0U;
  c1_info[55].mFileTimeHi = 0U;
  c1_info[56].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c1_info[56].name = "eml_zrot_rows";
  c1_info[56].dominantType = "int32";
  c1_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c1_info[56].fileTimeLo = 1286818826U;
  c1_info[56].fileTimeHi = 0U;
  c1_info[56].mFileTimeLo = 0U;
  c1_info[56].mFileTimeHi = 0U;
  c1_info[57].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c1_info[57].name = "eml_conjtimes";
  c1_info[57].dominantType = "double";
  c1_info[57].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_conjtimes.m";
  c1_info[57].fileTimeLo = 1286818696U;
  c1_info[57].fileTimeHi = 0U;
  c1_info[57].mFileTimeLo = 0U;
  c1_info[57].mFileTimeHi = 0U;
  c1_info[58].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c1_info[58].name = "eml_zrot_cols";
  c1_info[58].dominantType = "int32";
  c1_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m";
  c1_info[58].fileTimeLo = 1286818826U;
  c1_info[58].fileTimeHi = 0U;
  c1_info[58].mFileTimeLo = 0U;
  c1_info[58].mFileTimeHi = 0U;
  c1_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c1_info[59].name = "eml_matlab_zhgeqz";
  c1_info[59].dominantType = "int32";
  c1_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c1_info[59].fileTimeLo = 1286818820U;
  c1_info[59].fileTimeHi = 0U;
  c1_info[59].mFileTimeLo = 0U;
  c1_info[59].mFileTimeHi = 0U;
  c1_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c1_info[60].name = "eml_matlab_zlanhs";
  c1_info[60].dominantType = "int32";
  c1_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m";
  c1_info[60].fileTimeLo = 1286818820U;
  c1_info[60].fileTimeHi = 0U;
  c1_info[60].mFileTimeLo = 0U;
  c1_info[60].mFileTimeHi = 0U;
  c1_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c1_info[61].name = "mod";
  c1_info[61].dominantType = "int32";
  c1_info[61].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c1_info[61].fileTimeLo = 1286818744U;
  c1_info[61].fileTimeHi = 0U;
  c1_info[61].mFileTimeLo = 0U;
  c1_info[61].mFileTimeHi = 0U;
  c1_info[62].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c1_info[62].name = "eml_scalar_mod";
  c1_info[62].dominantType = "int32";
  c1_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_mod.m";
  c1_info[62].fileTimeLo = 1286818730U;
  c1_info[62].fileTimeHi = 0U;
  c1_info[62].mFileTimeLo = 0U;
  c1_info[62].mFileTimeHi = 0U;
  c1_info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c1_info[63].name = "eml_guarded_inf";
  c1_info[63].dominantType = "char";
  c1_info[63].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m";
  c1_info[63].fileTimeLo = 1286818776U;
  c1_info[63].fileTimeHi = 0U;
  c1_info[63].mFileTimeLo = 0U;
  c1_info[63].mFileTimeHi = 0U;
}

static void c1_b_info_helper(c1_ResolvedFunctionInfo c1_info[102])
{
  c1_info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c1_info[64].name = "eml_scalar_hypot";
  c1_info[64].dominantType = "double";
  c1_info[64].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m";
  c1_info[64].fileTimeLo = 1286818726U;
  c1_info[64].fileTimeHi = 0U;
  c1_info[64].mFileTimeLo = 0U;
  c1_info[64].mFileTimeHi = 0U;
  c1_info[65].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c1_info[65].name = "isequal";
  c1_info[65].dominantType = "double";
  c1_info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c1_info[65].fileTimeLo = 1286818758U;
  c1_info[65].fileTimeHi = 0U;
  c1_info[65].mFileTimeLo = 0U;
  c1_info[65].mFileTimeHi = 0U;
  c1_info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c1_info[66].name = "eml_isequal_core";
  c1_info[66].dominantType = "double";
  c1_info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c1_info[66].fileTimeLo = 1286818786U;
  c1_info[66].fileTimeHi = 0U;
  c1_info[66].mFileTimeLo = 0U;
  c1_info[66].mFileTimeHi = 0U;
  c1_info[67].context = "";
  c1_info[67].name = "sort";
  c1_info[67].dominantType = "double";
  c1_info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c1_info[67].fileTimeLo = 1303146208U;
  c1_info[67].fileTimeHi = 0U;
  c1_info[67].mFileTimeLo = 0U;
  c1_info[67].mFileTimeHi = 0U;
  c1_info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c1_info[68].name = "eml_sort";
  c1_info[68].dominantType = "double";
  c1_info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c1_info[68].fileTimeLo = 1305318002U;
  c1_info[68].fileTimeHi = 0U;
  c1_info[68].mFileTimeLo = 0U;
  c1_info[68].mFileTimeHi = 0U;
  c1_info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c1_info[69].name = "eml_nonsingleton_dim";
  c1_info[69].dominantType = "double";
  c1_info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_nonsingleton_dim.m";
  c1_info[69].fileTimeLo = 1286818788U;
  c1_info[69].fileTimeHi = 0U;
  c1_info[69].mFileTimeLo = 0U;
  c1_info[69].mFileTimeHi = 0U;
  c1_info[70].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c1_info[70].name = "eml_assert_valid_dim";
  c1_info[70].dominantType = "double";
  c1_info[70].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_dim.m";
  c1_info[70].fileTimeLo = 1286818694U;
  c1_info[70].fileTimeHi = 0U;
  c1_info[70].mFileTimeLo = 0U;
  c1_info[70].mFileTimeHi = 0U;
  c1_info[71].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c1_info[71].name = "eml_sort_idx";
  c1_info[71].dominantType = "double";
  c1_info[71].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c1_info[71].fileTimeLo = 1305318004U;
  c1_info[71].fileTimeHi = 0U;
  c1_info[71].mFileTimeLo = 0U;
  c1_info[71].mFileTimeHi = 0U;
  c1_info[72].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c1_info[72].name = "eml_size_ispow2";
  c1_info[72].dominantType = "int32";
  c1_info[72].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c1_info[72].fileTimeLo = 1286818798U;
  c1_info[72].fileTimeHi = 0U;
  c1_info[72].mFileTimeLo = 0U;
  c1_info[72].mFileTimeHi = 0U;
  c1_info[73].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c1_info[73].name = "eml_unsigned_class";
  c1_info[73].dominantType = "char";
  c1_info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c1_info[73].fileTimeLo = 1286818800U;
  c1_info[73].fileTimeHi = 0U;
  c1_info[73].mFileTimeLo = 0U;
  c1_info[73].mFileTimeHi = 0U;
  c1_info[74].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c1_info[74].name = "eml_sort_le";
  c1_info[74].dominantType = "int32";
  c1_info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m";
  c1_info[74].fileTimeLo = 1292190510U;
  c1_info[74].fileTimeHi = 0U;
  c1_info[74].mFileTimeLo = 0U;
  c1_info[74].mFileTimeHi = 0U;
  c1_info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m!eml_sort_ascending_le";
  c1_info[75].name = "eml_relop";
  c1_info[75].dominantType = "function_handle";
  c1_info[75].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m";
  c1_info[75].fileTimeLo = 1292190510U;
  c1_info[75].fileTimeHi = 0U;
  c1_info[75].mFileTimeLo = 0U;
  c1_info[75].mFileTimeHi = 0U;
  c1_info[76].context = "";
  c1_info[76].name = "any";
  c1_info[76].dominantType = "double";
  c1_info[76].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c1_info[76].fileTimeLo = 1286818834U;
  c1_info[76].fileTimeHi = 0U;
  c1_info[76].mFileTimeLo = 0U;
  c1_info[76].mFileTimeHi = 0U;
  c1_info[77].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c1_info[77].name = "eml_all_or_any";
  c1_info[77].dominantType = "char";
  c1_info[77].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c1_info[77].fileTimeLo = 1286818694U;
  c1_info[77].fileTimeHi = 0U;
  c1_info[77].mFileTimeLo = 0U;
  c1_info[77].mFileTimeHi = 0U;
  c1_info[78].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c1_info[78].name = "eml_const_nonsingleton_dim";
  c1_info[78].dominantType = "double";
  c1_info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  c1_info[78].fileTimeLo = 1286818696U;
  c1_info[78].fileTimeHi = 0U;
  c1_info[78].mFileTimeLo = 0U;
  c1_info[78].mFileTimeHi = 0U;
  c1_info[79].context = "";
  c1_info[79].name = "log";
  c1_info[79].dominantType = "double";
  c1_info[79].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c1_info[79].fileTimeLo = 1286818742U;
  c1_info[79].fileTimeHi = 0U;
  c1_info[79].mFileTimeLo = 0U;
  c1_info[79].mFileTimeHi = 0U;
  c1_info[80].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c1_info[80].name = "eml_scalar_log";
  c1_info[80].dominantType = "double";
  c1_info[80].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c1_info[80].fileTimeLo = 1286818728U;
  c1_info[80].fileTimeHi = 0U;
  c1_info[80].mFileTimeLo = 0U;
  c1_info[80].mFileTimeHi = 0U;
  c1_info[81].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c1_info[81].name = "eml_scalar_atan2";
  c1_info[81].dominantType = "double";
  c1_info[81].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m";
  c1_info[81].fileTimeLo = 1286818720U;
  c1_info[81].fileTimeHi = 0U;
  c1_info[81].mFileTimeLo = 0U;
  c1_info[81].mFileTimeHi = 0U;
  c1_info[82].context = "";
  c1_info[82].name = "exp";
  c1_info[82].dominantType = "double";
  c1_info[82].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c1_info[82].fileTimeLo = 1286818740U;
  c1_info[82].fileTimeHi = 0U;
  c1_info[82].mFileTimeLo = 0U;
  c1_info[82].mFileTimeHi = 0U;
  c1_info[83].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c1_info[83].name = "eml_scalar_exp";
  c1_info[83].dominantType = "double";
  c1_info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m";
  c1_info[83].fileTimeLo = 1301328464U;
  c1_info[83].fileTimeHi = 0U;
  c1_info[83].mFileTimeLo = 0U;
  c1_info[83].mFileTimeHi = 0U;
  c1_info[84].context = "";
  c1_info[84].name = "mldivide";
  c1_info[84].dominantType = "double";
  c1_info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c1_info[84].fileTimeLo = 1310137456U;
  c1_info[84].fileTimeHi = 0U;
  c1_info[84].mFileTimeLo = 1289519690U;
  c1_info[84].mFileTimeHi = 0U;
  c1_info[85].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c1_info[85].name = "eml_lusolve";
  c1_info[85].dominantType = "double";
  c1_info[85].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c1_info[85].fileTimeLo = 1305318000U;
  c1_info[85].fileTimeHi = 0U;
  c1_info[85].mFileTimeLo = 0U;
  c1_info[85].mFileTimeHi = 0U;
  c1_info[86].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[86].name = "eml_xgetrf";
  c1_info[86].dominantType = "int32";
  c1_info[86].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c1_info[86].fileTimeLo = 1286818806U;
  c1_info[86].fileTimeHi = 0U;
  c1_info[86].mFileTimeLo = 0U;
  c1_info[86].mFileTimeHi = 0U;
  c1_info[87].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c1_info[87].name = "eml_lapack_xgetrf";
  c1_info[87].dominantType = "int32";
  c1_info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c1_info[87].fileTimeLo = 1286818810U;
  c1_info[87].fileTimeHi = 0U;
  c1_info[87].mFileTimeLo = 0U;
  c1_info[87].mFileTimeHi = 0U;
  c1_info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c1_info[88].name = "eml_matlab_zgetrf";
  c1_info[88].dominantType = "int32";
  c1_info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[88].fileTimeLo = 1302688994U;
  c1_info[88].fileTimeHi = 0U;
  c1_info[88].mFileTimeLo = 0U;
  c1_info[88].mFileTimeHi = 0U;
  c1_info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c1_info[89].name = "intmin";
  c1_info[89].dominantType = "char";
  c1_info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c1_info[89].fileTimeLo = 1286818756U;
  c1_info[89].fileTimeHi = 0U;
  c1_info[89].mFileTimeLo = 0U;
  c1_info[89].mFileTimeHi = 0U;
  c1_info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[90].name = "eml_ixamax";
  c1_info[90].dominantType = "int32";
  c1_info[90].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c1_info[90].fileTimeLo = 1299076770U;
  c1_info[90].fileTimeHi = 0U;
  c1_info[90].mFileTimeLo = 0U;
  c1_info[90].mFileTimeHi = 0U;
  c1_info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c1_info[91].name = "eml_refblas_ixamax";
  c1_info[91].dominantType = "int32";
  c1_info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[91].fileTimeLo = 1299076770U;
  c1_info[91].fileTimeHi = 0U;
  c1_info[91].mFileTimeLo = 0U;
  c1_info[91].mFileTimeHi = 0U;
  c1_info[92].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[92].name = "eml_xcabs1";
  c1_info[92].dominantType = "double";
  c1_info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c1_info[92].fileTimeLo = 1286818706U;
  c1_info[92].fileTimeHi = 0U;
  c1_info[92].mFileTimeLo = 0U;
  c1_info[92].mFileTimeHi = 0U;
  c1_info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[93].name = "eml_xswap";
  c1_info[93].dominantType = "int32";
  c1_info[93].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c1_info[93].fileTimeLo = 1299076778U;
  c1_info[93].fileTimeHi = 0U;
  c1_info[93].mFileTimeLo = 0U;
  c1_info[93].mFileTimeHi = 0U;
  c1_info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c1_info[94].name = "eml_refblas_xswap";
  c1_info[94].dominantType = "int32";
  c1_info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[94].fileTimeLo = 1299076786U;
  c1_info[94].fileTimeHi = 0U;
  c1_info[94].mFileTimeLo = 0U;
  c1_info[94].mFileTimeHi = 0U;
  c1_info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[95].name = "eml_xgeru";
  c1_info[95].dominantType = "int32";
  c1_info[95].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c1_info[95].fileTimeLo = 1299076774U;
  c1_info[95].fileTimeHi = 0U;
  c1_info[95].mFileTimeLo = 0U;
  c1_info[95].mFileTimeHi = 0U;
  c1_info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgeru.m";
  c1_info[96].name = "eml_refblas_xgeru";
  c1_info[96].dominantType = "int32";
  c1_info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c1_info[96].fileTimeLo = 1299076776U;
  c1_info[96].fileTimeHi = 0U;
  c1_info[96].mFileTimeLo = 0U;
  c1_info[96].mFileTimeHi = 0U;
  c1_info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c1_info[97].name = "eml_refblas_xgerx";
  c1_info[97].dominantType = "int32";
  c1_info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[97].fileTimeLo = 1299076778U;
  c1_info[97].fileTimeHi = 0U;
  c1_info[97].mFileTimeLo = 0U;
  c1_info[97].mFileTimeHi = 0U;
  c1_info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[98].name = "eml_xtrsm";
  c1_info[98].dominantType = "int32";
  c1_info[98].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c1_info[98].fileTimeLo = 1299076778U;
  c1_info[98].fileTimeHi = 0U;
  c1_info[98].mFileTimeLo = 0U;
  c1_info[98].mFileTimeHi = 0U;
  c1_info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c1_info[99].name = "eml_refblas_xtrsm";
  c1_info[99].dominantType = "int32";
  c1_info[99].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[99].fileTimeLo = 1299076786U;
  c1_info[99].fileTimeHi = 0U;
  c1_info[99].mFileTimeLo = 0U;
  c1_info[99].mFileTimeHi = 0U;
  c1_info[100].context = "";
  c1_info[100].name = "log10";
  c1_info[100].dominantType = "double";
  c1_info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c1_info[100].fileTimeLo = 1286818744U;
  c1_info[100].fileTimeHi = 0U;
  c1_info[100].mFileTimeLo = 0U;
  c1_info[100].mFileTimeHi = 0U;
  c1_info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c1_info[101].name = "eml_scalar_log10";
  c1_info[101].dominantType = "double";
  c1_info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log10.m";
  c1_info[101].fileTimeLo = 1286818728U;
  c1_info[101].fileTimeHi = 0U;
  c1_info[101].mFileTimeLo = 0U;
  c1_info[101].mFileTimeHi = 0U;
}

static void c1_eml_error(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i5;
  static char_T c1_varargin_1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 's', 'q', 'r', 't', '_', 'd', 'o', 'm', 'a',
    'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c1_u[30];
  const mxArray *c1_y = NULL;
  for (c1_i5 = 0; c1_i5 < 30; c1_i5++) {
    c1_u[c1_i5] = c1_varargin_1[c1_i5];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 30), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static real_T c1_mrdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_A, real_T c1_B)
{
  return c1_A / c1_B;
}

static real_T c1_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y)
{
  return c1_x / c1_y;
}

static void c1_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   real_T c1_I[100])
{
  int32_T c1_i;
  int32_T c1_b_i;
  for (c1_i = 0; c1_i < 100; c1_i++) {
    c1_I[c1_i] = 0.0;
  }

  c1_i = 0;
  for (c1_b_i = 0; c1_b_i < 10; c1_b_i++) {
    c1_I[c1_i] = 1.0;
    c1_i += 11;
  }
}

static void c1_b_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[64])
{
  int32_T c1_i;
  int32_T c1_b_i;
  for (c1_i = 0; c1_i < 64; c1_i++) {
    c1_I[c1_i] = 0.0;
  }

  c1_i = 0;
  for (c1_b_i = 0; c1_b_i < 8; c1_b_i++) {
    c1_I[c1_i] = 1.0;
    c1_i += 9;
  }
}

static void c1_logspace(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_d1, real_T c1_d2, real_T c1_y[300])
{
  real_T c1_b_y[300];
  real_T c1_delta1;
  real_T c1_delta2;
  int32_T c1_b_k;
  if (muDoubleScalarAbs(c1_d2 - 3.1415926535897931) < 8.8817841970012523E-16) {
    c1_eml_warning(chartInstance);
  }

  c1_b_y[299] = c1_d2;
  c1_b_y[0] = c1_d1;
  if (((c1_d1 < 0.0 != c1_d2 < 0.0) && (muDoubleScalarAbs(c1_d1) >
        8.9884656743115785E+307)) || (muDoubleScalarAbs(c1_d2) >
       8.9884656743115785E+307)) {
    c1_delta1 = c1_d1 / 299.0;
    c1_delta2 = c1_d2 / 299.0;
    for (c1_b_k = 0; c1_b_k < 298; c1_b_k++) {
      c1_b_y[c1_b_k + 1] = (c1_d1 + c1_delta2 * (1.0 + (real_T)c1_b_k)) -
        c1_delta1 * (1.0 + (real_T)c1_b_k);
    }
  } else {
    c1_delta1 = (c1_d2 - c1_d1) / 299.0;
    for (c1_b_k = 0; c1_b_k < 298; c1_b_k++) {
      c1_b_y[c1_b_k + 1] = c1_d1 + (1.0 + (real_T)c1_b_k) * c1_delta1;
    }
  }

  for (c1_b_k = 0; c1_b_k < 300; c1_b_k++) {
    c1_y[c1_b_k] = muDoubleScalarPower(10.0, c1_b_y[c1_b_k]);
  }
}

static void c1_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i6;
  static char_T c1_varargin_1[32] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'p', 'i', 'S', 'p', 'e', 'c', 'i', 'a', 'l',
    'S', 'u', 'p', 'p', 'o', 'r', 't', 'e', 'd' };

  char_T c1_u[32];
  const mxArray *c1_y = NULL;
  for (c1_i6 = 0; c1_i6 < 32; c1_i6++) {
    c1_u[c1_i6] = c1_varargin_1[c1_i6];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 32), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static void c1_b_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y[2], real_T c1_z[2])
{
  int32_T c1_i7;
  for (c1_i7 = 0; c1_i7 < 2; c1_i7++) {
    c1_z[c1_i7] = c1_x / c1_y[c1_i7];
  }
}

static void c1_diag(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T c1_v[2], real_T c1_d[4])
{
  int32_T c1_j;
  int32_T c1_b_j;
  for (c1_j = 0; c1_j < 4; c1_j++) {
    c1_d[c1_j] = 0.0;
  }

  c1_j = 0;
  for (c1_b_j = 0; c1_b_j < 2; c1_b_j++) {
    c1_d[c1_j] = c1_v[c1_b_j];
    c1_j += 3;
  }
}

static void c1_b_diag(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, real_T c1_v[3], real_T c1_d[9])
{
  int32_T c1_j;
  int32_T c1_b_j;
  for (c1_j = 0; c1_j < 9; c1_j++) {
    c1_d[c1_j] = 0.0;
  }

  c1_j = 0;
  for (c1_b_j = 0; c1_b_j < 3; c1_b_j++) {
    c1_d[c1_j] = c1_v[c1_b_j];
    c1_j += 4;
  }
}

static void c1_c_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[9])
{
  int32_T c1_i;
  int32_T c1_b_i;
  for (c1_i = 0; c1_i < 9; c1_i++) {
    c1_I[c1_i] = 0.0;
  }

  c1_i = 0;
  for (c1_b_i = 0; c1_b_i < 3; c1_b_i++) {
    c1_I[c1_i] = 1.0;
    c1_i += 4;
  }
}

static void c1_c_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x[157], real_T c1_y, real_T c1_z[157])
{
  int32_T c1_i8;
  for (c1_i8 = 0; c1_i8 < 157; c1_i8++) {
    c1_z[c1_i8] = c1_x[c1_i8] / c1_y;
  }
}

static real_T c1_norm(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, real_T c1_x[3])
{
  real_T c1_y;
  real_T c1_scale;
  int32_T c1_b_k;
  real_T c1_absxk;
  real_T c1_t;
  c1_y = 0.0;
  c1_scale = 2.2250738585072014E-308;
  for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
    c1_absxk = muDoubleScalarAbs(c1_x[c1_b_k]);
    if (c1_absxk > c1_scale) {
      c1_t = c1_scale / c1_absxk;
      c1_y = 1.0 + c1_y * c1_t * c1_t;
      c1_scale = c1_absxk;
    } else {
      c1_t = c1_absxk / c1_scale;
      c1_y += c1_t * c1_t;
    }
  }

  return c1_scale * muDoubleScalarSqrt(c1_y);
}

static void c1_d_rdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_x[936], real_T c1_y, real_T c1_z[936])
{
  int32_T c1_i9;
  for (c1_i9 = 0; c1_i9 < 936; c1_i9++) {
    c1_z[c1_i9] = c1_x[c1_i9] / c1_y;
  }
}

static void c1_b_mrdivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_A[153], real_T c1_B, real_T c1_y[153])
{
  int32_T c1_i10;
  for (c1_i10 = 0; c1_i10 < 153; c1_i10++) {
    c1_y[c1_i10] = c1_A[c1_i10] / c1_B;
  }
}

static void c1_damp(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T c1_a[49], real_T c1_h, real_T c1_wn[7], real_T c1_z[7])
{
  creal_T c1_r[7];
  int32_T c1_i;
  real_T c1_x[7];
  int8_T c1_idx[7];
  int32_T c1_b_k;
  boolean_T c1_y;
  int8_T c1_idx0[7];
  int32_T c1_i2;
  int32_T c1_j;
  int32_T c1_pEnd;
  int32_T c1_p;
  int32_T c1_q;
  int32_T c1_qEnd;
  int32_T c1_kEnd;
  creal_T c1_b_r[7];
  boolean_T c1_exitg1;
  boolean_T c1_b0;
  real_T c1_b_a;
  real_T c1_b;
  c1_eig(chartInstance, c1_a, c1_r);
  for (c1_i = 0; c1_i < 7; c1_i++) {
    c1_x[c1_i] = c1_r[c1_i].re;
    c1_idx[c1_i] = (int8_T)(c1_i + 1);
  }

  for (c1_b_k = 0; c1_b_k < 5; c1_b_k += 2) {
    c1_i = c1_b_k + 1;
    if ((c1_x[c1_b_k] <= c1_x[c1_i]) || muDoubleScalarIsNaN(c1_x[c1_i])) {
      c1_y = TRUE;
    } else {
      c1_y = FALSE;
    }

    if (!c1_y) {
      c1_idx[c1_b_k] = (int8_T)(c1_b_k + 2);
      c1_idx[c1_b_k + 1] = (int8_T)(c1_b_k + 1);
    }
  }

  for (c1_i = 0; c1_i < 7; c1_i++) {
    c1_idx0[c1_i] = 1;
  }

  c1_i = 2;
  while (c1_i < 7) {
    c1_i2 = c1_i << 1;
    c1_j = 1;
    for (c1_pEnd = 1 + c1_i; c1_pEnd < 8; c1_pEnd = c1_qEnd + c1_i) {
      c1_p = c1_j;
      c1_q = c1_pEnd;
      c1_qEnd = c1_j + c1_i2;
      if (c1_qEnd > 8) {
        c1_qEnd = 8;
      }

      c1_b_k = 1;
      c1_kEnd = c1_qEnd - c1_j;
      while (c1_b_k <= c1_kEnd) {
        sf_mex_lw_bounds_check(c1_p, 1, 7);
        sf_mex_lw_bounds_check(c1_q, 1, 7);
        if ((c1_x[c1_idx[c1_p - 1] - 1] <= c1_x[c1_idx[c1_q - 1] - 1]) ||
            muDoubleScalarIsNaN(c1_x[c1_idx[c1_q - 1] - 1])) {
          c1_y = TRUE;
        } else {
          c1_y = FALSE;
        }

        if (c1_y) {
          c1_idx0[sf_mex_lw_bounds_check(c1_b_k, 1, 7) - 1] = c1_idx[c1_p - 1];
          c1_p++;
          if (c1_p == c1_pEnd) {
            while (c1_q < c1_qEnd) {
              c1_b_k++;
              c1_idx0[sf_mex_lw_bounds_check(c1_b_k, 1, 7) - 1] = c1_idx[c1_q -
                1];
              c1_q++;
            }
          }
        } else {
          c1_idx0[sf_mex_lw_bounds_check(c1_b_k, 1, 7) - 1] = c1_idx[c1_q - 1];
          c1_q++;
          if (c1_q == c1_qEnd) {
            while (c1_p < c1_pEnd) {
              c1_b_k++;
              c1_idx0[sf_mex_lw_bounds_check(c1_b_k, 1, 7) - 1] = c1_idx[c1_p -
                1];
              c1_p++;
            }
          }
        }

        c1_b_k++;
      }

      for (c1_b_k = 1; c1_b_k <= c1_kEnd; c1_b_k++) {
        c1_idx[sf_mex_lw_bounds_check((c1_j + c1_b_k) - 1, 1, 7) - 1] =
          c1_idx0[c1_b_k - 1];
      }

      c1_j = c1_qEnd;
    }

    c1_i = c1_i2;
  }

  for (c1_i = 0; c1_i < 7; c1_i++) {
    c1_b_r[c1_i] = c1_r[c1_idx[c1_i] - 1];
  }

  for (c1_i = 0; c1_i < 7; c1_i++) {
    c1_r[c1_i] = c1_b_r[c1_i];
  }

  c1_y = FALSE;
  c1_b_k = 0;
  c1_exitg1 = 0U;
  while ((c1_exitg1 == 0U) && (c1_b_k < 7)) {
    if (((c1_r[c1_b_k].re == 0.0) && (c1_r[c1_b_k].im == 0.0)) ||
        (muDoubleScalarIsNaN(c1_r[c1_b_k].re) || muDoubleScalarIsNaN(c1_r[c1_b_k]
          .im))) {
      c1_b0 = TRUE;
    } else {
      c1_b0 = FALSE;
    }

    if (!c1_b0) {
      c1_y = TRUE;
      c1_exitg1 = 1U;
    } else {
      c1_b_k++;
    }
  }

  if (c1_y == 0) {
    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_z[c1_i] = 0.0;
      c1_wn[c1_i] = 0.0;
    }
  } else {
    for (c1_b_k = 0; c1_b_k < 7; c1_b_k++) {
      c1_b_a = c1_r[c1_b_k].re;
      c1_b = c1_r[c1_b_k].im;
      if ((c1_r[c1_b_k].im == 0.0) && muDoubleScalarIsNaN(c1_r[c1_b_k].re)) {
      } else if ((muDoubleScalarAbs(c1_r[c1_b_k].re) > 8.9884656743115785E+307) ||
                 (muDoubleScalarAbs(c1_r[c1_b_k].im) > 8.9884656743115785E+307))
      {
        c1_b_a = muDoubleScalarAbs(c1_r[c1_b_k].re / 2.0);
        c1_b = muDoubleScalarAbs(c1_r[c1_b_k].im / 2.0);
        if (c1_b_a < c1_b) {
          c1_b_a /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_b_a * c1_b_a + 1.0);
        } else if (c1_b_a > c1_b) {
          c1_b /= c1_b_a;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_b_a * 1.4142135623730951;
          }
        }

        c1_b_a = muDoubleScalarLog(c1_b) + 0.69314718055994529;
        c1_b = muDoubleScalarAtan2(c1_r[c1_b_k].im, c1_r[c1_b_k].re);
      } else {
        c1_b_a = muDoubleScalarAbs(c1_r[c1_b_k].re);
        c1_b = muDoubleScalarAbs(c1_r[c1_b_k].im);
        if (c1_b_a < c1_b) {
          c1_b_a /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_b_a * c1_b_a + 1.0);
        } else if (c1_b_a > c1_b) {
          c1_b /= c1_b_a;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_b_a * 1.4142135623730951;
          }
        }

        c1_b_a = muDoubleScalarLog(c1_b);
        c1_b = muDoubleScalarAtan2(c1_r[c1_b_k].im, c1_r[c1_b_k].re);
      }

      c1_r[c1_b_k].re = c1_b_a;
      c1_r[c1_b_k].im = c1_b;
    }

    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_b_a = c1_r[c1_i].im;
      if (c1_r[c1_i].im == 0.0) {
        c1_r[c1_i].re /= c1_h;
        c1_r[c1_i].im = 0.0;
      } else if (c1_r[c1_i].re == 0.0) {
        c1_r[c1_i].re = 0.0;
        c1_r[c1_i].im = c1_b_a / c1_h;
      } else {
        c1_r[c1_i].re /= c1_h;
        c1_r[c1_i].im = c1_b_a / c1_h;
      }
    }

    for (c1_b_k = 0; c1_b_k < 7; c1_b_k++) {
      c1_b_a = muDoubleScalarAbs(c1_r[c1_b_k].re);
      c1_b = muDoubleScalarAbs(c1_r[c1_b_k].im);
      if (c1_b_a < c1_b) {
        c1_b_a /= c1_b;
        c1_b *= muDoubleScalarSqrt(c1_b_a * c1_b_a + 1.0);
      } else if (c1_b_a > c1_b) {
        c1_b /= c1_b_a;
        c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_b_a;
      } else {
        if (!muDoubleScalarIsNaN(c1_b)) {
          c1_b = c1_b_a * 1.4142135623730951;
        }
      }

      c1_wn[c1_b_k] = c1_b;
    }

    for (c1_i = 0; c1_i < 7; c1_i++) {
      c1_z[c1_i] = -c1_r[c1_i].re / c1_wn[c1_i];
    }
  }
}

static void c1_eig(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   real_T c1_A[49], creal_T c1_V[7])
{
  int32_T c1_ii;
  creal_T c1_b_A[49];
  real_T c1_info;
  real_T c1_anrm;
  boolean_T c1_exitg7;
  real_T c1_a;
  real_T c1_b;
  creal_T c1_beta1[7];
  boolean_T c1_ilascl;
  real_T c1_anrmto;
  boolean_T c1_notdone;
  real_T c1_cfrom1;
  real_T c1_cto1;
  real_T c1_mul;
  int32_T c1_ilo;
  int32_T c1_ihi;
  int32_T c1_exitg2;
  int32_T c1_i;
  int32_T c1_j;
  boolean_T c1_exitg5;
  int32_T c1_nzcount;
  int32_T c1_jrow;
  boolean_T c1_exitg6;
  boolean_T c1_c_A;
  boolean_T c1_guard2 = FALSE;
  creal_T c1_atmp;
  int32_T c1_exitg1;
  boolean_T c1_exitg3;
  boolean_T c1_exitg4;
  boolean_T c1_guard1 = FALSE;
  int32_T c1_jrowm1;
  creal_T c1_s;
  static creal_T c1_dc0 = { 0.0, 0.0 };

  creal_T c1_d_A;
  creal_T c1_e_A;
  for (c1_ii = 0; c1_ii < 49; c1_ii++) {
    c1_b_A[c1_ii].re = c1_A[c1_ii];
    c1_b_A[c1_ii].im = 0.0;
  }

  c1_info = 0.0;
  c1_anrm = 0.0;
  c1_ii = 0;
  c1_exitg7 = 0U;
  while ((c1_exitg7 == 0U) && (c1_ii < 49)) {
    c1_a = muDoubleScalarAbs(c1_b_A[c1_ii].re);
    if (c1_a > 0.0) {
      c1_b = 0.0 / c1_a;
      c1_a *= muDoubleScalarSqrt(c1_b * c1_b + 1.0);
    } else {
      c1_a *= 1.4142135623730951;
    }

    if (muDoubleScalarIsNaN(c1_a)) {
      c1_anrm = rtNaN;
      c1_exitg7 = 1U;
    } else {
      if (c1_a > c1_anrm) {
        c1_anrm = c1_a;
      }

      c1_ii++;
    }
  }

  if (!((!muDoubleScalarIsInf(c1_anrm)) && (!muDoubleScalarIsNaN(c1_anrm)))) {
    for (c1_ii = 0; c1_ii < 7; c1_ii++) {
      c1_V[c1_ii].re = rtNaN;
      c1_V[c1_ii].im = 0.0;
      c1_beta1[c1_ii].re = rtNaN;
      c1_beta1[c1_ii].im = 0.0;
    }
  } else {
    c1_ilascl = FALSE;
    c1_anrmto = c1_anrm;
    if ((c1_anrm > 0.0) && (c1_anrm < 6.7178761075670888E-139)) {
      c1_anrmto = 6.7178761075670888E-139;
      c1_ilascl = TRUE;
    } else {
      if (c1_anrm > 1.4885657073574029E+138) {
        c1_anrmto = 1.4885657073574029E+138;
        c1_ilascl = TRUE;
      }
    }

    if (c1_ilascl) {
      c1_a = c1_anrm;
      c1_b = c1_anrmto;
      c1_notdone = TRUE;
      while (c1_notdone) {
        c1_cfrom1 = c1_a * 2.0041683600089728E-292;
        c1_cto1 = c1_b / 4.9896007738368E+291;
        if ((c1_cfrom1 > c1_b) && (c1_b != 0.0)) {
          c1_mul = 2.0041683600089728E-292;
          c1_a = c1_cfrom1;
        } else if (c1_cto1 > c1_a) {
          c1_mul = 4.9896007738368E+291;
          c1_b = c1_cto1;
        } else {
          c1_mul = c1_b / c1_a;
          c1_notdone = FALSE;
        }

        for (c1_ii = 0; c1_ii < 49; c1_ii++) {
          c1_b_A[c1_ii].re *= c1_mul;
          c1_b_A[c1_ii].im *= c1_mul;
        }
      }
    }

    c1_ilo = 1;
    c1_ihi = 7;
    do {
      c1_exitg2 = 0U;
      c1_i = 0;
      c1_j = 0;
      c1_notdone = FALSE;
      c1_ii = c1_ihi;
      c1_exitg5 = 0U;
      while ((c1_exitg5 == 0U) && (c1_ii > 0)) {
        c1_nzcount = 0;
        c1_i = c1_ii;
        c1_j = c1_ihi;
        c1_jrow = 1;
        c1_exitg6 = 0U;
        while ((c1_exitg6 == 0U) && (c1_jrow <= c1_ihi)) {
          c1_c_A = ((c1_b_A[(sf_mex_lw_bounds_check(c1_ii, 1, 7) + 7 * (c1_jrow
            - 1)) - 1].re != 0.0) || (c1_b_A[(sf_mex_lw_bounds_check(c1_ii, 1, 7)
                      + 7 * (c1_jrow - 1)) - 1].im != 0.0));
          c1_guard2 = FALSE;
          if (c1_c_A || (c1_ii == c1_jrow)) {
            if (c1_nzcount == 0) {
              c1_j = c1_jrow;
              c1_nzcount = 1;
              c1_guard2 = TRUE;
            } else {
              c1_nzcount = 2;
              c1_exitg6 = 1U;
            }
          } else {
            c1_guard2 = TRUE;
          }

          if (c1_guard2 == TRUE) {
            c1_jrow++;
          }
        }

        if (c1_nzcount < 2) {
          c1_notdone = TRUE;
          c1_exitg5 = 1U;
        } else {
          c1_ii--;
        }
      }

      if (!c1_notdone) {
        c1_exitg2 = 2U;
      } else {
        if (c1_i != c1_ihi) {
          for (c1_ii = 0; c1_ii < 7; c1_ii++) {
            c1_atmp = c1_b_A[(sf_mex_lw_bounds_check(c1_i, 1, 7) + 7 * c1_ii) -
              1];
            c1_b_A[(c1_i + 7 * c1_ii) - 1] = c1_b_A[(sf_mex_lw_bounds_check
              (c1_ihi, 1, 7) + 7 * c1_ii) - 1];
            c1_b_A[(c1_ihi + 7 * c1_ii) - 1] = c1_atmp;
          }
        }

        if (c1_j != c1_ihi) {
          for (c1_ii = 0; c1_ii + 1 <= c1_ihi; c1_ii++) {
            c1_atmp = c1_b_A[c1_ii + 7 * (sf_mex_lw_bounds_check(c1_j, 1, 7) - 1)];
            c1_b_A[c1_ii + 7 * (c1_j - 1)] = c1_b_A[c1_ii + 7 * (c1_ihi - 1)];
            c1_b_A[c1_ii + 7 * (c1_ihi - 1)] = c1_atmp;
          }
        }

        sf_mex_lw_bounds_check(c1_ihi, 1, 7);
        c1_ihi--;
        if (c1_ihi == 1) {
          c1_exitg2 = 1U;
        }
      }
    } while (c1_exitg2 == 0U);

    if (c1_exitg2 == 1U) {
    } else {
      do {
        c1_exitg1 = 0U;
        c1_i = 0;
        c1_j = 0;
        c1_notdone = FALSE;
        c1_jrow = c1_ilo;
        c1_exitg3 = 0U;
        while ((c1_exitg3 == 0U) && (c1_jrow <= c1_ihi)) {
          c1_nzcount = 0;
          c1_i = c1_ihi;
          c1_j = c1_jrow;
          c1_ii = c1_ilo;
          c1_exitg4 = 0U;
          while ((c1_exitg4 == 0U) && (c1_ii <= c1_ihi)) {
            c1_c_A = ((c1_b_A[(c1_ii + 7 * (sf_mex_lw_bounds_check(c1_jrow, 1, 7)
              - 1)) - 1].re != 0.0) || (c1_b_A[(c1_ii + 7 *
                        (sf_mex_lw_bounds_check(c1_jrow, 1, 7) - 1)) - 1].im !=
                       0.0));
            c1_guard1 = FALSE;
            if (c1_c_A || (c1_ii == c1_jrow)) {
              if (c1_nzcount == 0) {
                c1_i = c1_ii;
                c1_nzcount = 1;
                c1_guard1 = TRUE;
              } else {
                c1_nzcount = 2;
                c1_exitg4 = 1U;
              }
            } else {
              c1_guard1 = TRUE;
            }

            if (c1_guard1 == TRUE) {
              c1_ii++;
            }
          }

          if (c1_nzcount < 2) {
            c1_notdone = TRUE;
            c1_exitg3 = 1U;
          } else {
            c1_jrow++;
          }
        }

        if (!c1_notdone) {
          c1_exitg1 = 1U;
        } else {
          if (c1_i != c1_ilo) {
            for (c1_ii = c1_ilo - 1; c1_ii + 1 < 8; c1_ii++) {
              c1_atmp = c1_b_A[(sf_mex_lw_bounds_check(c1_i, 1, 7) + 7 * c1_ii)
                - 1];
              c1_b_A[(c1_i + 7 * c1_ii) - 1] = c1_b_A[(sf_mex_lw_bounds_check
                (c1_ilo, 1, 7) + 7 * c1_ii) - 1];
              c1_b_A[(c1_ilo + 7 * c1_ii) - 1] = c1_atmp;
            }
          }

          if (c1_j != c1_ilo) {
            for (c1_ii = 0; c1_ii + 1 <= c1_ihi; c1_ii++) {
              c1_atmp = c1_b_A[c1_ii + 7 * (sf_mex_lw_bounds_check(c1_j, 1, 7) -
                1)];
              c1_b_A[c1_ii + 7 * (c1_j - 1)] = c1_b_A[c1_ii + 7 *
                (sf_mex_lw_bounds_check(c1_ilo, 1, 7) - 1)];
              c1_b_A[c1_ii + 7 * (c1_ilo - 1)] = c1_atmp;
            }
          }

          sf_mex_lw_bounds_check(c1_ilo, 1, 7);
          c1_ilo++;
          if (c1_ilo == c1_ihi) {
            c1_exitg1 = 1U;
          }
        }
      } while (c1_exitg1 == 0U);
    }

    if (!(c1_ihi < c1_ilo + 2)) {
      c1_ii = c1_ilo - 1;
      while (c1_ii + 1 < c1_ihi - 1) {
        c1_nzcount = c1_ii + 1;
        c1_jrow = c1_ihi - 1;
        while (c1_jrow + 1 > c1_nzcount + 1) {
          c1_jrowm1 = c1_jrow - 1;
          c1_eml_matlab_zlartg(chartInstance, c1_b_A[c1_jrowm1 + 7 * c1_ii],
                               c1_b_A[c1_jrow + 7 * c1_ii], &c1_a, &c1_s,
                               &c1_atmp);
          c1_b_A[c1_jrowm1 + 7 * c1_ii] = c1_atmp;
          c1_b_A[c1_jrow + 7 * c1_ii] = c1_dc0;
          for (c1_j = c1_nzcount; c1_j + 1 <= c1_ihi; c1_j++) {
            c1_atmp.re = c1_a * c1_b_A[c1_jrowm1 + 7 * c1_j].re;
            c1_atmp.im = c1_a * c1_b_A[c1_jrowm1 + 7 * c1_j].im;
            c1_b = c1_s.re * c1_b_A[c1_jrow + 7 * c1_j].re - c1_s.im *
              c1_b_A[c1_jrow + 7 * c1_j].im;
            c1_cfrom1 = c1_s.re * c1_b_A[c1_jrow + 7 * c1_j].im + c1_s.im *
              c1_b_A[c1_jrow + 7 * c1_j].re;
            c1_d_A = c1_b_A[c1_jrowm1 + 7 * c1_j];
            c1_e_A = c1_b_A[c1_jrowm1 + 7 * c1_j];
            c1_b_A[c1_jrow + 7 * c1_j].re = c1_a * c1_b_A[c1_jrow + 7 * c1_j].re
              - (c1_s.re * c1_b_A[c1_jrowm1 + 7 * c1_j].re + c1_s.im *
                 c1_b_A[c1_jrowm1 + 7 * c1_j].im);
            c1_b_A[c1_jrow + 7 * c1_j].im = c1_a * c1_b_A[c1_jrow + 7 * c1_j].im
              - (c1_s.re * c1_d_A.im - c1_s.im * c1_e_A.re);
            c1_b_A[c1_jrowm1 + 7 * c1_j].re = c1_atmp.re + c1_b;
            c1_b_A[c1_jrowm1 + 7 * c1_j].im = c1_atmp.im + c1_cfrom1;
          }

          c1_s.re = -c1_s.re;
          c1_s.im = -c1_s.im;
          for (c1_i = c1_ilo - 1; c1_i + 1 <= c1_ihi; c1_i++) {
            c1_atmp.re = c1_a * c1_b_A[c1_i + 7 * c1_jrow].re;
            c1_atmp.im = c1_a * c1_b_A[c1_i + 7 * c1_jrow].im;
            c1_b = c1_s.re * c1_b_A[c1_i + 7 * c1_jrowm1].re - c1_s.im *
              c1_b_A[c1_i + 7 * c1_jrowm1].im;
            c1_cfrom1 = c1_s.re * c1_b_A[c1_i + 7 * c1_jrowm1].im + c1_s.im *
              c1_b_A[c1_i + 7 * c1_jrowm1].re;
            c1_d_A = c1_b_A[c1_i + 7 * c1_jrow];
            c1_e_A = c1_b_A[c1_i + 7 * c1_jrow];
            c1_b_A[c1_i + 7 * c1_jrowm1].re = c1_a * c1_b_A[c1_i + 7 * c1_jrowm1]
              .re - (c1_s.re * c1_b_A[c1_i + 7 * c1_jrow].re + c1_s.im *
                     c1_b_A[c1_i + 7 * c1_jrow].im);
            c1_b_A[c1_i + 7 * c1_jrowm1].im = c1_a * c1_b_A[c1_i + 7 * c1_jrowm1]
              .im - (c1_s.re * c1_d_A.im - c1_s.im * c1_e_A.re);
            c1_b_A[c1_i + 7 * c1_jrow].re = c1_atmp.re + c1_b;
            c1_b_A[c1_i + 7 * c1_jrow].im = c1_atmp.im + c1_cfrom1;
          }

          c1_jrow = c1_jrowm1;
        }

        c1_ii = c1_nzcount;
      }
    }

    c1_eml_matlab_zhgeqz(chartInstance, c1_b_A, c1_ilo, c1_ihi, &c1_info, c1_V,
                         c1_beta1);
    if ((!(c1_info != 0.0)) && c1_ilascl) {
      c1_notdone = TRUE;
      while (c1_notdone) {
        c1_cfrom1 = c1_anrmto * 2.0041683600089728E-292;
        c1_cto1 = c1_anrm / 4.9896007738368E+291;
        if ((c1_cfrom1 > c1_anrm) && (c1_anrm != 0.0)) {
          c1_mul = 2.0041683600089728E-292;
          c1_anrmto = c1_cfrom1;
        } else if (c1_cto1 > c1_anrmto) {
          c1_mul = 4.9896007738368E+291;
          c1_anrm = c1_cto1;
        } else {
          c1_mul = c1_anrm / c1_anrmto;
          c1_notdone = FALSE;
        }

        for (c1_ii = 0; c1_ii < 7; c1_ii++) {
          c1_V[c1_ii].re *= c1_mul;
          c1_V[c1_ii].im *= c1_mul;
        }
      }
    }
  }

  for (c1_ii = 0; c1_ii < 7; c1_ii++) {
    c1_cto1 = c1_V[c1_ii].re;
    c1_mul = c1_V[c1_ii].im;
    if (c1_beta1[c1_ii].im == 0.0) {
      if (c1_V[c1_ii].im == 0.0) {
        c1_V[c1_ii].re /= c1_beta1[c1_ii].re;
        c1_V[c1_ii].im = 0.0;
      } else if (c1_V[c1_ii].re == 0.0) {
        c1_V[c1_ii].re = 0.0;
        c1_V[c1_ii].im = c1_mul / c1_beta1[c1_ii].re;
      } else {
        c1_V[c1_ii].re /= c1_beta1[c1_ii].re;
        c1_V[c1_ii].im = c1_mul / c1_beta1[c1_ii].re;
      }
    } else if (c1_beta1[c1_ii].re == 0.0) {
      if (c1_V[c1_ii].re == 0.0) {
        c1_V[c1_ii].re = c1_V[c1_ii].im / c1_beta1[c1_ii].im;
        c1_V[c1_ii].im = 0.0;
      } else if (c1_V[c1_ii].im == 0.0) {
        c1_V[c1_ii].re = 0.0;
        c1_V[c1_ii].im = -(c1_cto1 / c1_beta1[c1_ii].im);
      } else {
        c1_V[c1_ii].re = c1_V[c1_ii].im / c1_beta1[c1_ii].im;
        c1_V[c1_ii].im = -(c1_cto1 / c1_beta1[c1_ii].im);
      }
    } else {
      c1_cfrom1 = muDoubleScalarAbs(c1_beta1[c1_ii].re);
      c1_a = muDoubleScalarAbs(c1_beta1[c1_ii].im);
      if (c1_cfrom1 > c1_a) {
        c1_a = c1_beta1[c1_ii].im / c1_beta1[c1_ii].re;
        c1_b = c1_beta1[c1_ii].re + c1_a * c1_beta1[c1_ii].im;
        c1_V[c1_ii].re = (c1_V[c1_ii].re + c1_a * c1_V[c1_ii].im) / c1_b;
        c1_V[c1_ii].im = (c1_mul - c1_a * c1_cto1) / c1_b;
      } else if (c1_a == c1_cfrom1) {
        c1_a = c1_beta1[c1_ii].re > 0.0 ? 0.5 : -0.5;
        c1_b = c1_beta1[c1_ii].im > 0.0 ? 0.5 : -0.5;
        c1_V[c1_ii].re = (c1_V[c1_ii].re * c1_a + c1_V[c1_ii].im * c1_b) /
          c1_cfrom1;
        c1_V[c1_ii].im = (c1_mul * c1_a - c1_cto1 * c1_b) / c1_cfrom1;
      } else {
        c1_a = c1_beta1[c1_ii].re / c1_beta1[c1_ii].im;
        c1_b = c1_beta1[c1_ii].im + c1_a * c1_beta1[c1_ii].re;
        c1_V[c1_ii].re = (c1_a * c1_V[c1_ii].re + c1_V[c1_ii].im) / c1_b;
        c1_V[c1_ii].im = (c1_a * c1_mul - c1_cto1) / c1_b;
      }
    }
  }

  if (c1_info < 0.0) {
    c1_b_eml_warning(chartInstance);
  } else {
    if (c1_info > 0.0) {
      c1_c_eml_warning(chartInstance);
    }
  }
}

static void c1_eml_matlab_zlartg(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn,
  creal_T *c1_r)
{
  real_T c1_scale;
  real_T c1_dr;
  real_T c1_f2;
  real_T c1_fs_re;
  real_T c1_fs_im;
  real_T c1_gs_re;
  real_T c1_gs_im;
  int32_T c1_count;
  int32_T c1_rescaledir;
  boolean_T c1_guard1 = FALSE;
  real_T c1_b;
  real_T c1_g2s;
  c1_scale = muDoubleScalarAbs(c1_f.re);
  c1_dr = muDoubleScalarAbs(c1_f.im);
  if (c1_dr > c1_scale) {
    c1_scale = c1_dr;
  }

  c1_dr = muDoubleScalarAbs(c1_g.re);
  c1_f2 = muDoubleScalarAbs(c1_g.im);
  if (c1_f2 > c1_dr) {
    c1_dr = c1_f2;
  }

  if (c1_dr > c1_scale) {
    c1_scale = c1_dr;
  }

  c1_fs_re = c1_f.re;
  c1_fs_im = c1_f.im;
  c1_gs_re = c1_g.re;
  c1_gs_im = c1_g.im;
  c1_count = 0;
  c1_rescaledir = 0;
  c1_guard1 = FALSE;
  if (c1_scale >= 7.4428285367870146E+137) {
    do {
      c1_count++;
      c1_fs_re *= 1.3435752215134178E-138;
      c1_fs_im *= 1.3435752215134178E-138;
      c1_gs_re *= 1.3435752215134178E-138;
      c1_gs_im *= 1.3435752215134178E-138;
      c1_scale *= 1.3435752215134178E-138;
    } while (!(c1_scale < 7.4428285367870146E+137));

    c1_rescaledir = 1;
    c1_guard1 = TRUE;
  } else if (c1_scale <= 1.3435752215134178E-138) {
    if ((c1_g.re == 0.0) && (c1_g.im == 0.0)) {
      *c1_cs = 1.0;
      c1_sn->re = 0.0;
      c1_sn->im = 0.0;
      *c1_r = c1_f;
    } else {
      do {
        c1_count++;
        c1_fs_re *= 7.4428285367870146E+137;
        c1_fs_im *= 7.4428285367870146E+137;
        c1_gs_re *= 7.4428285367870146E+137;
        c1_gs_im *= 7.4428285367870146E+137;
        c1_scale *= 7.4428285367870146E+137;
      } while (!(c1_scale > 1.3435752215134178E-138));

      c1_rescaledir = -1;
      c1_guard1 = TRUE;
    }
  } else {
    c1_guard1 = TRUE;
  }

  if (c1_guard1 == TRUE) {
    c1_f2 = c1_fs_re * c1_fs_re + c1_fs_im * c1_fs_im;
    c1_scale = c1_gs_re * c1_gs_re + c1_gs_im * c1_gs_im;
    c1_dr = c1_scale;
    if (1.0 > c1_scale) {
      c1_dr = 1.0;
    }

    if (c1_f2 <= c1_dr * 2.0041683600089728E-292) {
      if ((c1_f.re == 0.0) && (c1_f.im == 0.0)) {
        *c1_cs = 0.0;
        c1_f2 = muDoubleScalarAbs(c1_g.re);
        c1_b = muDoubleScalarAbs(c1_g.im);
        if (c1_f2 < c1_b) {
          c1_f2 /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
        } else if (c1_f2 > c1_b) {
          c1_b /= c1_f2;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_f2 * 1.4142135623730951;
          }
        }

        c1_r->re = c1_b;
        c1_r->im = 0.0;
        c1_f2 = muDoubleScalarAbs(c1_gs_re);
        c1_b = muDoubleScalarAbs(c1_gs_im);
        if (c1_f2 < c1_b) {
          c1_f2 /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
        } else if (c1_f2 > c1_b) {
          c1_b /= c1_f2;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_f2 * 1.4142135623730951;
          }
        }

        c1_sn->re = c1_gs_re / c1_b;
        c1_sn->im = -c1_gs_im / c1_b;
      } else {
        c1_f2 = muDoubleScalarAbs(c1_fs_re);
        c1_b = muDoubleScalarAbs(c1_fs_im);
        if (c1_f2 < c1_b) {
          c1_f2 /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
        } else if (c1_f2 > c1_b) {
          c1_b /= c1_f2;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_f2 * 1.4142135623730951;
          }
        }

        if (c1_scale < 0.0) {
          c1_eml_error(chartInstance);
        }

        c1_g2s = muDoubleScalarSqrt(c1_scale);
        *c1_cs = c1_b / c1_g2s;
        c1_dr = muDoubleScalarAbs(c1_f.re);
        c1_f2 = muDoubleScalarAbs(c1_f.im);
        if (c1_f2 > c1_dr) {
          c1_dr = c1_f2;
        }

        if (c1_dr > 1.0) {
          c1_f2 = muDoubleScalarAbs(c1_f.re);
          c1_b = muDoubleScalarAbs(c1_f.im);
          if (c1_f2 < c1_b) {
            c1_f2 /= c1_b;
            c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
          } else if (c1_f2 > c1_b) {
            c1_b /= c1_f2;
            c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
          } else {
            if (!muDoubleScalarIsNaN(c1_b)) {
              c1_b = c1_f2 * 1.4142135623730951;
            }
          }

          c1_fs_re = c1_f.re / c1_b;
          c1_fs_im = c1_f.im / c1_b;
        } else {
          c1_dr = 7.4428285367870146E+137 * c1_f.re;
          c1_scale = 7.4428285367870146E+137 * c1_f.im;
          c1_f2 = muDoubleScalarAbs(c1_dr);
          c1_b = muDoubleScalarAbs(c1_scale);
          if (c1_f2 < c1_b) {
            c1_f2 /= c1_b;
            c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
          } else if (c1_f2 > c1_b) {
            c1_b /= c1_f2;
            c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
          } else {
            if (!muDoubleScalarIsNaN(c1_b)) {
              c1_b = c1_f2 * 1.4142135623730951;
            }
          }

          c1_fs_re = c1_dr / c1_b;
          c1_fs_im = c1_scale / c1_b;
        }

        c1_gs_re /= c1_g2s;
        c1_gs_im = -c1_gs_im / c1_g2s;
        c1_sn->re = c1_fs_re * c1_gs_re - c1_fs_im * c1_gs_im;
        c1_sn->im = c1_fs_re * c1_gs_im + c1_fs_im * c1_gs_re;
        c1_r->re = *c1_cs * c1_f.re + (c1_sn->re * c1_g.re - c1_sn->im * c1_g.im);
        c1_r->im = *c1_cs * c1_f.im + (c1_sn->re * c1_g.im + c1_sn->im * c1_g.re);
      }
    } else {
      c1_dr = 1.0 + c1_scale / c1_f2;
      if (c1_dr < 0.0) {
        c1_eml_error(chartInstance);
      }

      c1_b = muDoubleScalarSqrt(c1_dr);
      c1_r->re = c1_b * c1_fs_re;
      c1_r->im = c1_b * c1_fs_im;
      *c1_cs = 1.0 / c1_b;
      c1_b = c1_f2 + c1_scale;
      c1_f2 = c1_r->re / c1_b;
      c1_dr = c1_r->im / c1_b;
      c1_sn->re = c1_f2 * c1_gs_re - c1_dr * -c1_gs_im;
      c1_sn->im = c1_f2 * -c1_gs_im + c1_dr * c1_gs_re;
      if (c1_rescaledir > 0) {
        for (c1_rescaledir = 1; c1_rescaledir <= c1_count; c1_rescaledir++) {
          c1_r->re *= 7.4428285367870146E+137;
          c1_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c1_rescaledir < 0) {
          for (c1_rescaledir = 1; c1_rescaledir <= c1_count; c1_rescaledir++) {
            c1_r->re *= 1.3435752215134178E-138;
            c1_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c1_eml_matlab_zhgeqz(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_A[49], int32_T c1_ilo, int32_T c1_ihi, real_T
  *c1_info, creal_T c1_alpha1[7], creal_T c1_beta1[7])
{
  static creal_T c1_dc1 = { 0.0, 0.0 };

  int32_T c1_jm1;
  creal_T c1_b_A[49];
  real_T c1_eshift_re;
  real_T c1_eshift_im;
  static creal_T c1_dc2 = { 0.0, 0.0 };

  creal_T c1_ctemp;
  real_T c1_rho_re;
  real_T c1_rho_im;
  real_T c1_anorm;
  real_T c1_a;
  real_T c1_atol;
  real_T c1_ascale;
  boolean_T c1_failed;
  int32_T c1_j;
  boolean_T c1_guard1 = FALSE;
  boolean_T c1_guard2 = FALSE;
  int32_T c1_ifirst;
  int32_T c1_istart;
  int32_T c1_ilast;
  int32_T c1_ilastm1;
  int32_T c1_ifrstm;
  int32_T c1_iiter;
  int32_T c1_maxit;
  boolean_T c1_goto60;
  boolean_T c1_goto70;
  boolean_T c1_goto90;
  int32_T c1_jiter;
  int32_T c1_exitg1;
  boolean_T c1_exitg3;
  boolean_T c1_ilazro;
  boolean_T c1_b_guard1 = FALSE;
  creal_T c1_c_A;
  creal_T c1_t1;
  creal_T c1_d;
  creal_T c1_sigma1;
  real_T c1_sigma2_re;
  real_T c1_sigma2_im;
  real_T c1_b;
  int32_T c1_jp1;
  boolean_T c1_exitg2;
  creal_T c1_d_A;
  int32_T c1_i;
  c1_dc1.re = rtNaN;
  for (c1_jm1 = 0; c1_jm1 < 49; c1_jm1++) {
    c1_b_A[c1_jm1] = c1_A[c1_jm1];
  }

  for (c1_jm1 = 0; c1_jm1 < 7; c1_jm1++) {
    c1_alpha1[c1_jm1].re = 0.0;
    c1_alpha1[c1_jm1].im = 0.0;
    c1_beta1[c1_jm1].re = 1.0;
    c1_beta1[c1_jm1].im = 0.0;
  }

  c1_eshift_re = 0.0;
  c1_eshift_im = 0.0;
  c1_ctemp = c1_dc2;
  c1_rho_re = 0.0;
  c1_rho_im = 0.0;
  c1_anorm = c1_eml_matlab_zlanhs(chartInstance, c1_A, c1_ilo, c1_ihi);
  c1_a = 2.2204460492503131E-16 * c1_anorm;
  c1_atol = 2.2250738585072014E-308;
  if (c1_a > 2.2250738585072014E-308) {
    c1_atol = c1_a;
  }

  c1_a = 2.2250738585072014E-308;
  if (c1_anorm > 2.2250738585072014E-308) {
    c1_a = c1_anorm;
  }

  c1_ascale = 1.0 / c1_a;
  c1_failed = TRUE;
  for (c1_j = c1_ihi + 1; c1_j < 8; c1_j++) {
    c1_alpha1[sf_mex_lw_bounds_check(c1_j, 1, 7) - 1] = c1_A
      [(sf_mex_lw_bounds_check(c1_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c1_j, 1,
          7) - 1)) - 1];
  }

  c1_guard1 = FALSE;
  c1_guard2 = FALSE;
  if (c1_ihi >= c1_ilo) {
    c1_ifirst = c1_ilo;
    c1_istart = c1_ilo;
    c1_ilast = c1_ihi;
    c1_ilastm1 = c1_ihi - 2;
    c1_ifrstm = c1_ilo;
    c1_iiter = 0;
    c1_maxit = 30 * ((c1_ihi - c1_ilo) + 1);
    c1_goto60 = FALSE;
    c1_goto70 = FALSE;
    c1_goto90 = FALSE;
    c1_jiter = 1;
    do {
      c1_exitg1 = 0U;
      if (c1_jiter <= c1_maxit) {
        if (c1_ilast == c1_ilo) {
          c1_goto60 = TRUE;
        } else {
          sf_mex_lw_bounds_check(c1_ilastm1 + 1, 1, 7);
          sf_mex_lw_bounds_check(c1_ilast, 1, 7);
          if (muDoubleScalarAbs(c1_b_A[(c1_ilast + 7 * c1_ilastm1) - 1].re) +
              muDoubleScalarAbs(c1_b_A[(c1_ilast + 7 * c1_ilastm1) - 1].im) <=
              c1_atol) {
            c1_b_A[(c1_ilast + 7 * c1_ilastm1) - 1] = c1_dc2;
            c1_goto60 = TRUE;
          } else {
            c1_j = c1_ilastm1;
            c1_exitg3 = 0U;
            while ((c1_exitg3 == 0U) && (c1_j + 1 >= c1_ilo)) {
              c1_jm1 = c1_j - 1;
              if (c1_j + 1 == c1_ilo) {
                c1_ilazro = TRUE;
              } else {
                sf_mex_lw_bounds_check(c1_jm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c1_j + 1, 1, 7);
                if (muDoubleScalarAbs(c1_b_A[c1_j + 7 * c1_jm1].re) +
                    muDoubleScalarAbs(c1_b_A[c1_j + 7 * c1_jm1].im) <= c1_atol)
                {
                  c1_b_A[c1_j + 7 * c1_jm1] = c1_dc2;
                  c1_ilazro = TRUE;
                } else {
                  c1_ilazro = FALSE;
                }
              }

              if (c1_ilazro) {
                c1_ifirst = c1_j + 1;
                c1_goto70 = TRUE;
                c1_exitg3 = 1U;
              } else {
                c1_j = c1_jm1;
              }
            }
          }
        }

        if (c1_goto60 || c1_goto70) {
          c1_ilazro = TRUE;
        } else {
          c1_ilazro = FALSE;
        }

        if (!c1_ilazro) {
          for (c1_jm1 = 0; c1_jm1 < 7; c1_jm1++) {
            c1_alpha1[c1_jm1] = c1_dc1;
            c1_beta1[c1_jm1] = c1_dc1;
          }

          *c1_info = -1.0;
          c1_exitg1 = 1U;
        } else {
          c1_b_guard1 = FALSE;
          if (c1_goto60) {
            c1_goto60 = FALSE;
            c1_alpha1[sf_mex_lw_bounds_check(c1_ilast, 1, 7) - 1] = c1_b_A
              [(sf_mex_lw_bounds_check(c1_ilast, 1, 7) + 7 *
                (sf_mex_lw_bounds_check(c1_ilast, 1, 7) - 1)) - 1];
            c1_ilast = c1_ilastm1 + 1;
            c1_ilastm1--;
            if (c1_ilast < c1_ilo) {
              c1_failed = FALSE;
              c1_guard2 = TRUE;
              c1_exitg1 = 1U;
            } else {
              c1_iiter = 0;
              c1_eshift_re = 0.0;
              c1_eshift_im = 0.0;
              if (c1_ifrstm > c1_ilast) {
                c1_ifrstm = c1_ilo;
              }

              c1_b_guard1 = TRUE;
            }
          } else {
            if (c1_goto70) {
              c1_goto70 = FALSE;
              c1_iiter++;
              c1_ifrstm = c1_ifirst;
              if (c1_iiter - c1_div_s32_floor(chartInstance, c1_iiter, 10) * 10
                  != 0) {
                sf_mex_lw_bounds_check(c1_ilastm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c1_ilast, 1, 7);
                c1_c_A.re = -(c1_b_A[(c1_ilast + 7 * (c1_ilast - 1)) - 1].re -
                              c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].re);
                c1_c_A.im = -(c1_b_A[(c1_ilast + 7 * (c1_ilast - 1)) - 1].im -
                              c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].im);
                c1_t1 = c1_eml_div(chartInstance, c1_c_A, 2.0);
                c1_d.re = (c1_t1.re * c1_t1.re - c1_t1.im * c1_t1.im) +
                  (c1_b_A[c1_ilastm1 + 7 * (c1_ilast - 1)].re * c1_b_A[(c1_ilast
                    + 7 * c1_ilastm1) - 1].re - c1_b_A[c1_ilastm1 + 7 *
                   (c1_ilast - 1)].im * c1_b_A[(c1_ilast + 7 * c1_ilastm1) - 1].
                   im);
                c1_d.im = (c1_t1.re * c1_t1.im + c1_t1.im * c1_t1.re) +
                  (c1_b_A[c1_ilastm1 + 7 * (c1_ilast - 1)].re * c1_b_A[(c1_ilast
                    + 7 * c1_ilastm1) - 1].im + c1_b_A[c1_ilastm1 + 7 *
                   (c1_ilast - 1)].im * c1_b_A[(c1_ilast + 7 * c1_ilastm1) - 1].
                   re);
                c1_b_sqrt(chartInstance, &c1_d);
                c1_sigma1.re = c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].re -
                  (c1_t1.re - c1_d.re);
                c1_sigma1.im = c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].im -
                  (c1_t1.im - c1_d.im);
                c1_sigma2_re = c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].re -
                  (c1_t1.re + c1_d.re);
                c1_sigma2_im = c1_b_A[c1_ilastm1 + 7 * c1_ilastm1].im -
                  (c1_t1.im + c1_d.im);
                c1_a = muDoubleScalarAbs(c1_sigma1.re - c1_b_A[(c1_ilast + 7 *
                  (c1_ilast - 1)) - 1].re);
                c1_anorm = muDoubleScalarAbs(c1_sigma1.im - c1_b_A[(c1_ilast + 7
                  * (c1_ilast - 1)) - 1].im);
                if (c1_a < c1_anorm) {
                  c1_a /= c1_anorm;
                  c1_anorm *= muDoubleScalarSqrt(c1_a * c1_a + 1.0);
                } else if (c1_a > c1_anorm) {
                  c1_anorm /= c1_a;
                  c1_anorm = muDoubleScalarSqrt(c1_anorm * c1_anorm + 1.0) *
                    c1_a;
                } else {
                  if (!muDoubleScalarIsNaN(c1_anorm)) {
                    c1_anorm = c1_a * 1.4142135623730951;
                  }
                }

                c1_a = muDoubleScalarAbs(c1_sigma2_re - c1_b_A[(c1_ilast + 7 *
                  (c1_ilast - 1)) - 1].re);
                c1_b = muDoubleScalarAbs(c1_sigma2_im - c1_b_A[(c1_ilast + 7 *
                  (c1_ilast - 1)) - 1].im);
                if (c1_a < c1_b) {
                  c1_a /= c1_b;
                  c1_b *= muDoubleScalarSqrt(c1_a * c1_a + 1.0);
                } else if (c1_a > c1_b) {
                  c1_b /= c1_a;
                  c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_a;
                } else {
                  if (!muDoubleScalarIsNaN(c1_b)) {
                    c1_b = c1_a * 1.4142135623730951;
                  }
                }

                if (c1_anorm <= c1_b) {
                  c1_sigma2_re = c1_sigma1.re;
                  c1_sigma2_im = c1_sigma1.im;
                  c1_rho_re = c1_t1.re - c1_d.re;
                  c1_rho_im = c1_t1.im - c1_d.im;
                } else {
                  c1_rho_re = c1_t1.re + c1_d.re;
                  c1_rho_im = c1_t1.im + c1_d.im;
                }
              } else {
                c1_eshift_re += c1_b_A[(sf_mex_lw_bounds_check(c1_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c1_ilastm1 + 1, 1, 7) - 1)) - 1].
                  re;
                c1_eshift_im += c1_b_A[(sf_mex_lw_bounds_check(c1_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c1_ilastm1 + 1, 1, 7) - 1)) - 1].
                  im;
                c1_sigma2_re = c1_eshift_re;
                c1_sigma2_im = c1_eshift_im;
              }

              c1_j = c1_ilastm1;
              c1_jp1 = c1_ilastm1 + 1;
              c1_exitg2 = 0U;
              while ((c1_exitg2 == 0U) && (c1_j + 1 > c1_ifirst)) {
                c1_jm1 = c1_j - 1;
                c1_istart = c1_j + 1;
                c1_ctemp.re = c1_b_A[c1_j + 7 * c1_j].re - c1_sigma2_re;
                c1_ctemp.im = c1_b_A[c1_j + 7 * c1_j].im - c1_sigma2_im;
                c1_anorm = c1_ascale * (muDoubleScalarAbs(c1_ctemp.re) +
                  muDoubleScalarAbs(c1_ctemp.im));
                sf_mex_lw_bounds_check(c1_jp1 + 1, 1, 7);
                c1_a = c1_ascale * (muDoubleScalarAbs(c1_b_A[c1_jp1 + 7 * c1_j].
                  re) + muDoubleScalarAbs(c1_b_A[c1_jp1 + 7 * c1_j].im));
                c1_b = c1_anorm;
                if (c1_a > c1_anorm) {
                  c1_b = c1_a;
                }

                if ((c1_b < 1.0) && (c1_b != 0.0)) {
                  c1_anorm /= c1_b;
                  c1_a /= c1_b;
                }

                sf_mex_lw_bounds_check(c1_jm1 + 1, 1, 7);
                if ((muDoubleScalarAbs(c1_b_A[c1_j + 7 * c1_jm1].re) +
                     muDoubleScalarAbs(c1_b_A[c1_j + 7 * c1_jm1].im)) * c1_a <=
                    c1_anorm * c1_atol) {
                  c1_goto90 = TRUE;
                  c1_exitg2 = 1U;
                } else {
                  c1_jp1 = c1_j;
                  c1_j = c1_jm1;
                }
              }

              if (!c1_goto90) {
                c1_istart = c1_ifirst;
                if (c1_ifirst == c1_ilastm1 + 1) {
                  c1_ctemp.re = c1_rho_re;
                  c1_ctemp.im = c1_rho_im;
                } else {
                  c1_ctemp.re = c1_b_A[(sf_mex_lw_bounds_check(c1_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c1_ifirst,
                    1, 7) - 1)) - 1].re - c1_sigma2_re;
                  c1_ctemp.im = c1_b_A[(sf_mex_lw_bounds_check(c1_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c1_ifirst,
                    1, 7) - 1)) - 1].im - c1_sigma2_im;
                }

                c1_goto90 = TRUE;
              }
            }

            if (c1_goto90) {
              c1_goto90 = FALSE;
              sf_mex_lw_bounds_check(c1_istart, 1, 7);
              sf_mex_lw_bounds_check(c1_istart + 1, 1, 7);
              c1_b_eml_matlab_zlartg(chartInstance, c1_ctemp, c1_b_A[c1_istart +
                7 * (c1_istart - 1)], &c1_a, &c1_sigma1);
              c1_j = c1_istart - 1;
              c1_jm1 = c1_istart - 2;
              while (c1_j + 1 < c1_ilast) {
                c1_jp1 = c1_j + 1;
                if (c1_j + 1 > c1_istart) {
                  c1_c_A = c1_b_A[(sf_mex_lw_bounds_check(c1_j + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c1_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c1_d_A = c1_b_A[(sf_mex_lw_bounds_check(c1_jp1 + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c1_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c1_eml_matlab_zlartg(chartInstance, c1_c_A, c1_d_A, &c1_a,
                                       &c1_sigma1, &c1_t1);
                  c1_b_A[c1_j + 7 * c1_jm1] = c1_t1;
                  c1_b_A[c1_jp1 + 7 * c1_jm1] = c1_dc2;
                }

                for (c1_jm1 = c1_j; c1_jm1 + 1 <= c1_ilast; c1_jm1++) {
                  sf_mex_lw_bounds_check(c1_jm1 + 1, 1, 7);
                  c1_t1.re = c1_a * c1_b_A[c1_j + 7 * c1_jm1].re;
                  c1_t1.im = c1_a * c1_b_A[c1_j + 7 * c1_jm1].im;
                  sf_mex_lw_bounds_check(c1_jp1 + 1, 1, 7);
                  c1_d.re = c1_sigma1.re * c1_b_A[c1_jp1 + 7 * c1_jm1].re -
                    c1_sigma1.im * c1_b_A[c1_jp1 + 7 * c1_jm1].im;
                  c1_d.im = c1_sigma1.re * c1_b_A[c1_jp1 + 7 * c1_jm1].im +
                    c1_sigma1.im * c1_b_A[c1_jp1 + 7 * c1_jm1].re;
                  c1_c_A = c1_b_A[c1_j + 7 * c1_jm1];
                  c1_d_A = c1_b_A[c1_j + 7 * c1_jm1];
                  c1_b_A[c1_jp1 + 7 * c1_jm1].re = c1_a * c1_b_A[c1_jp1 + 7 *
                    c1_jm1].re - (c1_sigma1.re * c1_b_A[c1_j + 7 * c1_jm1].re +
                                  c1_sigma1.im * c1_b_A[c1_j + 7 * c1_jm1].im);
                  c1_b_A[c1_jp1 + 7 * c1_jm1].im = c1_a * c1_b_A[c1_jp1 + 7 *
                    c1_jm1].im - (c1_sigma1.re * c1_c_A.im - c1_sigma1.im *
                                  c1_d_A.re);
                  c1_b_A[c1_j + 7 * c1_jm1].re = c1_t1.re + c1_d.re;
                  c1_b_A[c1_j + 7 * c1_jm1].im = c1_t1.im + c1_d.im;
                }

                c1_sigma1.re = -c1_sigma1.re;
                c1_sigma1.im = -c1_sigma1.im;
                c1_jm1 = c1_jp1 + 2;
                if (c1_ilast < c1_jm1) {
                  c1_jm1 = c1_ilast;
                }

                for (c1_i = c1_ifrstm - 1; c1_i + 1 <= c1_jm1; c1_i++) {
                  sf_mex_lw_bounds_check(c1_jp1 + 1, 1, 7);
                  sf_mex_lw_bounds_check(c1_i + 1, 1, 7);
                  c1_t1.re = c1_a * c1_b_A[c1_i + 7 * c1_jp1].re;
                  c1_t1.im = c1_a * c1_b_A[c1_i + 7 * c1_jp1].im;
                  c1_d.re = c1_sigma1.re * c1_b_A[c1_i + 7 * c1_j].re -
                    c1_sigma1.im * c1_b_A[c1_i + 7 * c1_j].im;
                  c1_d.im = c1_sigma1.re * c1_b_A[c1_i + 7 * c1_j].im +
                    c1_sigma1.im * c1_b_A[c1_i + 7 * c1_j].re;
                  c1_c_A = c1_b_A[c1_i + 7 * c1_jp1];
                  c1_d_A = c1_b_A[c1_i + 7 * c1_jp1];
                  c1_b_A[c1_i + 7 * c1_j].re = c1_a * c1_b_A[c1_i + 7 * c1_j].re
                    - (c1_sigma1.re * c1_b_A[c1_i + 7 * c1_jp1].re +
                       c1_sigma1.im * c1_b_A[c1_i + 7 * c1_jp1].im);
                  c1_b_A[c1_i + 7 * c1_j].im = c1_a * c1_b_A[c1_i + 7 * c1_j].im
                    - (c1_sigma1.re * c1_c_A.im - c1_sigma1.im * c1_d_A.re);
                  c1_b_A[c1_i + 7 * c1_jp1].re = c1_t1.re + c1_d.re;
                  c1_b_A[c1_i + 7 * c1_jp1].im = c1_t1.im + c1_d.im;
                }

                c1_jm1 = c1_j;
                c1_j = c1_jp1;
              }
            }

            c1_b_guard1 = TRUE;
          }

          if (c1_b_guard1 == TRUE) {
            c1_jiter++;
          }
        }
      } else {
        c1_guard2 = TRUE;
        c1_exitg1 = 1U;
      }
    } while (c1_exitg1 == 0U);
  } else {
    c1_guard1 = TRUE;
  }

  if (c1_guard2 == TRUE) {
    if (c1_failed) {
      *c1_info = (real_T)c1_ilast;
      for (c1_jm1 = 1; c1_jm1 <= c1_ilast; c1_jm1++) {
        c1_alpha1[sf_mex_lw_bounds_check(c1_jm1, 1, 7) - 1] = c1_dc1;
        c1_beta1[c1_jm1 - 1] = c1_dc1;
      }
    } else {
      c1_guard1 = TRUE;
    }
  }

  if (c1_guard1 == TRUE) {
    for (c1_j = 1; c1_j <= c1_ilo - 1; c1_j++) {
      c1_alpha1[sf_mex_lw_bounds_check(c1_j, 1, 7) - 1] = c1_b_A
        [(sf_mex_lw_bounds_check(c1_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c1_j,
            1, 7) - 1)) - 1];
    }

    *c1_info = 0.0;
  }
}

static real_T c1_eml_matlab_zlanhs(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, creal_T c1_A[49], int32_T c1_ilo, int32_T c1_ihi)
{
  real_T c1_f;
  real_T c1_scale;
  real_T c1_sumsq;
  boolean_T c1_firstNonZero;
  int32_T c1_j;
  int32_T c1_c;
  int32_T c1_i;
  real_T c1_temp1;
  real_T c1_temp2;
  c1_f = 0.0;
  if (!(c1_ilo > c1_ihi)) {
    c1_scale = 0.0;
    c1_sumsq = 0.0;
    c1_firstNonZero = TRUE;
    for (c1_j = c1_ilo; c1_j <= c1_ihi; c1_j++) {
      c1_c = c1_j + 1;
      if (c1_ihi < c1_c) {
        c1_c = c1_ihi;
      }

      for (c1_i = c1_ilo; c1_i <= c1_c; c1_i++) {
        sf_mex_lw_bounds_check(c1_i, 1, 7);
        sf_mex_lw_bounds_check(c1_j, 1, 7);
        if (c1_A[(c1_i + 7 * (c1_j - 1)) - 1].re != 0.0) {
          c1_temp1 = muDoubleScalarAbs(c1_A[(c1_i + 7 * (c1_j - 1)) - 1].re);
          if (c1_firstNonZero) {
            c1_sumsq = 1.0;
            c1_scale = c1_temp1;
            c1_firstNonZero = FALSE;
          } else if (c1_scale < c1_temp1) {
            c1_temp2 = c1_scale / c1_temp1;
            c1_sumsq = 1.0 + c1_sumsq * c1_temp2 * c1_temp2;
            c1_scale = c1_temp1;
          } else {
            c1_temp2 = c1_temp1 / c1_scale;
            c1_sumsq += c1_temp2 * c1_temp2;
          }
        }

        if (c1_A[(c1_i + 7 * (c1_j - 1)) - 1].im != 0.0) {
          c1_temp1 = muDoubleScalarAbs(c1_A[(c1_i + 7 * (c1_j - 1)) - 1].im);
          if (c1_firstNonZero) {
            c1_sumsq = 1.0;
            c1_scale = c1_temp1;
            c1_firstNonZero = FALSE;
          } else if (c1_scale < c1_temp1) {
            c1_temp2 = c1_scale / c1_temp1;
            c1_sumsq = 1.0 + c1_sumsq * c1_temp2 * c1_temp2;
            c1_scale = c1_temp1;
          } else {
            c1_temp2 = c1_temp1 / c1_scale;
            c1_sumsq += c1_temp2 * c1_temp2;
          }
        }
      }
    }

    if (c1_sumsq < 0.0) {
      c1_eml_error(chartInstance);
    }

    c1_f = c1_scale * muDoubleScalarSqrt(c1_sumsq);
  }

  return c1_f;
}

static creal_T c1_eml_div(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_x, real_T c1_y)
{
  creal_T c1_z;
  if (c1_x.im == 0.0) {
    c1_z.re = c1_x.re / c1_y;
    c1_z.im = 0.0;
  } else if (c1_x.re == 0.0) {
    c1_z.re = 0.0;
    c1_z.im = c1_x.im / c1_y;
  } else {
    c1_z.re = c1_x.re / c1_y;
    c1_z.im = c1_x.im / c1_y;
  }

  return c1_z;
}

static void c1_b_eml_matlab_zlartg(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn)
{
  real_T c1_scale;
  real_T c1_b;
  real_T c1_f2;
  real_T c1_fs_re;
  real_T c1_fs_im;
  real_T c1_gs_re;
  real_T c1_gs_im;
  boolean_T c1_guard1 = FALSE;
  real_T c1_b_b;
  real_T c1_g2s;
  c1_scale = muDoubleScalarAbs(c1_f.re);
  c1_b = muDoubleScalarAbs(c1_f.im);
  if (c1_b > c1_scale) {
    c1_scale = c1_b;
  }

  c1_b = muDoubleScalarAbs(c1_g.re);
  c1_f2 = muDoubleScalarAbs(c1_g.im);
  if (c1_f2 > c1_b) {
    c1_b = c1_f2;
  }

  if (c1_b > c1_scale) {
    c1_scale = c1_b;
  }

  c1_fs_re = c1_f.re;
  c1_fs_im = c1_f.im;
  c1_gs_re = c1_g.re;
  c1_gs_im = c1_g.im;
  c1_guard1 = FALSE;
  if (c1_scale >= 7.4428285367870146E+137) {
    do {
      c1_fs_re *= 1.3435752215134178E-138;
      c1_fs_im *= 1.3435752215134178E-138;
      c1_gs_re *= 1.3435752215134178E-138;
      c1_gs_im *= 1.3435752215134178E-138;
      c1_scale *= 1.3435752215134178E-138;
    } while (!(c1_scale < 7.4428285367870146E+137));

    c1_guard1 = TRUE;
  } else if (c1_scale <= 1.3435752215134178E-138) {
    if ((c1_g.re == 0.0) && (c1_g.im == 0.0)) {
      *c1_cs = 1.0;
      c1_sn->re = 0.0;
      c1_sn->im = 0.0;
    } else {
      do {
        c1_fs_re *= 7.4428285367870146E+137;
        c1_fs_im *= 7.4428285367870146E+137;
        c1_gs_re *= 7.4428285367870146E+137;
        c1_gs_im *= 7.4428285367870146E+137;
        c1_scale *= 7.4428285367870146E+137;
      } while (!(c1_scale > 1.3435752215134178E-138));

      c1_guard1 = TRUE;
    }
  } else {
    c1_guard1 = TRUE;
  }

  if (c1_guard1 == TRUE) {
    c1_f2 = c1_fs_re * c1_fs_re + c1_fs_im * c1_fs_im;
    c1_scale = c1_gs_re * c1_gs_re + c1_gs_im * c1_gs_im;
    c1_b = c1_scale;
    if (1.0 > c1_scale) {
      c1_b = 1.0;
    }

    if (c1_f2 <= c1_b * 2.0041683600089728E-292) {
      if ((c1_f.re == 0.0) && (c1_f.im == 0.0)) {
        *c1_cs = 0.0;
        c1_f2 = muDoubleScalarAbs(c1_gs_re);
        c1_b_b = muDoubleScalarAbs(c1_gs_im);
        if (c1_f2 < c1_b_b) {
          c1_f2 /= c1_b_b;
          c1_b_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
        } else if (c1_f2 > c1_b_b) {
          c1_b_b /= c1_f2;
          c1_b_b = muDoubleScalarSqrt(c1_b_b * c1_b_b + 1.0) * c1_f2;
        } else {
          if (!muDoubleScalarIsNaN(c1_b_b)) {
            c1_b_b = c1_f2 * 1.4142135623730951;
          }
        }

        c1_sn->re = c1_gs_re / c1_b_b;
        c1_sn->im = -c1_gs_im / c1_b_b;
      } else {
        c1_f2 = muDoubleScalarAbs(c1_fs_re);
        c1_b = muDoubleScalarAbs(c1_fs_im);
        if (c1_f2 < c1_b) {
          c1_f2 /= c1_b;
          c1_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
        } else if (c1_f2 > c1_b) {
          c1_b /= c1_f2;
          c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_f2;
        } else {
          if (!muDoubleScalarIsNaN(c1_b)) {
            c1_b = c1_f2 * 1.4142135623730951;
          }
        }

        if (c1_scale < 0.0) {
          c1_eml_error(chartInstance);
        }

        c1_g2s = muDoubleScalarSqrt(c1_scale);
        *c1_cs = c1_b / c1_g2s;
        c1_b = muDoubleScalarAbs(c1_f.re);
        c1_f2 = muDoubleScalarAbs(c1_f.im);
        if (c1_f2 > c1_b) {
          c1_b = c1_f2;
        }

        if (c1_b > 1.0) {
          c1_f2 = muDoubleScalarAbs(c1_f.re);
          c1_b_b = muDoubleScalarAbs(c1_f.im);
          if (c1_f2 < c1_b_b) {
            c1_f2 /= c1_b_b;
            c1_b_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
          } else if (c1_f2 > c1_b_b) {
            c1_b_b /= c1_f2;
            c1_b_b = muDoubleScalarSqrt(c1_b_b * c1_b_b + 1.0) * c1_f2;
          } else {
            if (!muDoubleScalarIsNaN(c1_b_b)) {
              c1_b_b = c1_f2 * 1.4142135623730951;
            }
          }

          c1_fs_re = c1_f.re / c1_b_b;
          c1_fs_im = c1_f.im / c1_b_b;
        } else {
          c1_b = 7.4428285367870146E+137 * c1_f.re;
          c1_scale = 7.4428285367870146E+137 * c1_f.im;
          c1_f2 = muDoubleScalarAbs(c1_b);
          c1_b_b = muDoubleScalarAbs(c1_scale);
          if (c1_f2 < c1_b_b) {
            c1_f2 /= c1_b_b;
            c1_b_b *= muDoubleScalarSqrt(c1_f2 * c1_f2 + 1.0);
          } else if (c1_f2 > c1_b_b) {
            c1_b_b /= c1_f2;
            c1_b_b = muDoubleScalarSqrt(c1_b_b * c1_b_b + 1.0) * c1_f2;
          } else {
            if (!muDoubleScalarIsNaN(c1_b_b)) {
              c1_b_b = c1_f2 * 1.4142135623730951;
            }
          }

          c1_fs_re = c1_b / c1_b_b;
          c1_fs_im = c1_scale / c1_b_b;
        }

        c1_gs_re /= c1_g2s;
        c1_gs_im = -c1_gs_im / c1_g2s;
        c1_sn->re = c1_fs_re * c1_gs_re - c1_fs_im * c1_gs_im;
        c1_sn->im = c1_fs_re * c1_gs_im + c1_fs_im * c1_gs_re;
      }
    } else {
      c1_b = 1.0 + c1_scale / c1_f2;
      if (c1_b < 0.0) {
        c1_eml_error(chartInstance);
      }

      c1_b = muDoubleScalarSqrt(c1_b);
      *c1_cs = 1.0 / c1_b;
      c1_b_b = c1_f2 + c1_scale;
      c1_fs_re = c1_b * c1_fs_re / c1_b_b;
      c1_fs_im = c1_b * c1_fs_im / c1_b_b;
      c1_sn->re = c1_fs_re * c1_gs_re - c1_fs_im * -c1_gs_im;
      c1_sn->im = c1_fs_re * -c1_gs_im + c1_fs_im * c1_gs_re;
    }
  }
}

static void c1_b_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i11;
  static char_T c1_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c1_u[26];
  const mxArray *c1_y = NULL;
  for (c1_i11 = 0; c1_i11 < 26; c1_i11++) {
    c1_u[c1_i11] = c1_varargin_1[c1_i11];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 26), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static void c1_c_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i12;
  static char_T c1_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  for (c1_i12 = 0; c1_i12 < 34; c1_i12++) {
    c1_u[c1_i12] = c1_varargin_1[c1_i12];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static void c1_d_eye(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T c1_I[49])
{
  int32_T c1_i;
  int32_T c1_b_i;
  for (c1_i = 0; c1_i < 49; c1_i++) {
    c1_I[c1_i] = 0.0;
  }

  c1_i = 0;
  for (c1_b_i = 0; c1_b_i < 7; c1_b_i++) {
    c1_I[c1_i] = 1.0;
    c1_i += 8;
  }
}

static void c1_mldivide(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_A[49], real_T c1_B[7], creal_T c1_Y[7])
{
  int32_T c1_i13;
  creal_T c1_b_A[49];
  int8_T c1_ipiv[7];
  int32_T c1_info;
  int32_T c1_j;
  int32_T c1_jj;
  int32_T c1_jp1j;
  int32_T c1_iy;
  int32_T c1_ix;
  real_T c1_smax;
  int32_T c1_b_k;
  real_T c1_s;
  int32_T c1_jrow;
  creal_T c1_temp;
  int32_T c1_jA;
  for (c1_i13 = 0; c1_i13 < 49; c1_i13++) {
    c1_b_A[c1_i13] = c1_A[c1_i13];
  }

  for (c1_i13 = 0; c1_i13 < 7; c1_i13++) {
    c1_ipiv[c1_i13] = (int8_T)(1 + c1_i13);
  }

  c1_info = 0;
  for (c1_j = 0; c1_j < 6; c1_j++) {
    c1_jj = c1_j << 3;
    c1_jp1j = c1_jj + 2;
    c1_iy = 1;
    c1_ix = c1_jj;
    c1_smax = muDoubleScalarAbs(c1_b_A[c1_jj].re) + muDoubleScalarAbs
      (c1_b_A[c1_jj].im);
    for (c1_b_k = 2; c1_b_k <= 7 - c1_j; c1_b_k++) {
      c1_ix++;
      sf_mex_lw_bounds_check(c1_ix + 1, 1, 49);
      c1_s = muDoubleScalarAbs(c1_b_A[c1_ix].re) + muDoubleScalarAbs
        (c1_b_A[c1_ix].im);
      if (c1_s > c1_smax) {
        c1_iy = c1_b_k;
        c1_smax = c1_s;
      }
    }

    if ((c1_b_A[(c1_jj + c1_iy) - 1].re != 0.0) || (c1_b_A[(c1_jj + c1_iy) - 1].
         im != 0.0)) {
      if (c1_iy - 1 != 0) {
        c1_ipiv[c1_j] = (int8_T)(c1_j + c1_iy);
        c1_jrow = 1 + c1_j;
        c1_iy = (c1_jrow + c1_iy) - 1;
        for (c1_b_k = 0; c1_b_k < 7; c1_b_k++) {
          c1_temp = c1_b_A[sf_mex_lw_bounds_check(c1_jrow, 1, 49) - 1];
          c1_b_A[c1_jrow - 1] = c1_b_A[sf_mex_lw_bounds_check(c1_iy, 1, 49) - 1];
          c1_b_A[c1_iy - 1] = c1_temp;
          c1_jrow += 7;
          c1_iy += 7;
        }
      }

      c1_i13 = (c1_jp1j - c1_j) + 5;
      for (c1_jrow = c1_jp1j; c1_jrow <= c1_i13; c1_jrow++) {
        c1_b_A[c1_jrow - 1] = c1_b_eml_div(chartInstance, c1_b_A[c1_jrow - 1],
          c1_b_A[c1_jj]);
      }
    } else {
      c1_info = c1_j + 1;
    }

    c1_jA = c1_jj + 8;
    c1_iy = c1_jj + 7;
    for (c1_jrow = 1; c1_jrow <= 6 - c1_j; c1_jrow++) {
      sf_mex_lw_bounds_check(c1_iy + 1, 1, 49);
      if ((c1_b_A[c1_iy].re != 0.0) || (c1_b_A[c1_iy].im != 0.0)) {
        c1_temp.re = -c1_b_A[c1_iy].re - c1_b_A[c1_iy].im * 0.0;
        c1_temp.im = c1_b_A[c1_iy].re * 0.0 + -c1_b_A[c1_iy].im;
        c1_ix = c1_jp1j;
        c1_i13 = (c1_jA - c1_j) + 6;
        for (c1_b_k = 1 + c1_jA; c1_b_k <= c1_i13; c1_b_k++) {
          c1_smax = c1_b_A[sf_mex_lw_bounds_check(c1_ix, 1, 49) - 1].re *
            c1_temp.re - c1_b_A[sf_mex_lw_bounds_check(c1_ix, 1, 49) - 1].im *
            c1_temp.im;
          c1_s = c1_b_A[sf_mex_lw_bounds_check(c1_ix, 1, 49) - 1].re *
            c1_temp.im + c1_b_A[sf_mex_lw_bounds_check(c1_ix, 1, 49) - 1].im *
            c1_temp.re;
          c1_b_A[sf_mex_lw_bounds_check(c1_b_k, 1, 49) - 1].re =
            c1_b_A[sf_mex_lw_bounds_check(c1_b_k, 1, 49) - 1].re + c1_smax;
          c1_b_A[sf_mex_lw_bounds_check(c1_b_k, 1, 49) - 1].im =
            c1_b_A[sf_mex_lw_bounds_check(c1_b_k, 1, 49) - 1].im + c1_s;
          c1_ix++;
        }
      }

      c1_iy += 7;
      c1_jA += 7;
    }
  }

  if ((c1_info == 0) && (!((c1_b_A[48].re != 0.0) || (c1_b_A[48].im != 0.0)))) {
    c1_info = 7;
  }

  if (c1_info > 0) {
    c1_d_eml_warning(chartInstance);
  }

  for (c1_i13 = 0; c1_i13 < 7; c1_i13++) {
    c1_Y[c1_i13].re = c1_B[c1_i13];
    c1_Y[c1_i13].im = 0.0;
  }

  for (c1_jrow = 0; c1_jrow < 7; c1_jrow++) {
    if (c1_ipiv[c1_jrow] != c1_jrow + 1) {
      c1_temp = c1_Y[c1_jrow];
      c1_Y[c1_jrow] = c1_Y[sf_mex_lw_bounds_check((int32_T)c1_ipiv[c1_jrow], 1,
        7) - 1];
      c1_Y[sf_mex_lw_bounds_check((int32_T)c1_ipiv[c1_jrow], 1, 7) - 1] =
        c1_temp;
    }
  }

  for (c1_b_k = 0; c1_b_k < 7; c1_b_k++) {
    c1_iy = 7 * c1_b_k;
    if ((c1_Y[c1_b_k].re != 0.0) || (c1_Y[c1_b_k].im != 0.0)) {
      for (c1_jrow = c1_b_k + 2; c1_jrow < 8; c1_jrow++) {
        c1_smax = c1_Y[c1_b_k].re * c1_b_A[(c1_jrow + c1_iy) - 1].im +
          c1_Y[c1_b_k].im * c1_b_A[(c1_jrow + c1_iy) - 1].re;
        c1_Y[c1_jrow - 1].re -= c1_Y[c1_b_k].re * c1_b_A[(c1_jrow + c1_iy) - 1].
          re - c1_Y[c1_b_k].im * c1_b_A[(c1_jrow + c1_iy) - 1].im;
        c1_Y[c1_jrow - 1].im -= c1_smax;
      }
    }
  }

  for (c1_b_k = 6; c1_b_k > -1; c1_b_k += -1) {
    c1_iy = 7 * c1_b_k;
    if ((c1_Y[c1_b_k].re != 0.0) || (c1_Y[c1_b_k].im != 0.0)) {
      c1_Y[c1_b_k] = c1_b_eml_div(chartInstance, c1_Y[c1_b_k], c1_b_A[c1_b_k +
        c1_iy]);
      for (c1_jrow = 0; c1_jrow + 1 <= c1_b_k; c1_jrow++) {
        c1_smax = c1_Y[c1_b_k].re * c1_b_A[c1_jrow + c1_iy].im + c1_Y[c1_b_k].im
          * c1_b_A[c1_jrow + c1_iy].re;
        c1_Y[c1_jrow].re -= c1_Y[c1_b_k].re * c1_b_A[c1_jrow + c1_iy].re -
          c1_Y[c1_b_k].im * c1_b_A[c1_jrow + c1_iy].im;
        c1_Y[c1_jrow].im -= c1_smax;
      }
    }
  }
}

static creal_T c1_b_eml_div(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, creal_T c1_x, creal_T c1_y)
{
  creal_T c1_z;
  real_T c1_brm;
  real_T c1_bim;
  real_T c1_d;
  if (c1_y.im == 0.0) {
    if (c1_x.im == 0.0) {
      c1_z.re = c1_x.re / c1_y.re;
      c1_z.im = 0.0;
    } else if (c1_x.re == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = c1_x.im / c1_y.re;
    } else {
      c1_z.re = c1_x.re / c1_y.re;
      c1_z.im = c1_x.im / c1_y.re;
    }
  } else if (c1_y.re == 0.0) {
    if (c1_x.re == 0.0) {
      c1_z.re = c1_x.im / c1_y.im;
      c1_z.im = 0.0;
    } else if (c1_x.im == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = -(c1_x.re / c1_y.im);
    } else {
      c1_z.re = c1_x.im / c1_y.im;
      c1_z.im = -(c1_x.re / c1_y.im);
    }
  } else {
    c1_brm = muDoubleScalarAbs(c1_y.re);
    c1_bim = muDoubleScalarAbs(c1_y.im);
    if (c1_brm > c1_bim) {
      c1_bim = c1_y.im / c1_y.re;
      c1_d = c1_y.re + c1_bim * c1_y.im;
      c1_z.re = (c1_x.re + c1_bim * c1_x.im) / c1_d;
      c1_z.im = (c1_x.im - c1_bim * c1_x.re) / c1_d;
    } else if (c1_bim == c1_brm) {
      c1_bim = c1_y.re > 0.0 ? 0.5 : -0.5;
      c1_d = c1_y.im > 0.0 ? 0.5 : -0.5;
      c1_z.re = (c1_x.re * c1_bim + c1_x.im * c1_d) / c1_brm;
      c1_z.im = (c1_x.im * c1_bim - c1_x.re * c1_d) / c1_brm;
    } else {
      c1_bim = c1_y.re / c1_y.im;
      c1_d = c1_y.im + c1_bim * c1_y.re;
      c1_z.re = (c1_bim * c1_x.re + c1_x.im) / c1_d;
      c1_z.im = (c1_bim * c1_x.im - c1_x.re) / c1_d;
    }
  }

  return c1_z;
}

static void c1_d_eml_warning(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i14;
  static char_T c1_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c1_u[27];
  const mxArray *c1_y = NULL;
  for (c1_i14 = 0; c1_i14 < 27; c1_i14++) {
    c1_u[c1_i14] = c1_varargin_1[c1_i14];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static void c1_abs(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   creal_T c1_x[2], real_T c1_y[2])
{
  int32_T c1_b_k;
  real_T c1_a;
  real_T c1_b;
  for (c1_b_k = 0; c1_b_k < 2; c1_b_k++) {
    c1_a = muDoubleScalarAbs(c1_x[c1_b_k].re);
    c1_b = muDoubleScalarAbs(c1_x[c1_b_k].im);
    if (c1_a < c1_b) {
      c1_a /= c1_b;
      c1_b *= muDoubleScalarSqrt(c1_a * c1_a + 1.0);
    } else if (c1_a > c1_b) {
      c1_b /= c1_a;
      c1_b = muDoubleScalarSqrt(c1_b * c1_b + 1.0) * c1_a;
    } else {
      if (!muDoubleScalarIsNaN(c1_b)) {
        c1_b = c1_a * 1.4142135623730951;
      }
    }

    c1_y[c1_b_k] = c1_b;
  }
}

static void c1_b_eml_error(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
  int32_T c1_i15;
  static char_T c1_varargin_1[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'l', 'o', 'g', '1', '0', '_', 'd', 'o', 'm',
    'a', 'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c1_u[31];
  const mxArray *c1_y = NULL;
  for (c1_i15 = 0; c1_i15 < 31; c1_i15++) {
    c1_u[c1_i15] = c1_varargin_1[c1_i15];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 31), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c1_y));
}

static real_T c1_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_ERR, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_ERR), &c1_thisId);
  sf_mex_destroy(&c1_ERR);
  return c1_y;
}

static real_T c1_b_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d1;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d1, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d1;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_c_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_Wn, const char_T *c1_identifier, real_T
  c1_y[7])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Wn), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_Wn);
}

static void c1_d_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7])
{
  real_T c1_dv22[7];
  int32_T c1_i16;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv22, 1, 0, 0U, 1, 0U, 1, 7);
  for (c1_i16 = 0; c1_i16 < 7; c1_i16++) {
    c1_y[c1_i16] = c1_dv22[c1_i16];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_e_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_ABK, const char_T *c1_identifier, real_T
  c1_y[70])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_ABK), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_ABK);
}

static void c1_f_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[70])
{
  real_T c1_dv23[70];
  int32_T c1_i17;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_ABK_not_empty = FALSE;
  } else {
    chartInstance->c1_ABK_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv23, 1, 0, 0U, 1, 0U, 2, 7,
                  10);
    for (c1_i17 = 0; c1_i17 < 70; c1_i17++) {
      c1_y[c1_i17] = c1_dv23[c1_i17];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_g_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_CD, const char_T *c1_identifier, real_T
  c1_y[16])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_CD), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_CD);
}

static void c1_h_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[16])
{
  real_T c1_dv24[16];
  int32_T c1_i18;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_CD_not_empty = FALSE;
  } else {
    chartInstance->c1_CD_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv24, 1, 0, 0U, 1, 0U, 2, 2,
                  8);
    for (c1_i18 = 0; c1_i18 < 16; c1_i18++) {
      c1_y[c1_i18] = c1_dv24[c1_i18];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_i_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Glk, const char_T *c1_identifier, real_T
  c1_y[459])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Glk), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Glk);
}

static void c1_j_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[459])
{
  real_T c1_dv25[459];
  int32_T c1_i19;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Glk_not_empty = FALSE;
  } else {
    chartInstance->c1_Glk_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv25, 1, 0, 0U, 1, 0U, 2,
                  153, 3);
    for (c1_i19 = 0; c1_i19 < 459; c1_i19++) {
      c1_y[c1_i19] = c1_dv25[c1_i19];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_k_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_LK, const char_T *c1_identifier, real_T
  c1_y[15000])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_LK), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_LK);
}

static void c1_l_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[15000])
{
  real_T c1_dv26[15000];
  int32_T c1_i20;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_LK_not_empty = FALSE;
  } else {
    chartInstance->c1_LK_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv26, 1, 0, 0U, 1, 0U, 2,
                  100, 150);
    for (c1_i20 = 0; c1_i20 < 15000; c1_i20++) {
      c1_y[c1_i20] = c1_dv26[c1_i20];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_m_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Pabk, const char_T *c1_identifier, real_T
  c1_y[100])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Pabk), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Pabk);
}

static void c1_n_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[100])
{
  real_T c1_dv27[100];
  int32_T c1_i21;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Pabk_not_empty = FALSE;
  } else {
    chartInstance->c1_Pabk_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv27, 1, 0, 0U, 1, 0U, 2, 10,
                  10);
    for (c1_i21 = 0; c1_i21 < 100; c1_i21++) {
      c1_y[c1_i21] = c1_dv27[c1_i21];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_o_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Pcd, const char_T *c1_identifier, real_T
  c1_y[64])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Pcd), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Pcd);
}

static void c1_p_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[64])
{
  real_T c1_dv28[64];
  int32_T c1_i22;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Pcd_not_empty = FALSE;
  } else {
    chartInstance->c1_Pcd_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv28, 1, 0, 0U, 1, 0U, 2, 8,
                  8);
    for (c1_i22 = 0; c1_i22 < 64; c1_i22++) {
      c1_y[c1_i22] = c1_dv28[c1_i22];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_q_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_Plk, const char_T *c1_identifier, real_T
  c1_y[1099])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Plk), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Plk);
}

static void c1_r_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[1099])
{
  real_T c1_dv29[1099];
  int32_T c1_i23;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Plk_not_empty = FALSE;
  } else {
    chartInstance->c1_Plk_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv29, 1, 0, 0U, 1, 0U, 2,
                  157, 7);
    for (c1_i23 = 0; c1_i23 < 1099; c1_i23++) {
      c1_y[c1_i23] = c1_dv29[c1_i23];
    }
  }

  sf_mex_destroy(&c1_u);
}

static real_T c1_s_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_U1, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_U1), &c1_thisId);
  sf_mex_destroy(&c1_b_U1);
  return c1_y;
}

static real_T c1_t_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d2;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_U1_not_empty = FALSE;
  } else {
    chartInstance->c1_U1_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d2, 1, 0, 0U, 0, 0U, 0);
    c1_y = c1_d2;
  }

  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_u_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_VARX, const char_T *c1_identifier, real_T
  c1_y[306])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_v_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_VARX), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_VARX);
}

static void c1_v_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[306])
{
  real_T c1_dv30[306];
  int32_T c1_i24;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_VARX_not_empty = FALSE;
  } else {
    chartInstance->c1_VARX_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv30, 1, 0, 0U, 1, 0U, 2, 2,
                  153);
    for (c1_i24 = 0; c1_i24 < 306; c1_i24++) {
      c1_y[c1_i24] = c1_dv30[c1_i24];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_w_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_VARX1, const char_T *c1_identifier, real_T
  c1_y[300])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_x_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_VARX1), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_VARX1);
}

static void c1_x_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[300])
{
  real_T c1_dv31[300];
  int32_T c1_i25;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_VARX1_not_empty = FALSE;
  } else {
    chartInstance->c1_VARX1_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv31, 1, 0, 0U, 1, 0U, 2, 2,
                  150);
    for (c1_i25 = 0; c1_i25 < 300; c1_i25++) {
      c1_y[c1_i25] = c1_dv31[c1_i25];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_y_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, const mxArray *c1_b_X, const char_T *c1_identifier, real_T
  c1_y[7])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_ab_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_X), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_X);
}

static void c1_ab_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7])
{
  real_T c1_dv32[7];
  int32_T c1_i26;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_X_not_empty = FALSE;
  } else {
    chartInstance->c1_X_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv32, 1, 0, 0U, 1, 0U, 1, 7);
    for (c1_i26 = 0; c1_i26 < 7; c1_i26++) {
      c1_y[c1_i26] = c1_dv32[c1_i26];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_bb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_X1, const char_T *c1_identifier, real_T
  c1_y[7])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_cb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_X1), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_X1);
}

static void c1_cb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[7])
{
  real_T c1_dv33[7];
  int32_T c1_i27;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_X1_not_empty = FALSE;
  } else {
    chartInstance->c1_X1_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv33, 1, 0, 0U, 1, 0U, 1, 7);
    for (c1_i27 = 0; c1_i27 < 7; c1_i27++) {
      c1_y[c1_i27] = c1_dv33[c1_i27];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_db_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Y1, const char_T *c1_identifier, real_T
  c1_y[2])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_eb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Y1), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Y1);
}

static void c1_eb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[2])
{
  real_T c1_dv34[2];
  int32_T c1_i28;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Y1_not_empty = FALSE;
  } else {
    chartInstance->c1_Y1_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv34, 1, 0, 0U, 1, 0U, 1, 2);
    for (c1_i28 = 0; c1_i28 < 2; c1_i28++) {
      c1_y[c1_i28] = c1_dv34[c1_i28];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_fb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Ylk, const char_T *c1_identifier, real_T
  c1_y[3])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_gb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Ylk), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Ylk);
}

static void c1_gb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3])
{
  real_T c1_dv35[3];
  int32_T c1_i29;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Ylk_not_empty = FALSE;
  } else {
    chartInstance->c1_Ylk_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv35, 1, 0, 0U, 1, 0U, 2, 1,
                  3);
    for (c1_i29 = 0; c1_i29 < 3; c1_i29++) {
      c1_y[c1_i29] = c1_dv35[c1_i29];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_hb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_Z, const char_T *c1_identifier, real_T
  c1_y[156])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_ib_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Z), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Z);
}

static void c1_ib_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[156])
{
  real_T c1_dv36[156];
  int32_T c1_i30;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_Z_not_empty = FALSE;
  } else {
    chartInstance->c1_Z_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv36, 1, 0, 0U, 1, 0U, 1,
                  156);
    for (c1_i30 = 0; c1_i30 < 156; c1_i30++) {
      c1_y[c1_i30] = c1_dv36[c1_i30];
    }
  }

  sf_mex_destroy(&c1_u);
}

static void c1_jb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_eta, const char_T *c1_identifier, real_T
  c1_y[468])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_kb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_eta), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_eta);
}

static void c1_kb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[468])
{
  real_T c1_dv37[468];
  int32_T c1_i31;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_eta_not_empty = FALSE;
  } else {
    chartInstance->c1_eta_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv37, 1, 0, 0U, 1, 0U, 2,
                  156, 3);
    for (c1_i31 = 0; c1_i31 < 468; c1_i31++) {
      c1_y[c1_i31] = c1_dv37[c1_i31];
    }
  }

  sf_mex_destroy(&c1_u);
}

static real_T c1_lb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_k, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_mb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_k), &c1_thisId);
  sf_mex_destroy(&c1_b_k);
  return c1_y;
}

static real_T c1_mb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d3;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_k_not_empty = FALSE;
  } else {
    chartInstance->c1_k_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d3, 1, 0, 0U, 0, 0U, 0);
    c1_y = c1_d3;
  }

  sf_mex_destroy(&c1_u);
  return c1_y;
}

static real_T c1_nb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_saw, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_ob_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_saw), &c1_thisId);
  sf_mex_destroy(&c1_b_saw);
  return c1_y;
}

static real_T c1_ob_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d4;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_saw_not_empty = FALSE;
  } else {
    chartInstance->c1_saw_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d4, 1, 0, 0U, 0, 0U, 0);
    c1_y = c1_d4;
  }

  sf_mex_destroy(&c1_u);
  return c1_y;
}

static real_T c1_pb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_start, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_qb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_start),
    &c1_thisId);
  sf_mex_destroy(&c1_b_start);
  return c1_y;
}

static real_T c1_qb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d5;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_start_not_empty = FALSE;
  } else {
    chartInstance->c1_start_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d5, 1, 0, 0U, 0, 0U, 0);
    c1_y = c1_d5;
  }

  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_rb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_b_w, const char_T *c1_identifier, real_T
  c1_y[300])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_sb_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_w), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_w);
}

static void c1_sb_emlrt_marshallIn(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *
  chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[300])
{
  real_T c1_dv38[300];
  int32_T c1_i32;
  if (mxIsEmpty(c1_u)) {
    chartInstance->c1_w_not_empty = FALSE;
  } else {
    chartInstance->c1_w_not_empty = TRUE;
    sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv38, 1, 0, 0U, 1, 0U, 2, 1,
                  300);
    for (c1_i32 = 0; c1_i32 < 300; c1_i32++) {
      c1_y[c1_i32] = c1_dv38[c1_i32];
    }
  }

  sf_mex_destroy(&c1_u);
}

static uint8_T c1_tb_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_b_is_active_c1_sim1_rlti_varx_fast_gust2, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_ub_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_sim1_rlti_varx_fast_gust2), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_sim1_rlti_varx_fast_gust2);
  return c1_y;
}

static uint8_T c1_ub_emlrt_marshallIn
  (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance, const mxArray
   *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_sqrt(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                    real_T *c1_x)
{
  if (*c1_x < 0.0) {
    c1_eml_error(chartInstance);
  }

  *c1_x = muDoubleScalarSqrt(*c1_x);
}

static real_T c1_fastQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[156], real_T c1_u[3], real_T c1_y[2], real_T
  c1_theta[306], real_T c1_M[1099], real_T c1_lambda, real_T c1_tol, real_T
  c1_b_eta[468], real_T c1_ireg, real_T c1_reg, real_T c1_G[459], real_T c1_Y[3],
  real_T *c1_b_saw)
{
  real_T c1_ERR;
  real_T c1_c;
  real_T c1_rho;
  int32_T c1_i33;
  real_T c1_b_M[3];
  real_T c1_dv39[9];
  int32_T c1_m;
  int32_T c1_n;
  int32_T c1_b_k;
  real_T c1_c_eta[468];
  real_T c1_G1[153];
  real_T c1_b_Y1;
  int32_T c1_j;
  int32_T c1_exitg3;
  real_T c1_b_z[156];
  real_T c1_c_M[936];
  int32_T c1_lda;
  int32_T c1_ldb;
  int32_T c1_ldc;
  char_T c1_TRANSA;
  char_T c1_TRANSB;
  real_T c1_b_y[6];
  real_T c1_b_G[153];
  int32_T c1_exitg2;
  boolean_T c1_guard2 = FALSE;
  real_T c1_x1[157];
  real_T c1_y1[157];
  real_T c1_x2[157];
  real_T c1_dv40[936];
  int32_T c1_exitg1;
  boolean_T c1_guard1 = FALSE;
  real_T c1_c_y[2];
  if (c1_M[0] == 0.0) {
    c1_M[0] = 1.0;
    c1_c = muDoubleScalarPower(c1_lambda, 51.0);
    c1_sqrt(chartInstance, &c1_c);
    c1_rho = c1_ireg;
    c1_sqrt(chartInstance, &c1_rho);
    c1_c = c1_mrdivide(chartInstance, c1_c, c1_rho);
    for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
      c1_b_M[c1_i33] = c1_c;
    }

    c1_b_diag(chartInstance, c1_b_M, c1_dv39);
    c1_i33 = 0;
    c1_m = 0;
    for (c1_n = 0; c1_n < 3; c1_n++) {
      for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
        c1_M[(c1_b_k + c1_i33) + 158] = c1_dv39[c1_b_k + c1_m];
      }

      c1_i33 += 157;
      c1_m += 3;
    }

    c1_c = c1_ireg;
    c1_sqrt(chartInstance, &c1_c);
    c1_c = c1_mrdivide(chartInstance, 1.0, c1_c);
    for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
      c1_b_M[c1_i33] = c1_c;
    }

    c1_b_diag(chartInstance, c1_b_M, c1_dv39);
    c1_i33 = 0;
    c1_m = 0;
    for (c1_n = 0; c1_n < 3; c1_n++) {
      for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
        c1_M[(c1_b_k + c1_i33) + 782] = c1_dv39[c1_b_k + c1_m];
      }

      c1_i33 += 157;
      c1_m += 3;
    }
  }

  if (*c1_b_saw > 50.5) {
    *c1_b_saw = 1.0;
    c1_c = (1.0 - muDoubleScalarPower(c1_lambda, 51.0)) * c1_reg;
    c1_sqrt(chartInstance, &c1_c);
    c1_c_eye(chartInstance, c1_dv39);
    c1_i33 = 0;
    for (c1_m = 0; c1_m < 3; c1_m++) {
      for (c1_n = 0; c1_n < 153; c1_n++) {
        c1_c_eta[c1_n + c1_i33] = c1_b_eta[(c1_n + c1_i33) + 3];
      }

      c1_i33 += 156;
    }

    c1_i33 = 0;
    c1_m = 0;
    for (c1_n = 0; c1_n < 3; c1_n++) {
      for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
        c1_c_eta[(c1_b_k + c1_i33) + 153] = c1_c * c1_dv39[c1_b_k + c1_m];
      }

      c1_i33 += 156;
      c1_m += 3;
    }

    c1_i33 = 0;
    for (c1_m = 0; c1_m < 3; c1_m++) {
      for (c1_n = 0; c1_n < 156; c1_n++) {
        c1_b_eta[c1_n + c1_i33] = c1_c_eta[c1_n + c1_i33];
      }

      c1_i33 += 156;
    }
  } else {
    (*c1_b_saw)++;
    c1_i33 = 0;
    for (c1_m = 0; c1_m < 3; c1_m++) {
      for (c1_n = 0; c1_n < 153; c1_n++) {
        c1_c_eta[c1_n + c1_i33] = c1_b_eta[(c1_n + c1_i33) + 3];
      }

      c1_i33 += 156;
    }

    c1_i33 = 0;
    for (c1_m = 0; c1_m < 3; c1_m++) {
      for (c1_n = 0; c1_n < 3; c1_n++) {
        c1_c_eta[(c1_n + c1_i33) + 153] = 0.0;
      }

      c1_i33 += 156;
    }

    c1_i33 = 0;
    for (c1_m = 0; c1_m < 3; c1_m++) {
      for (c1_n = 0; c1_n < 156; c1_n++) {
        c1_b_eta[c1_n + c1_i33] = c1_c_eta[c1_n + c1_i33];
      }

      c1_i33 += 156;
    }
  }

  for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
    c1_G1[c1_i33] = c1_M[c1_i33 + 4];
  }

  c1_b_Y1 = c1_M[0];
  c1_j = 0;
  do {
    c1_exitg3 = 0U;
    if (c1_j < 3) {
      c1_M[0] = c1_Y[c1_j];
      for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
        c1_b_z[c1_i33] = c1_b_eta[c1_i33 + 156 * c1_j];
      }

      c1_i33 = 0;
      c1_m = 0;
      for (c1_n = 0; c1_n < 6; c1_n++) {
        for (c1_b_k = 0; c1_b_k < 156; c1_b_k++) {
          c1_c_M[c1_b_k + c1_i33] = c1_M[(c1_b_k + c1_m) + 158];
        }

        c1_i33 += 156;
        c1_m += 157;
      }

      c1_m = 1;
      c1_n = 6;
      c1_b_k = 156;
      c1_rho = 1.0;
      c1_lda = 1;
      c1_ldb = 156;
      c1_c = 0.0;
      c1_ldc = 1;
      c1_TRANSA = 'N';
      c1_TRANSB = 'N';
      for (c1_i33 = 0; c1_i33 < 6; c1_i33++) {
        c1_b_y[c1_i33] = 0.0;
      }

      dgemm32(&c1_TRANSA, &c1_TRANSB, &c1_m, &c1_n, &c1_b_k, &c1_rho, &c1_b_z[0],
              &c1_lda, &c1_c_M[0], &c1_ldb, &c1_c, &c1_b_y[0], &c1_ldc);
      c1_i33 = 0;
      for (c1_m = 0; c1_m < 6; c1_m++) {
        c1_M[c1_i33 + 157] = c1_b_y[c1_m];
        c1_i33 += 157;
      }

      for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
        c1_b_G[c1_i33] = c1_G[c1_i33 + 153 * c1_j];
      }

      for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
        c1_M[c1_i33 + 1] = c1_b_G[c1_i33];
      }

      for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
        c1_M[c1_i33 + 154] = 0.0;
      }

      c1_m = 1;
      do {
        c1_exitg2 = 0U;
        if (c1_m - 1 < 6) {
          c1_guard2 = FALSE;
          if (1 + c1_m > 4) {
            if (c1_M[157 * c1_m] != 0.0) {
              c1_rho = c1_mrdivide(chartInstance, c1_M[157 * c1_m], c1_M[0]);
              c1_c = 1.0 + muDoubleScalarPower(muDoubleScalarAbs(c1_rho), 2.0);
              c1_sqrt(chartInstance, &c1_c);
              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_x1[c1_i33] = c1_M[c1_i33] + c1_rho * c1_M[c1_i33 + 157 * c1_m];
              }

              c1_c_rdivide(chartInstance, c1_x1, c1_c, c1_y1);
              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_M[c1_i33] = c1_y1[c1_i33];
              }

              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_x1[c1_i33] = -c1_rho * c1_M[c1_i33] + c1_c * c1_M[c1_i33 +
                  157 * c1_m];
              }

              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_M[c1_i33 + 157 * c1_m] = c1_x1[c1_i33];
              }
            }

            c1_guard2 = TRUE;
          } else if (muDoubleScalarAbs(c1_M[c1_m + 157 * c1_m]) >
                     muDoubleScalarAbs(c1_M[c1_m])) {
            if (c1_M[c1_m] != 0.0) {
              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_x1[c1_i33] = c1_M[c1_i33] - c1_M[c1_i33 + 157 * c1_m];
              }

              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_y1[c1_i33] = c1_M[c1_i33] + c1_M[c1_i33 + 157 * c1_m];
              }

              c1_c = c1_mrdivide(chartInstance, -c1_y1[c1_m], c1_x1[c1_m]);
              c1_sqrt(chartInstance, &c1_c);
              c1_c *= 0.5;
              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_x2[c1_i33] = c1_c * c1_x1[c1_i33];
              }

              c1_c = c1_mrdivide(chartInstance, -c1_x1[c1_m], c1_y1[c1_m]);
              c1_sqrt(chartInstance, &c1_c);
              c1_c *= 0.5;
              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_y1[c1_i33] *= c1_c;
              }

              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_M[c1_i33] = c1_x2[c1_i33] + c1_y1[c1_i33];
              }

              for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
                c1_M[c1_i33 + 157 * c1_m] = c1_y1[c1_i33] - c1_x2[c1_i33];
              }
            }

            c1_guard2 = TRUE;
          } else {
            c1_i33 = 0;
            for (c1_m = 0; c1_m < 7; c1_m++) {
              for (c1_n = 0; c1_n < 157; c1_n++) {
                c1_M[c1_n + c1_i33] = 0.0;
              }

              c1_i33 += 157;
            }

            c1_M[0] = 1.0;
            c1_c = muDoubleScalarPower(c1_lambda, 51.0);
            c1_sqrt(chartInstance, &c1_c);
            c1_rho = c1_ireg;
            c1_sqrt(chartInstance, &c1_rho);
            c1_c = c1_mrdivide(chartInstance, c1_c, c1_rho);
            for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
              c1_b_M[c1_i33] = c1_c;
            }

            c1_b_diag(chartInstance, c1_b_M, c1_dv39);
            c1_i33 = 0;
            c1_m = 0;
            for (c1_n = 0; c1_n < 3; c1_n++) {
              for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
                c1_M[(c1_b_k + c1_i33) + 158] = c1_dv39[c1_b_k + c1_m];
              }

              c1_i33 += 157;
              c1_m += 3;
            }

            c1_c = c1_ireg;
            c1_sqrt(chartInstance, &c1_c);
            c1_c = c1_mrdivide(chartInstance, 1.0, c1_c);
            for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
              c1_b_M[c1_i33] = c1_c;
            }

            c1_b_diag(chartInstance, c1_b_M, c1_dv39);
            c1_i33 = 0;
            c1_m = 0;
            for (c1_n = 0; c1_n < 3; c1_n++) {
              for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
                c1_M[(c1_b_k + c1_i33) + 782] = c1_dv39[c1_b_k + c1_m];
              }

              c1_i33 += 157;
              c1_m += 3;
            }

            c1_i33 = 0;
            for (c1_m = 0; c1_m < 3; c1_m++) {
              for (c1_n = 0; c1_n < 153; c1_n++) {
                c1_G[c1_n + c1_i33] = 0.0;
              }

              c1_i33 += 153;
            }

            for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
              c1_Y[c1_i33] = 1.0;
            }

            c1_i33 = 0;
            for (c1_m = 0; c1_m < 3; c1_m++) {
              c1_b_M[c1_m] = c1_M[c1_i33 + 157];
              c1_i33 += 157;
            }

            c1_ERR = c1_norm(chartInstance, c1_b_M);
            for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
              c1_z[c1_i33] = 0.0;
            }

            c1_exitg2 = 1U;
          }

          if (c1_guard2 == TRUE) {
            c1_m++;
            sf_mex_listen_for_ctrl_c(chartInstance->S);
          }
        } else {
          c1_Y[c1_j] = c1_M[0];
          for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
            c1_b_G[c1_i33] = c1_M[c1_i33 + 4];
          }

          for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
            c1_G[c1_i33 + 153 * c1_j] = c1_b_G[c1_i33];
          }

          c1_j++;
          sf_mex_listen_for_ctrl_c(chartInstance->S);
          c1_exitg2 = 2U;
        }
      } while (c1_exitg2 == 0U);

      if (c1_exitg2 == 1U) {
        c1_exitg3 = 1U;
      }
    } else {
      for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
        c1_b_z[c1_i33] = c1_z[c1_i33 + 3];
      }

      for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
        c1_b_z[c1_i33 + 153] = c1_u[c1_i33];
      }

      for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
        c1_z[c1_i33] = c1_b_z[c1_i33];
      }

      c1_M[0] = c1_b_Y1;
      c1_c = c1_lambda;
      c1_sqrt(chartInstance, &c1_c);
      c1_i33 = 0;
      c1_m = 0;
      for (c1_n = 0; c1_n < 6; c1_n++) {
        for (c1_b_k = 0; c1_b_k < 156; c1_b_k++) {
          c1_c_M[c1_b_k + c1_i33] = c1_M[(c1_b_k + c1_m) + 158];
        }

        c1_i33 += 156;
        c1_m += 157;
      }

      c1_d_rdivide(chartInstance, c1_c_M, c1_c, c1_dv40);
      c1_i33 = 0;
      c1_m = 0;
      for (c1_n = 0; c1_n < 6; c1_n++) {
        for (c1_b_k = 0; c1_b_k < 156; c1_b_k++) {
          c1_M[(c1_b_k + c1_i33) + 158] = c1_dv40[c1_b_k + c1_m];
        }

        c1_i33 += 157;
        c1_m += 156;
      }

      for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
        c1_b_z[c1_i33] = c1_z[c1_i33];
      }

      c1_i33 = 0;
      c1_m = 0;
      for (c1_n = 0; c1_n < 6; c1_n++) {
        for (c1_b_k = 0; c1_b_k < 156; c1_b_k++) {
          c1_c_M[c1_b_k + c1_i33] = c1_M[(c1_b_k + c1_m) + 158];
        }

        c1_i33 += 156;
        c1_m += 157;
      }

      c1_m = 1;
      c1_n = 6;
      c1_b_k = 156;
      c1_rho = 1.0;
      c1_lda = 1;
      c1_ldb = 156;
      c1_c = 0.0;
      c1_ldc = 1;
      c1_TRANSA = 'N';
      c1_TRANSB = 'N';
      for (c1_i33 = 0; c1_i33 < 6; c1_i33++) {
        c1_b_y[c1_i33] = 0.0;
      }

      dgemm32(&c1_TRANSA, &c1_TRANSB, &c1_m, &c1_n, &c1_b_k, &c1_rho, &c1_b_z[0],
              &c1_lda, &c1_c_M[0], &c1_ldb, &c1_c, &c1_b_y[0], &c1_ldc);
      c1_i33 = 0;
      for (c1_m = 0; c1_m < 6; c1_m++) {
        c1_M[c1_i33 + 157] = c1_b_y[c1_m];
        c1_i33 += 157;
      }

      for (c1_i33 = 0; c1_i33 < 153; c1_i33++) {
        c1_M[c1_i33 + 1] = c1_G1[c1_i33];
      }

      for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
        c1_M[c1_i33 + 154] = 0.0;
      }

      c1_m = 1;
      c1_exitg3 = 2U;
    }
  } while (c1_exitg3 == 0U);

  if (c1_exitg3 == 1U) {
  } else {
    do {
      c1_exitg1 = 0U;
      if (c1_m - 1 < 6) {
        c1_guard1 = FALSE;
        if (1 + c1_m > 4) {
          if (c1_M[157 * c1_m] != 0.0) {
            c1_rho = c1_mrdivide(chartInstance, c1_M[157 * c1_m], c1_M[0]);
            c1_c = 1.0 + muDoubleScalarPower(muDoubleScalarAbs(c1_rho), 2.0);
            c1_sqrt(chartInstance, &c1_c);
            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_x1[c1_i33] = c1_M[c1_i33] + c1_rho * c1_M[c1_i33 + 157 * c1_m];
            }

            c1_c_rdivide(chartInstance, c1_x1, c1_c, c1_y1);
            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_M[c1_i33] = c1_y1[c1_i33];
            }

            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_x1[c1_i33] = -c1_rho * c1_M[c1_i33] + c1_c * c1_M[c1_i33 + 157 *
                c1_m];
            }

            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_M[c1_i33 + 157 * c1_m] = c1_x1[c1_i33];
            }
          }

          c1_guard1 = TRUE;
        } else if (muDoubleScalarAbs(c1_M[c1_m + 157 * c1_m]) >
                   muDoubleScalarAbs(c1_M[c1_m])) {
          if (c1_M[c1_m] != 0.0) {
            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_x1[c1_i33] = c1_M[c1_i33] - c1_M[c1_i33 + 157 * c1_m];
            }

            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_y1[c1_i33] = c1_M[c1_i33] + c1_M[c1_i33 + 157 * c1_m];
            }

            c1_c = c1_mrdivide(chartInstance, -c1_y1[c1_m], c1_x1[c1_m]);
            c1_sqrt(chartInstance, &c1_c);
            c1_c *= 0.5;
            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_x2[c1_i33] = c1_c * c1_x1[c1_i33];
            }

            c1_c = c1_mrdivide(chartInstance, -c1_x1[c1_m], c1_y1[c1_m]);
            c1_sqrt(chartInstance, &c1_c);
            c1_c *= 0.5;
            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_y1[c1_i33] *= c1_c;
            }

            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_M[c1_i33] = c1_x2[c1_i33] + c1_y1[c1_i33];
            }

            for (c1_i33 = 0; c1_i33 < 157; c1_i33++) {
              c1_M[c1_i33 + 157 * c1_m] = c1_y1[c1_i33] - c1_x2[c1_i33];
            }
          }

          c1_guard1 = TRUE;
        } else {
          c1_i33 = 0;
          for (c1_m = 0; c1_m < 7; c1_m++) {
            for (c1_n = 0; c1_n < 157; c1_n++) {
              c1_M[c1_n + c1_i33] = 0.0;
            }

            c1_i33 += 157;
          }

          c1_M[0] = 1.0;
          c1_c = muDoubleScalarPower(c1_lambda, 51.0);
          c1_sqrt(chartInstance, &c1_c);
          c1_rho = c1_ireg;
          c1_sqrt(chartInstance, &c1_rho);
          c1_c = c1_mrdivide(chartInstance, c1_c, c1_rho);
          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_b_M[c1_i33] = c1_c;
          }

          c1_b_diag(chartInstance, c1_b_M, c1_dv39);
          c1_i33 = 0;
          c1_m = 0;
          for (c1_n = 0; c1_n < 3; c1_n++) {
            for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
              c1_M[(c1_b_k + c1_i33) + 158] = c1_dv39[c1_b_k + c1_m];
            }

            c1_i33 += 157;
            c1_m += 3;
          }

          c1_c = c1_ireg;
          c1_sqrt(chartInstance, &c1_c);
          c1_c = c1_mrdivide(chartInstance, 1.0, c1_c);
          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_b_M[c1_i33] = c1_c;
          }

          c1_b_diag(chartInstance, c1_b_M, c1_dv39);
          c1_i33 = 0;
          c1_m = 0;
          for (c1_n = 0; c1_n < 3; c1_n++) {
            for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
              c1_M[(c1_b_k + c1_i33) + 782] = c1_dv39[c1_b_k + c1_m];
            }

            c1_i33 += 157;
            c1_m += 3;
          }

          c1_i33 = 0;
          for (c1_m = 0; c1_m < 3; c1_m++) {
            for (c1_n = 0; c1_n < 153; c1_n++) {
              c1_G[c1_n + c1_i33] = 0.0;
            }

            c1_i33 += 153;
          }

          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_Y[c1_i33] = 1.0;
          }

          c1_i33 = 0;
          for (c1_m = 0; c1_m < 3; c1_m++) {
            c1_b_M[c1_m] = c1_M[c1_i33 + 157];
            c1_i33 += 157;
          }

          c1_ERR = c1_norm(chartInstance, c1_b_M);
          for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
            c1_z[c1_i33] = 0.0;
          }

          c1_exitg1 = 1U;
        }

        if (c1_guard1 == TRUE) {
          c1_m++;
          sf_mex_listen_for_ctrl_c(chartInstance->S);
        }
      } else {
        c1_i33 = 0;
        for (c1_m = 0; c1_m < 3; c1_m++) {
          c1_b_M[c1_m] = c1_M[c1_i33 + 157];
          c1_i33 += 157;
        }

        c1_ERR = c1_norm(chartInstance, c1_b_M);
        if (c1_ERR < c1_tol) {
          c1_b_mrdivide(chartInstance, *(real_T (*)[153])&c1_M[4], c1_M[0],
                        c1_G1);
          for (c1_i33 = 0; c1_i33 < 2; c1_i33++) {
            c1_c = 0.0;
            c1_m = 0;
            for (c1_n = 0; c1_n < 153; c1_n++) {
              c1_c += c1_theta[c1_m + c1_i33] * c1_z[c1_n + 3];
              c1_m += 2;
            }

            c1_c_y[c1_i33] = c1_y[c1_i33] - c1_c;
          }

          c1_i33 = 0;
          for (c1_m = 0; c1_m < 153; c1_m++) {
            for (c1_n = 0; c1_n < 2; c1_n++) {
              c1_theta[c1_n + c1_i33] += c1_c_y[c1_n] * c1_G1[c1_m];
            }

            c1_i33 += 2;
          }
        } else {
          c1_i33 = 0;
          for (c1_m = 0; c1_m < 7; c1_m++) {
            for (c1_n = 0; c1_n < 157; c1_n++) {
              c1_M[c1_n + c1_i33] = 0.0;
            }

            c1_i33 += 157;
          }

          c1_M[0] = 1.0;
          c1_c = muDoubleScalarPower(c1_lambda, 51.0);
          c1_sqrt(chartInstance, &c1_c);
          c1_rho = c1_ireg;
          c1_sqrt(chartInstance, &c1_rho);
          c1_c = c1_mrdivide(chartInstance, c1_c, c1_rho);
          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_b_M[c1_i33] = c1_c;
          }

          c1_b_diag(chartInstance, c1_b_M, c1_dv39);
          c1_i33 = 0;
          c1_m = 0;
          for (c1_n = 0; c1_n < 3; c1_n++) {
            for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
              c1_M[(c1_b_k + c1_i33) + 158] = c1_dv39[c1_b_k + c1_m];
            }

            c1_i33 += 157;
            c1_m += 3;
          }

          c1_c = c1_ireg;
          c1_sqrt(chartInstance, &c1_c);
          c1_c = c1_mrdivide(chartInstance, 1.0, c1_c);
          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_b_M[c1_i33] = c1_c;
          }

          c1_b_diag(chartInstance, c1_b_M, c1_dv39);
          c1_i33 = 0;
          c1_m = 0;
          for (c1_n = 0; c1_n < 3; c1_n++) {
            for (c1_b_k = 0; c1_b_k < 3; c1_b_k++) {
              c1_M[(c1_b_k + c1_i33) + 782] = c1_dv39[c1_b_k + c1_m];
            }

            c1_i33 += 157;
            c1_m += 3;
          }

          c1_i33 = 0;
          for (c1_m = 0; c1_m < 3; c1_m++) {
            for (c1_n = 0; c1_n < 153; c1_n++) {
              c1_G[c1_n + c1_i33] = 0.0;
            }

            c1_i33 += 153;
          }

          for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
            c1_Y[c1_i33] = 1.0;
          }

          for (c1_i33 = 0; c1_i33 < 156; c1_i33++) {
            c1_z[c1_i33] = 0.0;
          }
        }

        c1_exitg1 = 1U;
      }
    } while (c1_exitg1 == 0U);
  }

  return c1_ERR;
}

static void c1_inverseQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[8], real_T c1_y[2], real_T c1_theta[16], real_T
  c1_P[64], real_T c1_lambda)
{
  real_T c1_rho;
  real_T c1_c;
  int32_T c1_i34;
  int32_T c1_i;
  real_T c1_d6;
  int32_T c1_i35;
  real_T c1_b_z[8];
  real_T c1_M[81];
  int32_T c1_i36;
  real_T c1_b_y[2];
  if (c1_lambda < 0.0) {
    c1_eml_error(chartInstance);
  }

  c1_rho = muDoubleScalarSqrt(c1_lambda);
  if (c1_lambda < 0.0) {
    c1_eml_error(chartInstance);
  }

  c1_c = muDoubleScalarSqrt(c1_lambda);
  c1_i34 = 0;
  for (c1_i = 0; c1_i < 8; c1_i++) {
    c1_d6 = 0.0;
    for (c1_i35 = 0; c1_i35 < 8; c1_i35++) {
      c1_d6 += c1_z[c1_i35] * c1_P[c1_i35 + c1_i34];
    }

    c1_b_z[c1_i] = c1_d6 / c1_rho;
    c1_i34 += 8;
  }

  c1_M[0] = 1.0;
  c1_i34 = 0;
  for (c1_i = 0; c1_i < 8; c1_i++) {
    c1_M[c1_i34 + 9] = c1_b_z[c1_i];
    c1_i34 += 9;
  }

  for (c1_i34 = 0; c1_i34 < 8; c1_i34++) {
    c1_M[c1_i34 + 1] = 0.0;
  }

  c1_i34 = 0;
  c1_i = 0;
  for (c1_i35 = 0; c1_i35 < 8; c1_i35++) {
    for (c1_i36 = 0; c1_i36 < 8; c1_i36++) {
      c1_M[(c1_i36 + c1_i34) + 10] = c1_P[c1_i36 + c1_i] / c1_c;
    }

    c1_i34 += 9;
    c1_i += 8;
  }

  for (c1_i = 0; c1_i < 8; c1_i++) {
    if (!(c1_M[9 * (8 - c1_i)] == 0.0)) {
      if (c1_M[0] == 0.0) {
        for (c1_i34 = 0; c1_i34 < 9; c1_i34++) {
          c1_rho = c1_M[c1_i34 + 9 * (8 - c1_i)];
          c1_M[c1_i34 + 9 * (8 - c1_i)] = -c1_M[c1_i34];
          c1_M[c1_i34] = c1_rho;
        }
      } else {
        c1_rho = c1_M[9 * (8 - c1_i)] / c1_M[0];
        c1_c = muDoubleScalarSqrt(1.0 + muDoubleScalarPower(muDoubleScalarAbs
          (c1_rho), 2.0));
        for (c1_i34 = 0; c1_i34 < 9; c1_i34++) {
          c1_M[c1_i34] = (c1_M[c1_i34] + c1_rho * c1_M[c1_i34 + 9 * (8 - c1_i)])
            / c1_c;
          c1_M[c1_i34 + 9 * (8 - c1_i)] = -c1_rho * c1_M[c1_i34] + c1_c *
            c1_M[c1_i34 + 9 * (8 - c1_i)];
        }
      }
    }

    sf_mex_listen_for_ctrl_c(chartInstance->S);
  }

  c1_i34 = 0;
  c1_i = 0;
  for (c1_i35 = 0; c1_i35 < 8; c1_i35++) {
    for (c1_i36 = 0; c1_i36 < 8; c1_i36++) {
      c1_P[c1_i36 + c1_i34] = c1_M[(c1_i36 + c1_i) + 10];
    }

    c1_i34 += 8;
    c1_i += 9;
  }

  for (c1_i34 = 0; c1_i34 < 2; c1_i34++) {
    c1_d6 = 0.0;
    c1_i = 0;
    for (c1_i35 = 0; c1_i35 < 8; c1_i35++) {
      c1_d6 += c1_theta[c1_i + c1_i34] * c1_z[c1_i35];
      c1_i += 2;
    }

    c1_b_y[c1_i34] = c1_y[c1_i34] - c1_d6;
  }

  for (c1_i34 = 0; c1_i34 < 8; c1_i34++) {
    c1_b_z[c1_i34] = c1_M[c1_i34 + 1] / c1_M[0];
  }

  c1_i34 = 0;
  for (c1_i = 0; c1_i < 8; c1_i++) {
    for (c1_i35 = 0; c1_i35 < 2; c1_i35++) {
      c1_theta[c1_i35 + c1_i34] += c1_b_y[c1_i35] * c1_b_z[c1_i];
    }

    c1_i34 += 2;
  }
}

static void c1_b_inverseQR(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, real_T c1_z[10], real_T c1_y[7], real_T c1_theta[70], real_T
  c1_P[100], real_T c1_lambda)
{
  real_T c1_rho;
  real_T c1_c;
  int32_T c1_i37;
  int32_T c1_i;
  real_T c1_d7;
  int32_T c1_i38;
  real_T c1_b_z[10];
  real_T c1_M[121];
  int32_T c1_i39;
  real_T c1_b_y[7];
  if (c1_lambda < 0.0) {
    c1_eml_error(chartInstance);
  }

  c1_rho = muDoubleScalarSqrt(c1_lambda);
  if (c1_lambda < 0.0) {
    c1_eml_error(chartInstance);
  }

  c1_c = muDoubleScalarSqrt(c1_lambda);
  c1_i37 = 0;
  for (c1_i = 0; c1_i < 10; c1_i++) {
    c1_d7 = 0.0;
    for (c1_i38 = 0; c1_i38 < 10; c1_i38++) {
      c1_d7 += c1_z[c1_i38] * c1_P[c1_i38 + c1_i37];
    }

    c1_b_z[c1_i] = c1_d7 / c1_rho;
    c1_i37 += 10;
  }

  c1_M[0] = 1.0;
  c1_i37 = 0;
  for (c1_i = 0; c1_i < 10; c1_i++) {
    c1_M[c1_i37 + 11] = c1_b_z[c1_i];
    c1_i37 += 11;
  }

  for (c1_i37 = 0; c1_i37 < 10; c1_i37++) {
    c1_M[c1_i37 + 1] = 0.0;
  }

  c1_i37 = 0;
  c1_i = 0;
  for (c1_i38 = 0; c1_i38 < 10; c1_i38++) {
    for (c1_i39 = 0; c1_i39 < 10; c1_i39++) {
      c1_M[(c1_i39 + c1_i37) + 12] = c1_P[c1_i39 + c1_i] / c1_c;
    }

    c1_i37 += 11;
    c1_i += 10;
  }

  for (c1_i = 0; c1_i < 10; c1_i++) {
    if (!(c1_M[11 * (10 - c1_i)] == 0.0)) {
      if (c1_M[0] == 0.0) {
        for (c1_i37 = 0; c1_i37 < 11; c1_i37++) {
          c1_rho = c1_M[c1_i37 + 11 * (10 - c1_i)];
          c1_M[c1_i37 + 11 * (10 - c1_i)] = -c1_M[c1_i37];
          c1_M[c1_i37] = c1_rho;
        }
      } else {
        c1_rho = c1_M[11 * (10 - c1_i)] / c1_M[0];
        c1_c = muDoubleScalarSqrt(1.0 + muDoubleScalarPower(muDoubleScalarAbs
          (c1_rho), 2.0));
        for (c1_i37 = 0; c1_i37 < 11; c1_i37++) {
          c1_M[c1_i37] = (c1_M[c1_i37] + c1_rho * c1_M[c1_i37 + 11 * (10 - c1_i)])
            / c1_c;
          c1_M[c1_i37 + 11 * (10 - c1_i)] = -c1_rho * c1_M[c1_i37] + c1_c *
            c1_M[c1_i37 + 11 * (10 - c1_i)];
        }
      }
    }

    sf_mex_listen_for_ctrl_c(chartInstance->S);
  }

  c1_i37 = 0;
  c1_i = 0;
  for (c1_i38 = 0; c1_i38 < 10; c1_i38++) {
    for (c1_i39 = 0; c1_i39 < 10; c1_i39++) {
      c1_P[c1_i39 + c1_i37] = c1_M[(c1_i39 + c1_i) + 12];
    }

    c1_i37 += 10;
    c1_i += 11;
  }

  for (c1_i37 = 0; c1_i37 < 7; c1_i37++) {
    c1_d7 = 0.0;
    c1_i = 0;
    for (c1_i38 = 0; c1_i38 < 10; c1_i38++) {
      c1_d7 += c1_theta[c1_i + c1_i37] * c1_z[c1_i38];
      c1_i += 7;
    }

    c1_b_y[c1_i37] = c1_y[c1_i37] - c1_d7;
  }

  for (c1_i37 = 0; c1_i37 < 10; c1_i37++) {
    c1_b_z[c1_i37] = c1_M[c1_i37 + 1] / c1_M[0];
  }

  c1_i37 = 0;
  for (c1_i = 0; c1_i < 10; c1_i++) {
    for (c1_i38 = 0; c1_i38 < 7; c1_i38++) {
      c1_theta[c1_i38 + c1_i37] += c1_b_y[c1_i38] * c1_b_z[c1_i];
    }

    c1_i37 += 7;
  }
}

static void c1_b_sqrt(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
                      *chartInstance, creal_T *c1_x)
{
  real_T c1_absxi;
  real_T c1_a;
  real_T c1_absxr;
  if (c1_x->im == 0.0) {
    if (c1_x->re < 0.0) {
      c1_absxi = 0.0;
      c1_a = muDoubleScalarSqrt(muDoubleScalarAbs(c1_x->re));
    } else {
      c1_absxi = muDoubleScalarSqrt(c1_x->re);
      c1_a = 0.0;
    }
  } else if (c1_x->re == 0.0) {
    if (c1_x->im < 0.0) {
      c1_absxi = muDoubleScalarSqrt(-c1_x->im / 2.0);
      c1_a = -c1_absxi;
    } else {
      c1_absxi = muDoubleScalarSqrt(c1_x->im / 2.0);
      c1_a = c1_absxi;
    }
  } else if (muDoubleScalarIsNaN(c1_x->re) || muDoubleScalarIsNaN(c1_x->im)) {
    c1_absxi = rtNaN;
    c1_a = rtNaN;
  } else if (muDoubleScalarIsInf(c1_x->im)) {
    c1_absxi = rtInf;
    c1_a = c1_x->im;
  } else if (muDoubleScalarIsInf(c1_x->re)) {
    if (c1_x->re < 0.0) {
      c1_absxi = 0.0;
      c1_a = rtInf;
    } else {
      c1_absxi = rtInf;
      c1_a = 0.0;
    }
  } else {
    c1_absxr = muDoubleScalarAbs(c1_x->re);
    c1_absxi = muDoubleScalarAbs(c1_x->im);
    if ((c1_absxr > 4.4942328371557893E+307) || (c1_absxi >
         4.4942328371557893E+307)) {
      c1_absxr *= 0.5;
      c1_absxi *= 0.5;
      if (c1_absxr < c1_absxi) {
        c1_a = c1_absxr / c1_absxi;
        c1_absxi *= muDoubleScalarSqrt(c1_a * c1_a + 1.0);
      } else if (c1_absxr > c1_absxi) {
        c1_absxi /= c1_absxr;
        c1_absxi = muDoubleScalarSqrt(c1_absxi * c1_absxi + 1.0) * c1_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c1_absxi)) {
          c1_absxi = c1_absxr * 1.4142135623730951;
        }
      }

      if (c1_absxi > c1_absxr) {
        c1_absxi = muDoubleScalarSqrt(c1_absxi) * muDoubleScalarSqrt(1.0 +
          c1_absxr / c1_absxi);
      } else {
        c1_absxi = muDoubleScalarSqrt(c1_absxi) * 1.4142135623730951;
      }
    } else {
      if (c1_absxr < c1_absxi) {
        c1_a = c1_absxr / c1_absxi;
        c1_absxi *= muDoubleScalarSqrt(c1_a * c1_a + 1.0);
      } else if (c1_absxr > c1_absxi) {
        c1_absxi /= c1_absxr;
        c1_absxi = muDoubleScalarSqrt(c1_absxi * c1_absxi + 1.0) * c1_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c1_absxi)) {
          c1_absxi = c1_absxr * 1.4142135623730951;
        }
      }

      c1_absxi = muDoubleScalarSqrt((c1_absxi + c1_absxr) * 0.5);
    }

    if (c1_x->re > 0.0) {
      c1_a = 0.5 * (c1_x->im / c1_absxi);
    } else {
      if (c1_x->im < 0.0) {
        c1_a = -c1_absxi;
      } else {
        c1_a = c1_absxi;
      }

      c1_absxi = 0.5 * (c1_x->im / c1_a);
    }
  }

  c1_x->re = c1_absxi;
  c1_x->im = c1_a;
}

static void c1_exp(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                   creal_T *c1_x)
{
  real_T c1_r;
  creal_T c1_b_x;
  creal_T c1_c_x;
  c1_r = muDoubleScalarExp(c1_x->re / 2.0);
  c1_b_x = *c1_x;
  c1_c_x = *c1_x;
  c1_x->re = c1_r * (c1_r * muDoubleScalarCos(c1_b_x.im));
  c1_x->im = c1_r * (c1_r * muDoubleScalarSin(c1_c_x.im));
}

static void c1_log10(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance,
                     real_T *c1_x)
{
  if (*c1_x < 0.0) {
    c1_b_eml_error(chartInstance);
  }

  *c1_x = muDoubleScalarLog10(*c1_x);
}

static int32_T c1_div_s32_floor(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance, int32_T c1_numerator, int32_T c1_denominator)
{
  int32_T c1_quotient;
  uint32_T c1_absNumerator;
  uint32_T c1_absDenominator;
  int32_T c1_quotientNeedsNegation;
  uint32_T c1_tempAbsQuotient;
  if (c1_denominator == 0) {
    c1_quotient = c1_numerator >= 0 ? MAX_int32_T : MIN_int32_T;
    sf_mex_dividebyzero_error();
  } else {
    c1_absNumerator = (uint32_T)(c1_numerator >= 0 ? c1_numerator :
      -c1_numerator);
    c1_absDenominator = (uint32_T)(c1_denominator >= 0 ? c1_denominator :
      -c1_denominator);
    c1_quotientNeedsNegation = (c1_numerator < 0 != c1_denominator < 0);
    c1_tempAbsQuotient = c1_absNumerator / c1_absDenominator;
    if ((uint32_T)c1_quotientNeedsNegation) {
      c1_absNumerator %= c1_absDenominator;
      if (c1_absNumerator > (uint32_T)0) {
        c1_tempAbsQuotient++;
      }
    }

    c1_quotient = (uint32_T)c1_quotientNeedsNegation ? -(int32_T)
      c1_tempAbsQuotient : (int32_T)c1_tempAbsQuotient;
  }

  return c1_quotient;
}

static void init_dsm_address_info(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
void sf_c1_sim1_rlti_varx_fast_gust2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(605195618U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2127708425U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(873205752U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(971990910U);
}

mxArray *sf_c1_sim1_rlti_varx_fast_gust2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("aUoHy6eAEvmBi34ptddF4C");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(2);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(100);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

static const mxArray *sf_get_sim_state_info_c1_sim1_rlti_varx_fast_gust2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x10'type','srcId','name','auxInfo'{{M[1],M[21],T\"ERR\",},{M[1],M[5],T\"Fs1\",},{M[1],M[6],T\"Fs2\",},{M[1],M[9],T\"Wn\",},{M[1],M[4],T\"Ws\",},{M[1],M[10],T\"Zn\",},{M[4],M[0],T\"ABK\",S'l','i','p'{{M1x2[1236 1239],M[0],}}},{M[4],M[0],T\"CD\",S'l','i','p'{{M1x2[1240 1242],M[0],}}},{M[4],M[0],T\"Glk\",S'l','i','p'{{M1x2[1205 1208],M[0],}}},{M[4],M[0],T\"LK\",S'l','i','p'{{M1x2[1224 1226],M[0],}}}}",
    "100 S1x10'type','srcId','name','auxInfo'{{M[4],M[0],T\"Pabk\",S'l','i','p'{{M1x2[1231 1235],M[0],}}},{M[4],M[0],T\"Pcd\",S'l','i','p'{{M1x2[1227 1230],M[0],}}},{M[4],M[0],T\"Plk\",S'l','i','p'{{M1x2[1201 1204],M[0],}}},{M[4],M[0],T\"U1\",S'l','i','p'{{M1x2[1243 1245],M[0],}}},{M[4],M[0],T\"VARX\",S'l','i','p'{{M1x2[1213 1217],M[0],}}},{M[4],M[0],T\"VARX1\",S'l','i','p'{{M1x2[1218 1223],M[0],}}},{M[4],M[0],T\"X\",S'l','i','p'{{M1x2[1249 1250],M[0],}}},{M[4],M[0],T\"X1\",S'l','i','p'{{M1x2[1251 1253],M[0],}}},{M[4],M[0],T\"Y1\",S'l','i','p'{{M1x2[1246 1248],M[0],}}},{M[4],M[0],T\"Ylk\",S'l','i','p'{{M1x2[1209 1212],M[0],}}}}",
    "100 S1x7'type','srcId','name','auxInfo'{{M[4],M[0],T\"Z\",S'l','i','p'{{M1x2[1254 1255],M[0],}}},{M[4],M[0],T\"eta\",S'l','i','p'{{M1x2[1256 1259],M[0],}}},{M[4],M[0],T\"k\",S'l','i','p'{{M1x2[1262 1263],M[0],}}},{M[4],M[0],T\"saw\",S'l','i','p'{{M1x2[1270 1273],M[0],}}},{M[4],M[0],T\"start\",S'l','i','p'{{M1x2[1264 1269],M[0],}}},{M[4],M[0],T\"w\",S'l','i','p'{{M1x2[1260 1261],M[0],}}},{M[8],M[0],T\"is_active_c1_sim1_rlti_varx_fast_gust2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 27, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_sim1_rlti_varx_fast_gust2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void sf_opaque_initialize_c1_sim1_rlti_varx_fast_gust2(void
  *chartInstanceVar)
{
  initialize_params_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
  initialize_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_sim1_rlti_varx_fast_gust2(void *chartInstanceVar)
{
  enable_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_sim1_rlti_varx_fast_gust2(void
  *chartInstanceVar)
{
  disable_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_sim1_rlti_varx_fast_gust2(void
  *chartInstanceVar)
{
  sf_c1_sim1_rlti_varx_fast_gust2((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_sim1_rlti_varx_fast_gust2();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c1_sim1_rlti_varx_fast_gust2(SimStruct* S,
  const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_sim1_rlti_varx_fast_gust2();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*)chartInfo->chartInstance,
     mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_sim1_rlti_varx_fast_gust2
  (SimStruct* S)
{
  return sf_internal_get_sim_state_c1_sim1_rlti_varx_fast_gust2(S);
}

static void sf_opaque_set_sim_state_c1_sim1_rlti_varx_fast_gust2(SimStruct* S,
  const mxArray *st)
{
  sf_internal_set_sim_state_c1_sim1_rlti_varx_fast_gust2(S, st);
}

static void sf_opaque_terminate_c1_sim1_rlti_varx_fast_gust2(void
  *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*)
                    chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
    }

    finalize_c1_sim1_rlti_varx_fast_gust2
      ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
    free((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }

  unload_sim1_rlti_varx_fast_gust2_optimization_info();
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_sim1_rlti_varx_fast_gust2
    ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_sim1_rlti_varx_fast_gust2(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c1_sim1_rlti_varx_fast_gust2
      ((SFc1_sim1_rlti_varx_fast_gust2InstanceStruct*)(((ChartInfoStruct *)
         ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_sim1_rlti_varx_fast_gust2(SimStruct *S)
{
  /* Actual parameters from chart:
     Du Dy W
   */
  const char_T *rtParamNames[] = { "p1", "p2", "p3" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for Du*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for Dy*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for W*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_sim1_rlti_varx_fast_gust2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,infoStruct,1,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,infoStruct,1,4);
      sf_mark_chart_reusable_outputs(S,infoStruct,1,6);
    }

    sf_set_rtw_dwork_info(S,infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3797036093U));
  ssSetChecksum1(S,(4190760057U));
  ssSetChecksum2(S,(1920555931U));
  ssSetChecksum3(S,(2587818289U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
}

static void mdlRTW_c1_sim1_rlti_varx_fast_gust2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_sim1_rlti_varx_fast_gust2(SimStruct *S)
{
  SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *chartInstance;
  chartInstance = (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct *)malloc(sizeof
    (SFc1_sim1_rlti_varx_fast_gust2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_sim1_rlti_varx_fast_gust2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c1_sim1_rlti_varx_fast_gust2;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
}

void c1_sim1_rlti_varx_fast_gust2_method_dispatcher(SimStruct *S, int_T method,
  void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_sim1_rlti_varx_fast_gust2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_sim1_rlti_varx_fast_gust2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_sim1_rlti_varx_fast_gust2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_sim1_rlti_varx_fast_gust2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
