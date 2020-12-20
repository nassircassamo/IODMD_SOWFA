/* Include files */

#include "blascompat32.h"
#include "sim1_rlti_varx_gust2_sfun.h"
#include "c2_sim1_rlti_varx_gust2.h"
#include "mwmathutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void initialize_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void initialize_params_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void enable_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void disable_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void set_sim_state_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance, const mxArray *c2_st);
static void finalize_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void sf_c2_sim1_rlti_varx_gust2(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static void c2_chartstep_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void initSimStructsc2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber);
static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[106]);
static void c2_b_info_helper(c2_ResolvedFunctionInfo c2_info[106]);
static real_T c2_mrdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_A, real_T c2_B);
static real_T c2_rdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_x, real_T c2_y);
static void c2_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   real_T c2_I[23104]);
static void c2_b_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[100]);
static void c2_c_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[64]);
static void c2_logspace(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_d1, real_T c2_d2, real_T c2_y[300]);
static void c2_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static void c2_b_rdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_x, real_T c2_y[2], real_T c2_z[2]);
static void c2_diag(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    real_T c2_v[2], real_T c2_d[4]);
static void c2_eml_error(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance);
static void c2_damp(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    real_T c2_a[49], real_T c2_h, real_T c2_wn[7], real_T c2_z[7]);
static void c2_eig(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   real_T c2_A[49], creal_T c2_V[7]);
static void c2_eml_matlab_zlartg(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_f, creal_T c2_g, real_T *c2_cs, creal_T *c2_sn,
  creal_T *c2_r);
static void c2_eml_matlab_zhgeqz(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_A[49], int32_T c2_ilo, int32_T c2_ihi, real_T
  *c2_info, creal_T c2_alpha1[7], creal_T c2_beta1[7]);
static real_T c2_eml_matlab_zlanhs(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_A[49], int32_T c2_ilo, int32_T c2_ihi);
static creal_T c2_eml_div(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  creal_T c2_x, real_T c2_y);
static void c2_b_eml_matlab_zlartg(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_f, creal_T c2_g, real_T *c2_cs, creal_T *c2_sn);
static void c2_b_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static void c2_c_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static void c2_d_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[49]);
static void c2_mldivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  creal_T c2_A[49], real_T c2_B[7], creal_T c2_Y[7]);
static creal_T c2_b_eml_div(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_x, creal_T c2_y);
static void c2_d_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static void c2_abs(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   creal_T c2_x[2], real_T c2_y[2]);
static void c2_b_eml_error(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);
static real_T c2_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_Fs1, const char_T *c2_identifier);
static real_T c2_b_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_c_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_Wn, const char_T *c2_identifier, real_T
  c2_y[7]);
static void c2_d_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7]);
static void c2_e_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_ABK, const char_T *c2_identifier, real_T
  c2_y[70]);
static void c2_f_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[70]);
static void c2_g_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_CD, const char_T *c2_identifier, real_T
  c2_y[16]);
static void c2_h_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[16]);
static void c2_i_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_LK, const char_T *c2_identifier, real_T
  c2_y[15000]);
static void c2_j_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[15000]);
static void c2_k_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Pabk, const char_T *c2_identifier, real_T
  c2_y[100]);
static void c2_l_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[100]);
static void c2_m_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Pcd, const char_T *c2_identifier, real_T
  c2_y[64]);
static void c2_n_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[64]);
static void c2_o_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Plk, const char_T *c2_identifier, real_T
  c2_y[23104]);
static void c2_p_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[23104]);
static real_T c2_q_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_U1, const char_T *c2_identifier);
static real_T c2_r_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_s_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_VARX, const char_T *c2_identifier, real_T
  c2_y[304]);
static void c2_t_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[304]);
static void c2_u_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_X, const char_T *c2_identifier, real_T
  c2_y[7]);
static void c2_v_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7]);
static void c2_w_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_X1, const char_T *c2_identifier, real_T
  c2_y[7]);
static void c2_x_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7]);
static void c2_y_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Y1, const char_T *c2_identifier, real_T
  c2_y[2]);
static void c2_ab_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[2]);
static void c2_bb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Zp, const char_T *c2_identifier, real_T
  c2_y[150]);
static void c2_cb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[150]);
static real_T c2_db_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_k, const char_T *c2_identifier);
static real_T c2_eb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static real_T c2_fb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_saw, const char_T *c2_identifier);
static real_T c2_gb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static real_T c2_hb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_start, const char_T *c2_identifier);
static real_T c2_ib_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_jb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_w, const char_T *c2_identifier, real_T
  c2_y[300]);
static void c2_kb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[300]);
static uint8_T c2_lb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_sim1_rlti_varx_gust2, const
  char_T *c2_identifier);
static uint8_T c2_mb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_rls_ew_reg(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_z[152], real_T c2_y[2], real_T c2_theta[304], real_T c2_P[23104],
  real_T c2_lambda, real_T c2_reg, real_T *c2_b_saw);
static void c2_rls_ew(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                      real_T c2_z[8], real_T c2_y[2], real_T c2_theta[16],
                      real_T c2_P[64], real_T c2_lambda);
static void c2_b_rls_ew(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_z[10], real_T c2_y[7], real_T c2_theta[70], real_T c2_P[100], real_T
  c2_lambda);
static void c2_sqrt(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    creal_T *c2_x);
static void c2_exp(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   creal_T *c2_x);
static void c2_log10(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T *c2_x);
static int32_T c2_div_s32_floor(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, int32_T c2_numerator, int32_T c2_denominator);
static void init_dsm_address_info(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c2_Plk_not_empty = FALSE;
  chartInstance->c2_VARX_not_empty = FALSE;
  chartInstance->c2_LK_not_empty = FALSE;
  chartInstance->c2_Pcd_not_empty = FALSE;
  chartInstance->c2_Pabk_not_empty = FALSE;
  chartInstance->c2_ABK_not_empty = FALSE;
  chartInstance->c2_CD_not_empty = FALSE;
  chartInstance->c2_U1_not_empty = FALSE;
  chartInstance->c2_Y1_not_empty = FALSE;
  chartInstance->c2_X_not_empty = FALSE;
  chartInstance->c2_X1_not_empty = FALSE;
  chartInstance->c2_Zp_not_empty = FALSE;
  chartInstance->c2_w_not_empty = FALSE;
  chartInstance->c2_k_not_empty = FALSE;
  chartInstance->c2_start_not_empty = FALSE;
  chartInstance->c2_saw_not_empty = FALSE;
  chartInstance->c2_is_active_c2_sim1_rlti_varx_gust2 = 0U;
}

static void initialize_params_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  real_T c2_dv0[700];
  int32_T c2_i0;
  real_T c2_d0;
  real_T c2_dv1[4];
  sf_set_error_prefix_string(
    "Error evaluating data 'W' in the parent workspace.\n");
  sf_mex_import_named("W", sf_mex_get_sfun_param(chartInstance->S, 2, 0), c2_dv0,
                      0, 0, 0U, 1, 0U, 2, 7, 100);
  for (c2_i0 = 0; c2_i0 < 700; c2_i0++) {
    chartInstance->c2_W[c2_i0] = c2_dv0[c2_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Du' in the parent workspace.\n");
  sf_mex_import_named("Du", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c2_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c2_Du = c2_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Dy' in the parent workspace.\n");
  sf_mex_import_named("Dy", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      c2_dv1, 0, 0, 0U, 1, 0U, 2, 2, 2);
  for (c2_i0 = 0; c2_i0 < 4; c2_i0++) {
    chartInstance->c2_Dy[c2_i0] = c2_dv1[c2_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static const mxArray *get_sim_state_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  int32_T c2_i1;
  real_T c2_c_u[7];
  const mxArray *c2_d_y = NULL;
  real_T c2_d_u;
  const mxArray *c2_e_y = NULL;
  real_T c2_e_u[7];
  const mxArray *c2_f_y = NULL;
  real_T c2_f_u[70];
  const mxArray *c2_g_y = NULL;
  real_T c2_g_u[16];
  const mxArray *c2_h_y = NULL;
  real_T c2_h_u[15000];
  const mxArray *c2_i_y = NULL;
  real_T c2_i_u[100];
  const mxArray *c2_j_y = NULL;
  real_T c2_j_u[64];
  const mxArray *c2_k_y = NULL;
  static real_T c2_k_u[23104];
  const mxArray *c2_l_y = NULL;
  real_T c2_l_u;
  const mxArray *c2_m_y = NULL;
  real_T c2_m_u[304];
  const mxArray *c2_n_y = NULL;
  real_T c2_n_u[7];
  const mxArray *c2_o_y = NULL;
  real_T c2_o_u[7];
  const mxArray *c2_p_y = NULL;
  real_T c2_p_u[2];
  const mxArray *c2_q_y = NULL;
  real_T c2_q_u[150];
  const mxArray *c2_r_y = NULL;
  real_T c2_r_u;
  const mxArray *c2_s_y = NULL;
  real_T c2_s_u;
  const mxArray *c2_t_y = NULL;
  real_T c2_t_u;
  const mxArray *c2_u_y = NULL;
  real_T c2_u_u[300];
  const mxArray *c2_v_y = NULL;
  uint8_T c2_v_u;
  const mxArray *c2_w_y = NULL;
  real_T *c2_Fs1;
  real_T *c2_Fs2;
  real_T *c2_Ws;
  real_T (*c2_Zn)[7];
  real_T (*c2_Wn)[7];
  c2_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellarray(22), FALSE);
  c2_u = *c2_Fs1;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_b_u = *c2_Fs2;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  for (c2_i1 = 0; c2_i1 < 7; c2_i1++) {
    c2_c_u[c2_i1] = (*c2_Wn)[c2_i1];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  c2_d_u = *c2_Ws;
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 3, c2_e_y);
  for (c2_i1 = 0; c2_i1 < 7; c2_i1++) {
    c2_e_u[c2_i1] = (*c2_Zn)[c2_i1];
  }

  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", c2_e_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c2_y, 4, c2_f_y);
  for (c2_i1 = 0; c2_i1 < 70; c2_i1++) {
    c2_f_u[c2_i1] = chartInstance->c2_ABK[c2_i1];
  }

  c2_g_y = NULL;
  if (!chartInstance->c2_ABK_not_empty) {
    sf_mex_assign(&c2_g_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_g_y, sf_mex_create("y", c2_f_u, 0, 0U, 1U, 0U, 2, 7, 10),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 5, c2_g_y);
  for (c2_i1 = 0; c2_i1 < 16; c2_i1++) {
    c2_g_u[c2_i1] = chartInstance->c2_CD[c2_i1];
  }

  c2_h_y = NULL;
  if (!chartInstance->c2_CD_not_empty) {
    sf_mex_assign(&c2_h_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_h_y, sf_mex_create("y", c2_g_u, 0, 0U, 1U, 0U, 2, 2, 8),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 6, c2_h_y);
  for (c2_i1 = 0; c2_i1 < 15000; c2_i1++) {
    c2_h_u[c2_i1] = chartInstance->c2_LK[c2_i1];
  }

  c2_i_y = NULL;
  if (!chartInstance->c2_LK_not_empty) {
    sf_mex_assign(&c2_i_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_i_y, sf_mex_create("y", c2_h_u, 0, 0U, 1U, 0U, 2, 100, 150),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 7, c2_i_y);
  for (c2_i1 = 0; c2_i1 < 100; c2_i1++) {
    c2_i_u[c2_i1] = chartInstance->c2_Pabk[c2_i1];
  }

  c2_j_y = NULL;
  if (!chartInstance->c2_Pabk_not_empty) {
    sf_mex_assign(&c2_j_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_j_y, sf_mex_create("y", c2_i_u, 0, 0U, 1U, 0U, 2, 10, 10),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 8, c2_j_y);
  for (c2_i1 = 0; c2_i1 < 64; c2_i1++) {
    c2_j_u[c2_i1] = chartInstance->c2_Pcd[c2_i1];
  }

  c2_k_y = NULL;
  if (!chartInstance->c2_Pcd_not_empty) {
    sf_mex_assign(&c2_k_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_k_y, sf_mex_create("y", c2_j_u, 0, 0U, 1U, 0U, 2, 8, 8),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 9, c2_k_y);
  for (c2_i1 = 0; c2_i1 < 23104; c2_i1++) {
    c2_k_u[c2_i1] = chartInstance->c2_Plk[c2_i1];
  }

  c2_l_y = NULL;
  if (!chartInstance->c2_Plk_not_empty) {
    sf_mex_assign(&c2_l_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_l_y, sf_mex_create("y", c2_k_u, 0, 0U, 1U, 0U, 2, 152, 152),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 10, c2_l_y);
  c2_l_u = chartInstance->c2_U1;
  c2_m_y = NULL;
  if (!chartInstance->c2_U1_not_empty) {
    sf_mex_assign(&c2_m_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_m_y, sf_mex_create("y", &c2_l_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c2_y, 11, c2_m_y);
  for (c2_i1 = 0; c2_i1 < 304; c2_i1++) {
    c2_m_u[c2_i1] = chartInstance->c2_VARX[c2_i1];
  }

  c2_n_y = NULL;
  if (!chartInstance->c2_VARX_not_empty) {
    sf_mex_assign(&c2_n_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_n_y, sf_mex_create("y", c2_m_u, 0, 0U, 1U, 0U, 2, 2, 152),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 12, c2_n_y);
  for (c2_i1 = 0; c2_i1 < 7; c2_i1++) {
    c2_n_u[c2_i1] = chartInstance->c2_X[c2_i1];
  }

  c2_o_y = NULL;
  if (!chartInstance->c2_X_not_empty) {
    sf_mex_assign(&c2_o_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_o_y, sf_mex_create("y", c2_n_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 13, c2_o_y);
  for (c2_i1 = 0; c2_i1 < 7; c2_i1++) {
    c2_o_u[c2_i1] = chartInstance->c2_X1[c2_i1];
  }

  c2_p_y = NULL;
  if (!chartInstance->c2_X1_not_empty) {
    sf_mex_assign(&c2_p_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_p_y, sf_mex_create("y", c2_o_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 14, c2_p_y);
  for (c2_i1 = 0; c2_i1 < 2; c2_i1++) {
    c2_p_u[c2_i1] = chartInstance->c2_Y1[c2_i1];
  }

  c2_q_y = NULL;
  if (!chartInstance->c2_Y1_not_empty) {
    sf_mex_assign(&c2_q_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_q_y, sf_mex_create("y", c2_p_u, 0, 0U, 1U, 0U, 1, 2),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 15, c2_q_y);
  for (c2_i1 = 0; c2_i1 < 150; c2_i1++) {
    c2_q_u[c2_i1] = chartInstance->c2_Zp[c2_i1];
  }

  c2_r_y = NULL;
  if (!chartInstance->c2_Zp_not_empty) {
    sf_mex_assign(&c2_r_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_r_y, sf_mex_create("y", c2_q_u, 0, 0U, 1U, 0U, 1, 150),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 16, c2_r_y);
  c2_r_u = chartInstance->c2_k;
  c2_s_y = NULL;
  if (!chartInstance->c2_k_not_empty) {
    sf_mex_assign(&c2_s_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_s_y, sf_mex_create("y", &c2_r_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c2_y, 17, c2_s_y);
  c2_s_u = chartInstance->c2_saw;
  c2_t_y = NULL;
  if (!chartInstance->c2_saw_not_empty) {
    sf_mex_assign(&c2_t_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_t_y, sf_mex_create("y", &c2_s_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c2_y, 18, c2_t_y);
  c2_t_u = chartInstance->c2_start;
  c2_u_y = NULL;
  if (!chartInstance->c2_start_not_empty) {
    sf_mex_assign(&c2_u_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_u_y, sf_mex_create("y", &c2_t_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c2_y, 19, c2_u_y);
  for (c2_i1 = 0; c2_i1 < 300; c2_i1++) {
    c2_u_u[c2_i1] = chartInstance->c2_w[c2_i1];
  }

  c2_v_y = NULL;
  if (!chartInstance->c2_w_not_empty) {
    sf_mex_assign(&c2_v_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c2_v_y, sf_mex_create("y", c2_u_u, 0, 0U, 1U, 0U, 2, 1, 300),
                  FALSE);
  }

  sf_mex_setcell(c2_y, 20, c2_v_y);
  c2_v_u = chartInstance->c2_is_active_c2_sim1_rlti_varx_gust2;
  c2_w_y = NULL;
  sf_mex_assign(&c2_w_y, sf_mex_create("y", &c2_v_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 21, c2_w_y);
  sf_mex_assign(&c2_st, c2_y, FALSE);
  return c2_st;
}

static void set_sim_state_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv2[7];
  int32_T c2_i2;
  real_T c2_dv3[70];
  real_T c2_dv4[16];
  static real_T c2_dv5[15000];
  real_T c2_dv6[100];
  real_T c2_dv7[64];
  static real_T c2_dv8[23104];
  real_T c2_dv9[304];
  real_T c2_dv10[2];
  real_T c2_dv11[150];
  real_T c2_dv12[300];
  real_T *c2_Fs1;
  real_T *c2_Fs2;
  real_T *c2_Ws;
  real_T (*c2_Wn)[7];
  real_T (*c2_Zn)[7];
  c2_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_u = sf_mex_dup(c2_st);
  *c2_Fs1 = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 0)),
    "Fs1");
  *c2_Fs2 = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 1)),
    "Fs2");
  c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 2)), "Wn",
                        c2_dv2);
  for (c2_i2 = 0; c2_i2 < 7; c2_i2++) {
    (*c2_Wn)[c2_i2] = c2_dv2[c2_i2];
  }

  *c2_Ws = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 3)),
    "Ws");
  c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 4)), "Zn",
                        c2_dv2);
  for (c2_i2 = 0; c2_i2 < 7; c2_i2++) {
    (*c2_Zn)[c2_i2] = c2_dv2[c2_i2];
  }

  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 5)),
                        "ABK", c2_dv3);
  for (c2_i2 = 0; c2_i2 < 70; c2_i2++) {
    chartInstance->c2_ABK[c2_i2] = c2_dv3[c2_i2];
  }

  c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 6)), "CD",
                        c2_dv4);
  for (c2_i2 = 0; c2_i2 < 16; c2_i2++) {
    chartInstance->c2_CD[c2_i2] = c2_dv4[c2_i2];
  }

  c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 7)), "LK",
                        c2_dv5);
  for (c2_i2 = 0; c2_i2 < 15000; c2_i2++) {
    chartInstance->c2_LK[c2_i2] = c2_dv5[c2_i2];
  }

  c2_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 8)),
                        "Pabk", c2_dv6);
  for (c2_i2 = 0; c2_i2 < 100; c2_i2++) {
    chartInstance->c2_Pabk[c2_i2] = c2_dv6[c2_i2];
  }

  c2_m_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 9)),
                        "Pcd", c2_dv7);
  for (c2_i2 = 0; c2_i2 < 64; c2_i2++) {
    chartInstance->c2_Pcd[c2_i2] = c2_dv7[c2_i2];
  }

  c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 10)),
                        "Plk", c2_dv8);
  for (c2_i2 = 0; c2_i2 < 23104; c2_i2++) {
    chartInstance->c2_Plk[c2_i2] = c2_dv8[c2_i2];
  }

  chartInstance->c2_U1 = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 11)), "U1");
  c2_s_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 12)),
                        "VARX", c2_dv9);
  for (c2_i2 = 0; c2_i2 < 304; c2_i2++) {
    chartInstance->c2_VARX[c2_i2] = c2_dv9[c2_i2];
  }

  c2_u_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 13)), "X",
                        c2_dv2);
  for (c2_i2 = 0; c2_i2 < 7; c2_i2++) {
    chartInstance->c2_X[c2_i2] = c2_dv2[c2_i2];
  }

  c2_w_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 14)),
                        "X1", c2_dv2);
  for (c2_i2 = 0; c2_i2 < 7; c2_i2++) {
    chartInstance->c2_X1[c2_i2] = c2_dv2[c2_i2];
  }

  c2_y_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 15)),
                        "Y1", c2_dv10);
  for (c2_i2 = 0; c2_i2 < 2; c2_i2++) {
    chartInstance->c2_Y1[c2_i2] = c2_dv10[c2_i2];
  }

  c2_bb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 16)),
    "Zp", c2_dv11);
  for (c2_i2 = 0; c2_i2 < 150; c2_i2++) {
    chartInstance->c2_Zp[c2_i2] = c2_dv11[c2_i2];
  }

  chartInstance->c2_k = c2_db_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 17)), "k");
  chartInstance->c2_saw = c2_fb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 18)), "saw");
  chartInstance->c2_start = c2_hb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 19)), "start");
  c2_jb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 20)),
    "w", c2_dv12);
  for (c2_i2 = 0; c2_i2 < 300; c2_i2++) {
    chartInstance->c2_w[c2_i2] = c2_dv12[c2_i2];
  }

  chartInstance->c2_is_active_c2_sim1_rlti_varx_gust2 = c2_lb_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 21)),
     "is_active_c2_sim1_rlti_varx_gust2");
  sf_mex_destroy(&c2_u);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
}

static void sf_c2_sim1_rlti_varx_gust2(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  c2_chartstep_c2_sim1_rlti_varx_gust2(chartInstance);
}

static void c2_chartstep_c2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  real_T c2_ON;
  real_T c2_U;
  int32_T c2_i3;
  real_T c2_Fs[2];
  real_T c2_b_W[700];
  real_T c2_b_Du;
  real_T c2_b_Dy[4];
  real_T c2_alpha1;
  int32_T c2_ldc;
  int32_T c2_i;
  real_T c2_dv13[300];
  real_T c2_c_Dy[2];
  real_T c2_dv14[2];
  real_T c2_a[4];
  real_T c2_dv15[147];
  real_T c2_b_VARX[304];
  static real_T c2_b_Plk[23104];
  real_T c2_dv16[152];
  int32_T c2_n;
  int32_T c2_b_k;
  int32_T c2_lda;
  int32_T c2_ldb;
  char_T c2_TRANSA;
  char_T c2_TRANSB;
  real_T c2_y[1050];
  real_T c2_dv17[8];
  real_T c2_dv18[10];
  real_T c2_A[49];
  real_T c2_b_a[7];
  real_T c2_C[14];
  real_T c2_dv19[49];
  real_T c2_Zn[7];
  real_T c2_Wn[7];
  creal_T c2_b_y;
  creal_T c2_c_y[49];
  real_T c2_c_a[7];
  creal_T c2_b[7];
  creal_T c2_b_C[14];
  creal_T c2_c_C[2];
  real_T *c2_Ws;
  real_T *c2_Fs1;
  real_T *c2_Fs2;
  real_T *c2_b_ON;
  real_T *c2_b_U;
  real_T (*c2_b_Wn)[7];
  real_T (*c2_b_Zn)[7];
  real_T (*c2_Y)[2];
  c2_b_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_b_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Y = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c2_b_U = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c2_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_b_ON = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c2_ON = *c2_b_ON;
  c2_U = *c2_b_U;
  for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
    c2_Fs[c2_i3] = (*c2_Y)[c2_i3];
  }

  for (c2_i3 = 0; c2_i3 < 700; c2_i3++) {
    c2_b_W[c2_i3] = chartInstance->c2_W[c2_i3];
  }

  c2_b_Du = chartInstance->c2_Du;
  for (c2_i3 = 0; c2_i3 < 4; c2_i3++) {
    c2_b_Dy[c2_i3] = chartInstance->c2_Dy[c2_i3];
  }

  if (!chartInstance->c2_Plk_not_empty) {
    c2_eye(chartInstance, chartInstance->c2_Plk);
    c2_alpha1 = c2_mrdivide(chartInstance, 1.0, 100.0);
    for (c2_i3 = 0; c2_i3 < 23104; c2_i3++) {
      chartInstance->c2_Plk[c2_i3] *= c2_alpha1;
    }

    chartInstance->c2_Plk_not_empty = TRUE;
    for (c2_i3 = 0; c2_i3 < 304; c2_i3++) {
      chartInstance->c2_VARX[c2_i3] = 0.0;
    }

    chartInstance->c2_VARX_not_empty = TRUE;
    chartInstance->c2_saw = 0.0;
    chartInstance->c2_saw_not_empty = TRUE;
    for (c2_i3 = 0; c2_i3 < 15000; c2_i3++) {
      chartInstance->c2_LK[c2_i3] = 0.0;
    }

    chartInstance->c2_LK_not_empty = TRUE;
    c2_b_eye(chartInstance, chartInstance->c2_Pabk);
    c2_alpha1 = c2_mrdivide(chartInstance, 1.0, 1.0E-6);
    for (c2_i3 = 0; c2_i3 < 100; c2_i3++) {
      chartInstance->c2_Pabk[c2_i3] *= c2_alpha1;
    }

    chartInstance->c2_Pabk_not_empty = TRUE;
    c2_i3 = 0;
    for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
      for (c2_i = 0; c2_i < 7; c2_i++) {
        chartInstance->c2_ABK[c2_i + c2_i3] = 0.0;
      }

      c2_i3 += 7;
    }

    for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
      chartInstance->c2_ABK[c2_i3 + 49] = 0.0;
    }

    c2_i3 = 0;
    for (c2_ldc = 0; c2_ldc < 2; c2_ldc++) {
      for (c2_i = 0; c2_i < 7; c2_i++) {
        chartInstance->c2_ABK[(c2_i + c2_i3) + 56] = 0.0;
      }

      c2_i3 += 7;
    }

    chartInstance->c2_ABK_not_empty = TRUE;
    c2_c_eye(chartInstance, chartInstance->c2_Pcd);
    c2_alpha1 = c2_mrdivide(chartInstance, 1.0, 1.0E-6);
    for (c2_i3 = 0; c2_i3 < 64; c2_i3++) {
      chartInstance->c2_Pcd[c2_i3] *= c2_alpha1;
    }

    chartInstance->c2_Pcd_not_empty = TRUE;
    c2_i3 = 0;
    for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
      for (c2_i = 0; c2_i < 2; c2_i++) {
        chartInstance->c2_CD[c2_i + c2_i3] = 0.0;
      }

      c2_i3 += 2;
    }

    chartInstance->c2_CD_not_empty = TRUE;
    for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
      chartInstance->c2_CD[c2_i3 + 14] = 0.0;
      chartInstance->c2_Y1[c2_i3] = 0.0;
    }

    chartInstance->c2_Y1_not_empty = TRUE;
    chartInstance->c2_U1 = 0.0;
    chartInstance->c2_U1_not_empty = TRUE;
    chartInstance->c2_X1_not_empty = TRUE;
    for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
      chartInstance->c2_X1[c2_i3] = 0.0;
      chartInstance->c2_X[c2_i3] = 0.0;
    }

    chartInstance->c2_X_not_empty = TRUE;
    for (c2_i3 = 0; c2_i3 < 150; c2_i3++) {
      chartInstance->c2_Zp[c2_i3] = 0.0;
    }

    chartInstance->c2_Zp_not_empty = TRUE;
    c2_logspace(chartInstance, -1.0, 1.0, c2_dv13);
    for (c2_i3 = 0; c2_i3 < 300; c2_i3++) {
      chartInstance->c2_w[c2_i3] = c2_dv13[c2_i3];
    }

    chartInstance->c2_w_not_empty = TRUE;
    chartInstance->c2_k = 1.0;
    chartInstance->c2_k_not_empty = TRUE;
    chartInstance->c2_start = 1.0;
    chartInstance->c2_start_not_empty = TRUE;
  }

  if (c2_ON > 0.5) {
    c2_U *= c2_rdivide(chartInstance, 1.0, c2_b_Du);
    c2_i3 = 0;
    for (c2_ldc = 0; c2_ldc < 2; c2_ldc++) {
      c2_c_Dy[c2_ldc] = c2_b_Dy[c2_i3];
      c2_i3 += 3;
    }

    c2_b_rdivide(chartInstance, 1.0, c2_c_Dy, c2_dv14);
    c2_diag(chartInstance, c2_dv14, c2_a);
    for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
      c2_c_Dy[c2_i3] = 0.0;
      c2_ldc = 0;
      for (c2_i = 0; c2_i < 2; c2_i++) {
        c2_c_Dy[c2_i3] += c2_a[c2_ldc + c2_i3] * c2_Fs[c2_i];
        c2_ldc += 2;
      }
    }

    for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
      c2_Fs[c2_i3] = c2_c_Dy[c2_i3];
    }

    for (c2_i3 = 0; c2_i3 < 147; c2_i3++) {
      c2_dv15[c2_i3] = chartInstance->c2_Zp[c2_i3 + 3];
    }

    for (c2_i3 = 0; c2_i3 < 147; c2_i3++) {
      chartInstance->c2_Zp[c2_i3] = c2_dv15[c2_i3];
    }

    chartInstance->c2_Zp[147] = chartInstance->c2_U1;
    for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
      chartInstance->c2_Zp[c2_i3 + 148] = chartInstance->c2_Y1[c2_i3];
    }

    if ((chartInstance->c2_start < 0.5) || (chartInstance->c2_k >= 50.0)) {
      for (c2_i3 = 0; c2_i3 < 304; c2_i3++) {
        c2_b_VARX[c2_i3] = chartInstance->c2_VARX[c2_i3];
      }

      for (c2_i3 = 0; c2_i3 < 23104; c2_i3++) {
        c2_b_Plk[c2_i3] = chartInstance->c2_Plk[c2_i3];
      }

      c2_ON = chartInstance->c2_saw;
      for (c2_i3 = 0; c2_i3 < 150; c2_i3++) {
        c2_dv16[c2_i3] = chartInstance->c2_Zp[c2_i3];
      }

      c2_dv16[150] = c2_U;
      c2_dv16[151] = 1.0;
      c2_rls_ew_reg(chartInstance, c2_dv16, c2_Fs, c2_b_VARX, c2_b_Plk, 0.999,
                    1.0, &c2_ON);
      for (c2_i3 = 0; c2_i3 < 304; c2_i3++) {
        chartInstance->c2_VARX[c2_i3] = c2_b_VARX[c2_i3];
      }

      for (c2_i3 = 0; c2_i3 < 23104; c2_i3++) {
        chartInstance->c2_Plk[c2_i3] = c2_b_Plk[c2_i3];
      }

      chartInstance->c2_saw = c2_ON;
    }

    if ((chartInstance->c2_start < 0.5) || (chartInstance->c2_k >= 150.0)) {
      for (c2_i3 = 0; c2_i3 < 15000; c2_i3++) {
        chartInstance->c2_LK[c2_i3] = 0.0;
      }

      for (c2_i = 0; c2_i < 50; c2_i++) {
        c2_n = 0;
        while (c2_n <= 49 - c2_i) {
          c2_b_k = (c2_i << 1) - 1;
          c2_lda = (c2_i + c2_n) * 3;
          c2_ldb = c2_n * 3 - 1;
          for (c2_i3 = 0; c2_i3 < 3; c2_i3++) {
            for (c2_ldc = 0; c2_ldc < 2; c2_ldc++) {
              chartInstance->c2_LK[((c2_ldc + c2_b_k) + 100 *
                                    (sf_mex_lw_bounds_check((c2_i3 + c2_lda) + 1,
                1, 150) - 1)) + 1] = chartInstance->c2_VARX[c2_ldc + (((c2_i3 +
                c2_ldb) + 1) << 1)];
            }
          }

          c2_n++;
          sf_mex_listen_for_ctrl_c(chartInstance->S);
        }

        sf_mex_listen_for_ctrl_c(chartInstance->S);
      }

      c2_i = 7;
      c2_n = 150;
      c2_b_k = 100;
      c2_alpha1 = 1.0;
      c2_lda = 7;
      c2_ldb = 100;
      c2_ON = 0.0;
      c2_ldc = 7;
      c2_TRANSA = 'N';
      c2_TRANSB = 'N';
      for (c2_i3 = 0; c2_i3 < 1050; c2_i3++) {
        c2_y[c2_i3] = 0.0;
      }

      dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_i, &c2_n, &c2_b_k, &c2_alpha1,
              &c2_b_W[0], &c2_lda, &chartInstance->c2_LK[0], &c2_ldb, &c2_ON,
              &c2_y[0], &c2_ldc);
      c2_i = 7;
      c2_n = 1;
      c2_b_k = 150;
      c2_alpha1 = 1.0;
      c2_lda = 7;
      c2_ldb = 150;
      c2_ON = 0.0;
      c2_ldc = 7;
      c2_TRANSA = 'N';
      c2_TRANSB = 'N';
      for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
        chartInstance->c2_X[c2_i3] = 0.0;
      }

      dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_i, &c2_n, &c2_b_k, &c2_alpha1, &c2_y[0],
              &c2_lda, &chartInstance->c2_Zp[0], &c2_ldb, &c2_ON,
              &chartInstance->c2_X[0], &c2_ldc);
    }

    if ((chartInstance->c2_start < 0.5) || (chartInstance->c2_k >= 151.0)) {
      for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
        c2_dv17[c2_i3] = chartInstance->c2_X1[c2_i3];
      }

      c2_dv17[7] = chartInstance->c2_U1;
      c2_rls_ew(chartInstance, c2_dv17, chartInstance->c2_Y1,
                chartInstance->c2_CD, chartInstance->c2_Pcd, 0.999);
      for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
        c2_dv17[c2_i3] = chartInstance->c2_X1[c2_i3];
      }

      c2_dv17[7] = chartInstance->c2_U1;
      for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
        c2_alpha1 = 0.0;
        c2_ldc = 0;
        for (c2_i = 0; c2_i < 8; c2_i++) {
          c2_alpha1 += chartInstance->c2_CD[c2_ldc + c2_i3] * c2_dv17[c2_i];
          c2_ldc += 2;
        }

        c2_c_Dy[c2_i3] = chartInstance->c2_Y1[c2_i3] - c2_alpha1;
      }

      for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
        c2_dv18[c2_i3] = chartInstance->c2_X1[c2_i3];
      }

      c2_dv18[7] = chartInstance->c2_U1;
      for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
        c2_dv18[c2_i3 + 8] = c2_c_Dy[c2_i3];
      }

      c2_b_rls_ew(chartInstance, c2_dv18, chartInstance->c2_X,
                  chartInstance->c2_ABK, chartInstance->c2_Pabk, 0.999);
    }

    if ((chartInstance->c2_start < 0.5) || (chartInstance->c2_k >= 150.0)) {
      for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
        chartInstance->c2_X1[c2_i3] = chartInstance->c2_X[c2_i3];
      }
    }

    for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
      chartInstance->c2_Y1[c2_i3] = c2_Fs[c2_i3];
    }

    chartInstance->c2_U1 = c2_U;
  }

  c2_ON = c2_rdivide(chartInstance, 1.0, c2_b_Du);
  c2_i3 = 0;
  for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_A[c2_i + c2_i3] = chartInstance->c2_ABK[c2_i + c2_i3];
    }

    c2_b_a[c2_ldc] = chartInstance->c2_ABK[c2_ldc + 49];
    c2_i3 += 7;
  }

  for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
    for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
      c2_C[c2_i3 + (c2_ldc << 1)] = 0.0;
      for (c2_i = 0; c2_i < 2; c2_i++) {
        c2_C[c2_i3 + (c2_ldc << 1)] += c2_b_Dy[c2_i3 + (c2_i << 1)] *
          chartInstance->c2_CD[c2_i + (c2_ldc << 1)];
      }
    }

    c2_Fs[c2_i3] = 0.0;
    for (c2_ldc = 0; c2_ldc < 2; c2_ldc++) {
      c2_Fs[c2_i3] += c2_b_Dy[c2_i3 + (c2_ldc << 1)] * chartInstance->c2_CD[14 +
        c2_ldc];
    }
  }

  c2_U = c2_rdivide(chartInstance, 1.0, c2_b_Du);
  c2_i3 = 0;
  for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_dv19[c2_i + c2_i3] = chartInstance->c2_ABK[c2_i + c2_i3];
    }

    c2_i3 += 7;
  }

  c2_damp(chartInstance, c2_dv19, 0.05, c2_Wn, c2_Zn);
  c2_b_Du = chartInstance->c2_w[sf_mex_lw_bounds_check((int32_T)
    chartInstance->c2_k, 1, 300) - 1];
  c2_b_y.re = 3.1415926535897931 * (2.0 * (chartInstance->c2_w[(int32_T)
    chartInstance->c2_k - 1] * 0.0));
  c2_b_y.im = 3.1415926535897931 * (2.0 * (chartInstance->c2_w[(int32_T)
    chartInstance->c2_k - 1] * 0.05));
  c2_exp(chartInstance, &c2_b_y);
  c2_d_eye(chartInstance, c2_dv19);
  for (c2_i3 = 0; c2_i3 < 49; c2_i3++) {
    c2_c_y[c2_i3].re = c2_b_y.re * c2_dv19[c2_i3] - c2_A[c2_i3];
    c2_c_y[c2_i3].im = c2_b_y.im * c2_dv19[c2_i3];
  }

  for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
    c2_c_a[c2_i3] = c2_b_a[c2_i3] * c2_ON;
  }

  c2_mldivide(chartInstance, c2_c_y, c2_c_a, c2_b);
  c2_i3 = 0;
  for (c2_ldc = 0; c2_ldc < 7; c2_ldc++) {
    for (c2_i = 0; c2_i < 2; c2_i++) {
      c2_b_C[c2_i + c2_i3].re = c2_C[c2_i + c2_i3];
      c2_b_C[c2_i + c2_i3].im = 0.0;
    }

    c2_i3 += 2;
  }

  for (c2_i3 = 0; c2_i3 < 2; c2_i3++) {
    c2_ON = 0.0;
    c2_alpha1 = 0.0;
    c2_ldc = 0;
    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_ON += c2_b_C[c2_ldc + c2_i3].re * c2_b[c2_i].re - 0.0 * c2_b[c2_i].im;
      c2_alpha1 += c2_b_C[c2_ldc + c2_i3].re * c2_b[c2_i].im + 0.0 * c2_b[c2_i].
        re;
      c2_ldc += 2;
    }

    c2_c_C[c2_i3].re = c2_ON + c2_Fs[c2_i3] * c2_U;
    c2_c_C[c2_i3].im = c2_alpha1;
  }

  c2_abs(chartInstance, c2_c_C, c2_Fs);
  c2_log10(chartInstance, &c2_b_Du);
  if (c2_Fs[0] < 1.0E-5) {
    c2_alpha1 = -100.0;
  } else {
    c2_ON = c2_Fs[0];
    c2_log10(chartInstance, &c2_ON);
    c2_alpha1 = 20.0 * c2_ON;
  }

  if (c2_Fs[1] < 1.0E-5) {
    *c2_Fs2 = -100.0;
  } else {
    c2_ON = c2_Fs[1];
    c2_log10(chartInstance, &c2_ON);
    *c2_Fs2 = 20.0 * c2_ON;
  }

  if (chartInstance->c2_k >= 300.0) {
    chartInstance->c2_start = 0.0;
    chartInstance->c2_k = 1.0;
  } else {
    chartInstance->c2_k++;
  }

  *c2_Ws = c2_b_Du;
  *c2_Fs1 = c2_alpha1;
  for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
    (*c2_b_Wn)[c2_i3] = c2_Wn[c2_i3];
    (*c2_b_Zn)[c2_i3] = c2_Zn[c2_i3];
  }
}

static void initSimStructsc2_sim1_rlti_varx_gust2
  (SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber)
{
}

const mxArray *sf_c2_sim1_rlti_varx_gust2_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo;
  c2_ResolvedFunctionInfo c2_info[106];
  const mxArray *c2_m0 = NULL;
  int32_T c2_i4;
  c2_ResolvedFunctionInfo *c2_r0;
  c2_nameCaptureInfo = NULL;
  c2_info_helper(c2_info);
  c2_b_info_helper(c2_info);
  sf_mex_assign(&c2_m0, sf_mex_createstruct("nameCaptureInfo", 1, 106), FALSE);
  for (c2_i4 = 0; c2_i4 < 106; c2_i4++) {
    c2_r0 = &c2_info[c2_i4];
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->context)), "context", "nameCaptureInfo",
                    c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c2_r0->name)), "name", "nameCaptureInfo", c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c2_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->resolved)), "resolved", "nameCaptureInfo",
                    c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c2_i4);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c2_i4);
  }

  sf_mex_assign(&c2_nameCaptureInfo, c2_m0, FALSE);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[106])
{
  c2_info[0].context = "";
  c2_info[0].name = "mrdivide";
  c2_info[0].dominantType = "double";
  c2_info[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[0].fileTimeLo = 1310137456U;
  c2_info[0].fileTimeHi = 0U;
  c2_info[0].mFileTimeLo = 1289519692U;
  c2_info[0].mFileTimeHi = 0U;
  c2_info[1].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[1].name = "rdivide";
  c2_info[1].dominantType = "double";
  c2_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[1].fileTimeLo = 1286818844U;
  c2_info[1].fileTimeHi = 0U;
  c2_info[1].mFileTimeLo = 0U;
  c2_info[1].mFileTimeHi = 0U;
  c2_info[2].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[2].name = "eml_div";
  c2_info[2].dominantType = "double";
  c2_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[2].fileTimeLo = 1305318000U;
  c2_info[2].fileTimeHi = 0U;
  c2_info[2].mFileTimeLo = 0U;
  c2_info[2].mFileTimeHi = 0U;
  c2_info[3].context = "";
  c2_info[3].name = "mtimes";
  c2_info[3].dominantType = "double";
  c2_info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[3].fileTimeLo = 1289519692U;
  c2_info[3].fileTimeHi = 0U;
  c2_info[3].mFileTimeLo = 0U;
  c2_info[3].mFileTimeHi = 0U;
  c2_info[4].context = "";
  c2_info[4].name = "eye";
  c2_info[4].dominantType = "double";
  c2_info[4].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m";
  c2_info[4].fileTimeLo = 1286818688U;
  c2_info[4].fileTimeHi = 0U;
  c2_info[4].mFileTimeLo = 0U;
  c2_info[4].mFileTimeHi = 0U;
  c2_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[5].name = "eml_assert_valid_size_arg";
  c2_info[5].dominantType = "double";
  c2_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[5].fileTimeLo = 1286818694U;
  c2_info[5].fileTimeHi = 0U;
  c2_info[5].mFileTimeLo = 0U;
  c2_info[5].mFileTimeHi = 0U;
  c2_info[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral";
  c2_info[6].name = "isinf";
  c2_info[6].dominantType = "double";
  c2_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m";
  c2_info[6].fileTimeLo = 1286818760U;
  c2_info[6].fileTimeHi = 0U;
  c2_info[6].mFileTimeLo = 0U;
  c2_info[6].mFileTimeHi = 0U;
  c2_info[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[7].name = "eml_index_class";
  c2_info[7].dominantType = "";
  c2_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[7].fileTimeLo = 1286818778U;
  c2_info[7].fileTimeHi = 0U;
  c2_info[7].mFileTimeLo = 0U;
  c2_info[7].mFileTimeHi = 0U;
  c2_info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[8].name = "intmax";
  c2_info[8].dominantType = "char";
  c2_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[8].fileTimeLo = 1286818756U;
  c2_info[8].fileTimeHi = 0U;
  c2_info[8].mFileTimeLo = 0U;
  c2_info[8].mFileTimeHi = 0U;
  c2_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[9].name = "eml_is_float_class";
  c2_info[9].dominantType = "char";
  c2_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c2_info[9].fileTimeLo = 1286818782U;
  c2_info[9].fileTimeHi = 0U;
  c2_info[9].mFileTimeLo = 0U;
  c2_info[9].mFileTimeHi = 0U;
  c2_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[10].name = "min";
  c2_info[10].dominantType = "double";
  c2_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[10].fileTimeLo = 1308747330U;
  c2_info[10].fileTimeHi = 0U;
  c2_info[10].mFileTimeLo = 0U;
  c2_info[10].mFileTimeHi = 0U;
  c2_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[11].name = "eml_min_or_max";
  c2_info[11].dominantType = "char";
  c2_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  c2_info[11].fileTimeLo = 1303146212U;
  c2_info[11].fileTimeHi = 0U;
  c2_info[11].mFileTimeLo = 0U;
  c2_info[11].mFileTimeHi = 0U;
  c2_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[12].name = "eml_scalar_eg";
  c2_info[12].dominantType = "double";
  c2_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[12].fileTimeLo = 1286818796U;
  c2_info[12].fileTimeHi = 0U;
  c2_info[12].mFileTimeLo = 0U;
  c2_info[12].mFileTimeHi = 0U;
  c2_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[13].name = "eml_scalexp_alloc";
  c2_info[13].dominantType = "double";
  c2_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c2_info[13].fileTimeLo = 1286818796U;
  c2_info[13].fileTimeHi = 0U;
  c2_info[13].mFileTimeLo = 0U;
  c2_info[13].mFileTimeHi = 0U;
  c2_info[14].context = "";
  c2_info[14].name = "logspace";
  c2_info[14].dominantType = "double";
  c2_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[14].fileTimeLo = 1305318000U;
  c2_info[14].fileTimeHi = 0U;
  c2_info[14].mFileTimeLo = 0U;
  c2_info[14].mFileTimeHi = 0U;
  c2_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[15].name = "eml_scalar_floor";
  c2_info[15].dominantType = "double";
  c2_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c2_info[15].fileTimeLo = 1286818726U;
  c2_info[15].fileTimeHi = 0U;
  c2_info[15].mFileTimeLo = 0U;
  c2_info[15].mFileTimeHi = 0U;
  c2_info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[16].name = "abs";
  c2_info[16].dominantType = "double";
  c2_info[16].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[16].fileTimeLo = 1286818694U;
  c2_info[16].fileTimeHi = 0U;
  c2_info[16].mFileTimeLo = 0U;
  c2_info[16].mFileTimeHi = 0U;
  c2_info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[17].name = "eml_scalar_abs";
  c2_info[17].dominantType = "double";
  c2_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c2_info[17].fileTimeLo = 1286818712U;
  c2_info[17].fileTimeHi = 0U;
  c2_info[17].mFileTimeLo = 0U;
  c2_info[17].mFileTimeHi = 0U;
  c2_info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[18].name = "eps";
  c2_info[18].dominantType = "char";
  c2_info[18].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c2_info[18].fileTimeLo = 1286818686U;
  c2_info[18].fileTimeHi = 0U;
  c2_info[18].mFileTimeLo = 0U;
  c2_info[18].mFileTimeHi = 0U;
  c2_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[19].name = "eml_warning";
  c2_info[19].dominantType = "char";
  c2_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c2_info[19].fileTimeLo = 1286818802U;
  c2_info[19].fileTimeHi = 0U;
  c2_info[19].mFileTimeLo = 0U;
  c2_info[19].mFileTimeHi = 0U;
  c2_info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c2_info[20].name = "linspace";
  c2_info[20].dominantType = "double";
  c2_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c2_info[20].fileTimeLo = 1286818762U;
  c2_info[20].fileTimeHi = 0U;
  c2_info[20].mFileTimeLo = 0U;
  c2_info[20].mFileTimeHi = 0U;
  c2_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c2_info[21].name = "realmax";
  c2_info[21].dominantType = "char";
  c2_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c2_info[21].fileTimeLo = 1286818766U;
  c2_info[21].fileTimeHi = 0U;
  c2_info[21].mFileTimeLo = 0U;
  c2_info[21].mFileTimeHi = 0U;
  c2_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c2_info[22].name = "mpower";
  c2_info[22].dominantType = "double";
  c2_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c2_info[22].fileTimeLo = 1286818842U;
  c2_info[22].fileTimeHi = 0U;
  c2_info[22].mFileTimeLo = 0U;
  c2_info[22].mFileTimeHi = 0U;
  c2_info[23].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c2_info[23].name = "power";
  c2_info[23].dominantType = "double";
  c2_info[23].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  c2_info[23].fileTimeLo = 1294067944U;
  c2_info[23].fileTimeHi = 0U;
  c2_info[23].mFileTimeLo = 0U;
  c2_info[23].mFileTimeHi = 0U;
  c2_info[24].context = "";
  c2_info[24].name = "diag";
  c2_info[24].dominantType = "double";
  c2_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c2_info[24].fileTimeLo = 1286818686U;
  c2_info[24].fileTimeHi = 0U;
  c2_info[24].mFileTimeLo = 0U;
  c2_info[24].mFileTimeHi = 0U;
  c2_info[25].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c2_info[25].name = "eml_index_plus";
  c2_info[25].dominantType = "int32";
  c2_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[25].fileTimeLo = 1286818778U;
  c2_info[25].fileTimeHi = 0U;
  c2_info[25].mFileTimeLo = 0U;
  c2_info[25].mFileTimeHi = 0U;
  c2_info[26].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c2_info[26].name = "eml_index_times";
  c2_info[26].dominantType = "int32";
  c2_info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[26].fileTimeLo = 1286818780U;
  c2_info[26].fileTimeHi = 0U;
  c2_info[26].mFileTimeLo = 0U;
  c2_info[26].mFileTimeHi = 0U;
  c2_info[27].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c2_info[27].name = "eml_index_minus";
  c2_info[27].dominantType = "int32";
  c2_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[27].fileTimeLo = 1286818778U;
  c2_info[27].fileTimeHi = 0U;
  c2_info[27].mFileTimeLo = 0U;
  c2_info[27].mFileTimeHi = 0U;
  c2_info[28].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[28].name = "eml_xgemm";
  c2_info[28].dominantType = "int32";
  c2_info[28].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[28].fileTimeLo = 1299076772U;
  c2_info[28].fileTimeHi = 0U;
  c2_info[28].mFileTimeLo = 0U;
  c2_info[28].mFileTimeHi = 0U;
  c2_info[29].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[29].name = "eml_blas_inline";
  c2_info[29].dominantType = "";
  c2_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[29].fileTimeLo = 1299076768U;
  c2_info[29].fileTimeHi = 0U;
  c2_info[29].mFileTimeLo = 0U;
  c2_info[29].mFileTimeHi = 0U;
  c2_info[30].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[30].name = "eml_refblas_xgemm";
  c2_info[30].dominantType = "int32";
  c2_info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c2_info[30].fileTimeLo = 1299076774U;
  c2_info[30].fileTimeHi = 0U;
  c2_info[30].mFileTimeLo = 0U;
  c2_info[30].mFileTimeHi = 0U;
  c2_info[31].context = "";
  c2_info[31].name = "sqrt";
  c2_info[31].dominantType = "double";
  c2_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[31].fileTimeLo = 1286818752U;
  c2_info[31].fileTimeHi = 0U;
  c2_info[31].mFileTimeLo = 0U;
  c2_info[31].mFileTimeHi = 0U;
  c2_info[32].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[32].name = "eml_error";
  c2_info[32].dominantType = "char";
  c2_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c2_info[32].fileTimeLo = 1305318000U;
  c2_info[32].fileTimeHi = 0U;
  c2_info[32].mFileTimeLo = 0U;
  c2_info[32].mFileTimeHi = 0U;
  c2_info[33].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[33].name = "eml_scalar_sqrt";
  c2_info[33].dominantType = "double";
  c2_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c2_info[33].fileTimeLo = 1286818738U;
  c2_info[33].fileTimeHi = 0U;
  c2_info[33].mFileTimeLo = 0U;
  c2_info[33].mFileTimeHi = 0U;
  c2_info[34].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[34].name = "eml_xdotu";
  c2_info[34].dominantType = "int32";
  c2_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c2_info[34].fileTimeLo = 1299076772U;
  c2_info[34].fileTimeHi = 0U;
  c2_info[34].mFileTimeLo = 0U;
  c2_info[34].mFileTimeHi = 0U;
  c2_info[35].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c2_info[35].name = "eml_xdot";
  c2_info[35].dominantType = "int32";
  c2_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c2_info[35].fileTimeLo = 1299076772U;
  c2_info[35].fileTimeHi = 0U;
  c2_info[35].mFileTimeLo = 0U;
  c2_info[35].mFileTimeHi = 0U;
  c2_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c2_info[36].name = "eml_refblas_xdot";
  c2_info[36].dominantType = "int32";
  c2_info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[36].fileTimeLo = 1299076772U;
  c2_info[36].fileTimeHi = 0U;
  c2_info[36].mFileTimeLo = 0U;
  c2_info[36].mFileTimeHi = 0U;
  c2_info[37].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[37].name = "eml_refblas_xdotx";
  c2_info[37].dominantType = "int32";
  c2_info[37].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[37].fileTimeLo = 1299076774U;
  c2_info[37].fileTimeHi = 0U;
  c2_info[37].mFileTimeLo = 0U;
  c2_info[37].mFileTimeHi = 0U;
  c2_info[38].context = "";
  c2_info[38].name = "norm";
  c2_info[38].dominantType = "double";
  c2_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m";
  c2_info[38].fileTimeLo = 1286818826U;
  c2_info[38].fileTimeHi = 0U;
  c2_info[38].mFileTimeLo = 0U;
  c2_info[38].mFileTimeHi = 0U;
  c2_info[39].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm";
  c2_info[39].name = "eml_xnrm2";
  c2_info[39].dominantType = "int32";
  c2_info[39].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m";
  c2_info[39].fileTimeLo = 1299076776U;
  c2_info[39].fileTimeHi = 0U;
  c2_info[39].mFileTimeLo = 0U;
  c2_info[39].mFileTimeHi = 0U;
  c2_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m";
  c2_info[40].name = "eml_refblas_xnrm2";
  c2_info[40].dominantType = "int32";
  c2_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c2_info[40].fileTimeLo = 1299076784U;
  c2_info[40].fileTimeHi = 0U;
  c2_info[40].mFileTimeLo = 0U;
  c2_info[40].mFileTimeHi = 0U;
  c2_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c2_info[41].name = "realmin";
  c2_info[41].dominantType = "char";
  c2_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c2_info[41].fileTimeLo = 1286818768U;
  c2_info[41].fileTimeHi = 0U;
  c2_info[41].mFileTimeLo = 0U;
  c2_info[41].mFileTimeHi = 0U;
  c2_info[42].context = "";
  c2_info[42].name = "colon";
  c2_info[42].dominantType = "double";
  c2_info[42].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c2_info[42].fileTimeLo = 1286818838U;
  c2_info[42].fileTimeHi = 0U;
  c2_info[42].mFileTimeLo = 0U;
  c2_info[42].mFileTimeHi = 0U;
  c2_info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c2_info[43].name = "isfinite";
  c2_info[43].dominantType = "double";
  c2_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c2_info[43].fileTimeLo = 1286818758U;
  c2_info[43].fileTimeHi = 0U;
  c2_info[43].mFileTimeLo = 0U;
  c2_info[43].mFileTimeHi = 0U;
  c2_info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c2_info[44].name = "isnan";
  c2_info[44].dominantType = "double";
  c2_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  c2_info[44].fileTimeLo = 1286818760U;
  c2_info[44].fileTimeHi = 0U;
  c2_info[44].mFileTimeLo = 0U;
  c2_info[44].mFileTimeHi = 0U;
  c2_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c2_info[45].name = "floor";
  c2_info[45].dominantType = "double";
  c2_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c2_info[45].fileTimeLo = 1286818742U;
  c2_info[45].fileTimeHi = 0U;
  c2_info[45].mFileTimeLo = 0U;
  c2_info[45].mFileTimeHi = 0U;
  c2_info[46].context = "";
  c2_info[46].name = "eig";
  c2_info[46].dominantType = "double";
  c2_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c2_info[46].fileTimeLo = 1305318000U;
  c2_info[46].fileTimeHi = 0U;
  c2_info[46].mFileTimeLo = 0U;
  c2_info[46].mFileTimeHi = 0U;
  c2_info[47].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c2_info[47].name = "eml_xgeev";
  c2_info[47].dominantType = "double";
  c2_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c2_info[47].fileTimeLo = 1286818804U;
  c2_info[47].fileTimeHi = 0U;
  c2_info[47].mFileTimeLo = 0U;
  c2_info[47].mFileTimeHi = 0U;
  c2_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c2_info[48].name = "eml_lapack_xgeev";
  c2_info[48].dominantType = "double";
  c2_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c2_info[48].fileTimeLo = 1301328468U;
  c2_info[48].fileTimeHi = 0U;
  c2_info[48].mFileTimeLo = 0U;
  c2_info[48].mFileTimeHi = 0U;
  c2_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c2_info[49].name = "eml_matlab_zggev";
  c2_info[49].dominantType = "double";
  c2_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[49].fileTimeLo = 1286818818U;
  c2_info[49].fileTimeHi = 0U;
  c2_info[49].mFileTimeLo = 0U;
  c2_info[49].mFileTimeHi = 0U;
  c2_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[50].name = "eml_matlab_zlangeM";
  c2_info[50].dominantType = "double";
  c2_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c2_info[50].fileTimeLo = 1286818820U;
  c2_info[50].fileTimeHi = 0U;
  c2_info[50].mFileTimeLo = 0U;
  c2_info[50].mFileTimeHi = 0U;
  c2_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c2_info[51].name = "eml_dlapy2";
  c2_info[51].dominantType = "double";
  c2_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m";
  c2_info[51].fileTimeLo = 1286818698U;
  c2_info[51].fileTimeHi = 0U;
  c2_info[51].mFileTimeLo = 0U;
  c2_info[51].mFileTimeHi = 0U;
  c2_info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c2_info[52].name = "eml_guarded_nan";
  c2_info[52].dominantType = "char";
  c2_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m";
  c2_info[52].fileTimeLo = 1286818776U;
  c2_info[52].fileTimeHi = 0U;
  c2_info[52].mFileTimeLo = 0U;
  c2_info[52].mFileTimeHi = 0U;
  c2_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[53].name = "eml_matlab_zlascl";
  c2_info[53].dominantType = "double";
  c2_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m";
  c2_info[53].fileTimeLo = 1286818822U;
  c2_info[53].fileTimeHi = 0U;
  c2_info[53].mFileTimeLo = 0U;
  c2_info[53].mFileTimeHi = 0U;
  c2_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[54].name = "eml_matlab_zggbal";
  c2_info[54].dominantType = "double";
  c2_info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m";
  c2_info[54].fileTimeLo = 1286818818U;
  c2_info[54].fileTimeHi = 0U;
  c2_info[54].mFileTimeLo = 0U;
  c2_info[54].mFileTimeHi = 0U;
  c2_info[55].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[55].name = "eml_matlab_zgghrd";
  c2_info[55].dominantType = "int32";
  c2_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c2_info[55].fileTimeLo = 1286818820U;
  c2_info[55].fileTimeHi = 0U;
  c2_info[55].mFileTimeLo = 0U;
  c2_info[55].mFileTimeHi = 0U;
  c2_info[56].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c2_info[56].name = "eml_matlab_zlartg";
  c2_info[56].dominantType = "double";
  c2_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c2_info[56].fileTimeLo = 1286818822U;
  c2_info[56].fileTimeHi = 0U;
  c2_info[56].mFileTimeLo = 0U;
  c2_info[56].mFileTimeHi = 0U;
  c2_info[57].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c2_info[57].name = "fix";
  c2_info[57].dominantType = "double";
  c2_info[57].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c2_info[57].fileTimeLo = 1286818742U;
  c2_info[57].fileTimeHi = 0U;
  c2_info[57].mFileTimeLo = 0U;
  c2_info[57].mFileTimeHi = 0U;
  c2_info[58].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c2_info[58].name = "eml_scalar_fix";
  c2_info[58].dominantType = "double";
  c2_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m";
  c2_info[58].fileTimeLo = 1286818726U;
  c2_info[58].fileTimeHi = 0U;
  c2_info[58].mFileTimeLo = 0U;
  c2_info[58].mFileTimeHi = 0U;
  c2_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c2_info[59].name = "eml_zrot_rows";
  c2_info[59].dominantType = "int32";
  c2_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c2_info[59].fileTimeLo = 1286818826U;
  c2_info[59].fileTimeHi = 0U;
  c2_info[59].mFileTimeLo = 0U;
  c2_info[59].mFileTimeHi = 0U;
  c2_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c2_info[60].name = "eml_conjtimes";
  c2_info[60].dominantType = "double";
  c2_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_conjtimes.m";
  c2_info[60].fileTimeLo = 1286818696U;
  c2_info[60].fileTimeHi = 0U;
  c2_info[60].mFileTimeLo = 0U;
  c2_info[60].mFileTimeHi = 0U;
  c2_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c2_info[61].name = "eml_zrot_cols";
  c2_info[61].dominantType = "int32";
  c2_info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m";
  c2_info[61].fileTimeLo = 1286818826U;
  c2_info[61].fileTimeHi = 0U;
  c2_info[61].mFileTimeLo = 0U;
  c2_info[61].mFileTimeHi = 0U;
  c2_info[62].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c2_info[62].name = "eml_matlab_zhgeqz";
  c2_info[62].dominantType = "int32";
  c2_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c2_info[62].fileTimeLo = 1286818820U;
  c2_info[62].fileTimeHi = 0U;
  c2_info[62].mFileTimeLo = 0U;
  c2_info[62].mFileTimeHi = 0U;
  c2_info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c2_info[63].name = "eml_matlab_zlanhs";
  c2_info[63].dominantType = "int32";
  c2_info[63].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m";
  c2_info[63].fileTimeLo = 1286818820U;
  c2_info[63].fileTimeHi = 0U;
  c2_info[63].mFileTimeLo = 0U;
  c2_info[63].mFileTimeHi = 0U;
}

static void c2_b_info_helper(c2_ResolvedFunctionInfo c2_info[106])
{
  c2_info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c2_info[64].name = "mod";
  c2_info[64].dominantType = "int32";
  c2_info[64].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c2_info[64].fileTimeLo = 1286818744U;
  c2_info[64].fileTimeHi = 0U;
  c2_info[64].mFileTimeLo = 0U;
  c2_info[64].mFileTimeHi = 0U;
  c2_info[65].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c2_info[65].name = "eml_scalar_mod";
  c2_info[65].dominantType = "int32";
  c2_info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_mod.m";
  c2_info[65].fileTimeLo = 1286818730U;
  c2_info[65].fileTimeHi = 0U;
  c2_info[65].mFileTimeLo = 0U;
  c2_info[65].mFileTimeHi = 0U;
  c2_info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c2_info[66].name = "eml_guarded_inf";
  c2_info[66].dominantType = "char";
  c2_info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m";
  c2_info[66].fileTimeLo = 1286818776U;
  c2_info[66].fileTimeHi = 0U;
  c2_info[66].mFileTimeLo = 0U;
  c2_info[66].mFileTimeHi = 0U;
  c2_info[67].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c2_info[67].name = "eml_scalar_hypot";
  c2_info[67].dominantType = "double";
  c2_info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m";
  c2_info[67].fileTimeLo = 1286818726U;
  c2_info[67].fileTimeHi = 0U;
  c2_info[67].mFileTimeLo = 0U;
  c2_info[67].mFileTimeHi = 0U;
  c2_info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c2_info[68].name = "isequal";
  c2_info[68].dominantType = "double";
  c2_info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[68].fileTimeLo = 1286818758U;
  c2_info[68].fileTimeHi = 0U;
  c2_info[68].mFileTimeLo = 0U;
  c2_info[68].mFileTimeHi = 0U;
  c2_info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[69].name = "eml_isequal_core";
  c2_info[69].dominantType = "double";
  c2_info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c2_info[69].fileTimeLo = 1286818786U;
  c2_info[69].fileTimeHi = 0U;
  c2_info[69].mFileTimeLo = 0U;
  c2_info[69].mFileTimeHi = 0U;
  c2_info[70].context = "";
  c2_info[70].name = "sort";
  c2_info[70].dominantType = "double";
  c2_info[70].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c2_info[70].fileTimeLo = 1303146208U;
  c2_info[70].fileTimeHi = 0U;
  c2_info[70].mFileTimeLo = 0U;
  c2_info[70].mFileTimeHi = 0U;
  c2_info[71].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c2_info[71].name = "eml_sort";
  c2_info[71].dominantType = "double";
  c2_info[71].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c2_info[71].fileTimeLo = 1305318002U;
  c2_info[71].fileTimeHi = 0U;
  c2_info[71].mFileTimeLo = 0U;
  c2_info[71].mFileTimeHi = 0U;
  c2_info[72].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c2_info[72].name = "eml_nonsingleton_dim";
  c2_info[72].dominantType = "double";
  c2_info[72].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_nonsingleton_dim.m";
  c2_info[72].fileTimeLo = 1286818788U;
  c2_info[72].fileTimeHi = 0U;
  c2_info[72].mFileTimeLo = 0U;
  c2_info[72].mFileTimeHi = 0U;
  c2_info[73].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c2_info[73].name = "eml_assert_valid_dim";
  c2_info[73].dominantType = "double";
  c2_info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_dim.m";
  c2_info[73].fileTimeLo = 1286818694U;
  c2_info[73].fileTimeHi = 0U;
  c2_info[73].mFileTimeLo = 0U;
  c2_info[73].mFileTimeHi = 0U;
  c2_info[74].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c2_info[74].name = "eml_sort_idx";
  c2_info[74].dominantType = "double";
  c2_info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c2_info[74].fileTimeLo = 1305318004U;
  c2_info[74].fileTimeHi = 0U;
  c2_info[74].mFileTimeLo = 0U;
  c2_info[74].mFileTimeHi = 0U;
  c2_info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c2_info[75].name = "eml_size_ispow2";
  c2_info[75].dominantType = "int32";
  c2_info[75].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c2_info[75].fileTimeLo = 1286818798U;
  c2_info[75].fileTimeHi = 0U;
  c2_info[75].mFileTimeLo = 0U;
  c2_info[75].mFileTimeHi = 0U;
  c2_info[76].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c2_info[76].name = "eml_unsigned_class";
  c2_info[76].dominantType = "char";
  c2_info[76].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c2_info[76].fileTimeLo = 1286818800U;
  c2_info[76].fileTimeHi = 0U;
  c2_info[76].mFileTimeLo = 0U;
  c2_info[76].mFileTimeHi = 0U;
  c2_info[77].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c2_info[77].name = "eml_sort_le";
  c2_info[77].dominantType = "int32";
  c2_info[77].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m";
  c2_info[77].fileTimeLo = 1292190510U;
  c2_info[77].fileTimeHi = 0U;
  c2_info[77].mFileTimeLo = 0U;
  c2_info[77].mFileTimeHi = 0U;
  c2_info[78].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m!eml_sort_ascending_le";
  c2_info[78].name = "eml_relop";
  c2_info[78].dominantType = "function_handle";
  c2_info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m";
  c2_info[78].fileTimeLo = 1292190510U;
  c2_info[78].fileTimeHi = 0U;
  c2_info[78].mFileTimeLo = 0U;
  c2_info[78].mFileTimeHi = 0U;
  c2_info[79].context = "";
  c2_info[79].name = "any";
  c2_info[79].dominantType = "double";
  c2_info[79].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c2_info[79].fileTimeLo = 1286818834U;
  c2_info[79].fileTimeHi = 0U;
  c2_info[79].mFileTimeLo = 0U;
  c2_info[79].mFileTimeHi = 0U;
  c2_info[80].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c2_info[80].name = "eml_all_or_any";
  c2_info[80].dominantType = "char";
  c2_info[80].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c2_info[80].fileTimeLo = 1286818694U;
  c2_info[80].fileTimeHi = 0U;
  c2_info[80].mFileTimeLo = 0U;
  c2_info[80].mFileTimeHi = 0U;
  c2_info[81].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c2_info[81].name = "eml_const_nonsingleton_dim";
  c2_info[81].dominantType = "double";
  c2_info[81].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  c2_info[81].fileTimeLo = 1286818696U;
  c2_info[81].fileTimeHi = 0U;
  c2_info[81].mFileTimeLo = 0U;
  c2_info[81].mFileTimeHi = 0U;
  c2_info[82].context = "";
  c2_info[82].name = "log";
  c2_info[82].dominantType = "double";
  c2_info[82].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c2_info[82].fileTimeLo = 1286818742U;
  c2_info[82].fileTimeHi = 0U;
  c2_info[82].mFileTimeLo = 0U;
  c2_info[82].mFileTimeHi = 0U;
  c2_info[83].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c2_info[83].name = "eml_scalar_log";
  c2_info[83].dominantType = "double";
  c2_info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c2_info[83].fileTimeLo = 1286818728U;
  c2_info[83].fileTimeHi = 0U;
  c2_info[83].mFileTimeLo = 0U;
  c2_info[83].mFileTimeHi = 0U;
  c2_info[84].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c2_info[84].name = "eml_scalar_atan2";
  c2_info[84].dominantType = "double";
  c2_info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m";
  c2_info[84].fileTimeLo = 1286818720U;
  c2_info[84].fileTimeHi = 0U;
  c2_info[84].mFileTimeLo = 0U;
  c2_info[84].mFileTimeHi = 0U;
  c2_info[85].context = "";
  c2_info[85].name = "exp";
  c2_info[85].dominantType = "double";
  c2_info[85].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c2_info[85].fileTimeLo = 1286818740U;
  c2_info[85].fileTimeHi = 0U;
  c2_info[85].mFileTimeLo = 0U;
  c2_info[85].mFileTimeHi = 0U;
  c2_info[86].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c2_info[86].name = "eml_scalar_exp";
  c2_info[86].dominantType = "double";
  c2_info[86].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m";
  c2_info[86].fileTimeLo = 1301328464U;
  c2_info[86].fileTimeHi = 0U;
  c2_info[86].mFileTimeLo = 0U;
  c2_info[86].mFileTimeHi = 0U;
  c2_info[87].context = "";
  c2_info[87].name = "mldivide";
  c2_info[87].dominantType = "double";
  c2_info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[87].fileTimeLo = 1310137456U;
  c2_info[87].fileTimeHi = 0U;
  c2_info[87].mFileTimeLo = 1289519690U;
  c2_info[87].mFileTimeHi = 0U;
  c2_info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[88].name = "eml_lusolve";
  c2_info[88].dominantType = "double";
  c2_info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c2_info[88].fileTimeLo = 1305318000U;
  c2_info[88].fileTimeHi = 0U;
  c2_info[88].mFileTimeLo = 0U;
  c2_info[88].mFileTimeHi = 0U;
  c2_info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[89].name = "eml_xgetrf";
  c2_info[89].dominantType = "int32";
  c2_info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c2_info[89].fileTimeLo = 1286818806U;
  c2_info[89].fileTimeHi = 0U;
  c2_info[89].mFileTimeLo = 0U;
  c2_info[89].mFileTimeHi = 0U;
  c2_info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c2_info[90].name = "eml_lapack_xgetrf";
  c2_info[90].dominantType = "int32";
  c2_info[90].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c2_info[90].fileTimeLo = 1286818810U;
  c2_info[90].fileTimeHi = 0U;
  c2_info[90].mFileTimeLo = 0U;
  c2_info[90].mFileTimeHi = 0U;
  c2_info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c2_info[91].name = "eml_matlab_zgetrf";
  c2_info[91].dominantType = "int32";
  c2_info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[91].fileTimeLo = 1302688994U;
  c2_info[91].fileTimeHi = 0U;
  c2_info[91].mFileTimeLo = 0U;
  c2_info[91].mFileTimeHi = 0U;
  c2_info[92].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c2_info[92].name = "intmin";
  c2_info[92].dominantType = "char";
  c2_info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c2_info[92].fileTimeLo = 1286818756U;
  c2_info[92].fileTimeHi = 0U;
  c2_info[92].mFileTimeLo = 0U;
  c2_info[92].mFileTimeHi = 0U;
  c2_info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[93].name = "eml_ixamax";
  c2_info[93].dominantType = "int32";
  c2_info[93].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c2_info[93].fileTimeLo = 1299076770U;
  c2_info[93].fileTimeHi = 0U;
  c2_info[93].mFileTimeLo = 0U;
  c2_info[93].mFileTimeHi = 0U;
  c2_info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m!below_threshold";
  c2_info[94].name = "length";
  c2_info[94].dominantType = "double";
  c2_info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c2_info[94].fileTimeLo = 1303146206U;
  c2_info[94].fileTimeHi = 0U;
  c2_info[94].mFileTimeLo = 0U;
  c2_info[94].mFileTimeHi = 0U;
  c2_info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c2_info[95].name = "eml_refblas_ixamax";
  c2_info[95].dominantType = "int32";
  c2_info[95].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[95].fileTimeLo = 1299076770U;
  c2_info[95].fileTimeHi = 0U;
  c2_info[95].mFileTimeLo = 0U;
  c2_info[95].mFileTimeHi = 0U;
  c2_info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[96].name = "eml_xcabs1";
  c2_info[96].dominantType = "double";
  c2_info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c2_info[96].fileTimeLo = 1286818706U;
  c2_info[96].fileTimeHi = 0U;
  c2_info[96].mFileTimeLo = 0U;
  c2_info[96].mFileTimeHi = 0U;
  c2_info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[97].name = "eml_xswap";
  c2_info[97].dominantType = "int32";
  c2_info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c2_info[97].fileTimeLo = 1299076778U;
  c2_info[97].fileTimeHi = 0U;
  c2_info[97].mFileTimeLo = 0U;
  c2_info[97].mFileTimeHi = 0U;
  c2_info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c2_info[98].name = "eml_refblas_xswap";
  c2_info[98].dominantType = "int32";
  c2_info[98].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[98].fileTimeLo = 1299076786U;
  c2_info[98].fileTimeHi = 0U;
  c2_info[98].mFileTimeLo = 0U;
  c2_info[98].mFileTimeHi = 0U;
  c2_info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[99].name = "eml_xgeru";
  c2_info[99].dominantType = "int32";
  c2_info[99].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c2_info[99].fileTimeLo = 1299076774U;
  c2_info[99].fileTimeHi = 0U;
  c2_info[99].mFileTimeLo = 0U;
  c2_info[99].mFileTimeHi = 0U;
  c2_info[100].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgeru.m";
  c2_info[100].name = "eml_refblas_xgeru";
  c2_info[100].dominantType = "int32";
  c2_info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c2_info[100].fileTimeLo = 1299076776U;
  c2_info[100].fileTimeHi = 0U;
  c2_info[100].mFileTimeLo = 0U;
  c2_info[100].mFileTimeHi = 0U;
  c2_info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c2_info[101].name = "eml_refblas_xgerx";
  c2_info[101].dominantType = "int32";
  c2_info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[101].fileTimeLo = 1299076778U;
  c2_info[101].fileTimeHi = 0U;
  c2_info[101].mFileTimeLo = 0U;
  c2_info[101].mFileTimeHi = 0U;
  c2_info[102].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[102].name = "eml_xtrsm";
  c2_info[102].dominantType = "int32";
  c2_info[102].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c2_info[102].fileTimeLo = 1299076778U;
  c2_info[102].fileTimeHi = 0U;
  c2_info[102].mFileTimeLo = 0U;
  c2_info[102].mFileTimeHi = 0U;
  c2_info[103].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c2_info[103].name = "eml_refblas_xtrsm";
  c2_info[103].dominantType = "int32";
  c2_info[103].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[103].fileTimeLo = 1299076786U;
  c2_info[103].fileTimeHi = 0U;
  c2_info[103].mFileTimeLo = 0U;
  c2_info[103].mFileTimeHi = 0U;
  c2_info[104].context = "";
  c2_info[104].name = "log10";
  c2_info[104].dominantType = "double";
  c2_info[104].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c2_info[104].fileTimeLo = 1286818744U;
  c2_info[104].fileTimeHi = 0U;
  c2_info[104].mFileTimeLo = 0U;
  c2_info[104].mFileTimeHi = 0U;
  c2_info[105].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c2_info[105].name = "eml_scalar_log10";
  c2_info[105].dominantType = "double";
  c2_info[105].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log10.m";
  c2_info[105].fileTimeLo = 1286818728U;
  c2_info[105].fileTimeHi = 0U;
  c2_info[105].mFileTimeLo = 0U;
  c2_info[105].mFileTimeHi = 0U;
}

static real_T c2_mrdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_A, real_T c2_B)
{
  return c2_A / c2_B;
}

static real_T c2_rdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_x, real_T c2_y)
{
  return c2_x / c2_y;
}

static void c2_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   real_T c2_I[23104])
{
  int32_T c2_i;
  int32_T c2_b_i;
  for (c2_i = 0; c2_i < 23104; c2_i++) {
    c2_I[c2_i] = 0.0;
  }

  c2_i = 0;
  for (c2_b_i = 0; c2_b_i < 152; c2_b_i++) {
    c2_I[c2_i] = 1.0;
    c2_i += 153;
  }
}

static void c2_b_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[100])
{
  int32_T c2_i;
  int32_T c2_b_i;
  for (c2_i = 0; c2_i < 100; c2_i++) {
    c2_I[c2_i] = 0.0;
  }

  c2_i = 0;
  for (c2_b_i = 0; c2_b_i < 10; c2_b_i++) {
    c2_I[c2_i] = 1.0;
    c2_i += 11;
  }
}

static void c2_c_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[64])
{
  int32_T c2_i;
  int32_T c2_b_i;
  for (c2_i = 0; c2_i < 64; c2_i++) {
    c2_I[c2_i] = 0.0;
  }

  c2_i = 0;
  for (c2_b_i = 0; c2_b_i < 8; c2_b_i++) {
    c2_I[c2_i] = 1.0;
    c2_i += 9;
  }
}

static void c2_logspace(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_d1, real_T c2_d2, real_T c2_y[300])
{
  real_T c2_b_y[300];
  real_T c2_delta1;
  real_T c2_delta2;
  int32_T c2_b_k;
  if (muDoubleScalarAbs(c2_d2 - 3.1415926535897931) < 8.8817841970012523E-16) {
    c2_eml_warning(chartInstance);
  }

  c2_b_y[299] = c2_d2;
  c2_b_y[0] = c2_d1;
  if (((c2_d1 < 0.0 != c2_d2 < 0.0) && (muDoubleScalarAbs(c2_d1) >
        8.9884656743115785E+307)) || (muDoubleScalarAbs(c2_d2) >
       8.9884656743115785E+307)) {
    c2_delta1 = c2_d1 / 299.0;
    c2_delta2 = c2_d2 / 299.0;
    for (c2_b_k = 0; c2_b_k < 298; c2_b_k++) {
      c2_b_y[c2_b_k + 1] = (c2_d1 + c2_delta2 * (1.0 + (real_T)c2_b_k)) -
        c2_delta1 * (1.0 + (real_T)c2_b_k);
    }
  } else {
    c2_delta1 = (c2_d2 - c2_d1) / 299.0;
    for (c2_b_k = 0; c2_b_k < 298; c2_b_k++) {
      c2_b_y[c2_b_k + 1] = c2_d1 + (1.0 + (real_T)c2_b_k) * c2_delta1;
    }
  }

  for (c2_b_k = 0; c2_b_k < 300; c2_b_k++) {
    c2_y[c2_b_k] = muDoubleScalarPower(10.0, c2_b_y[c2_b_k]);
  }
}

static void c2_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  int32_T c2_i5;
  static char_T c2_varargin_1[32] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'p', 'i', 'S', 'p', 'e', 'c', 'i', 'a', 'l',
    'S', 'u', 'p', 'p', 'o', 'r', 't', 'e', 'd' };

  char_T c2_u[32];
  const mxArray *c2_y = NULL;
  for (c2_i5 = 0; c2_i5 < 32; c2_i5++) {
    c2_u[c2_i5] = c2_varargin_1[c2_i5];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 32), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static void c2_b_rdivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_x, real_T c2_y[2], real_T c2_z[2])
{
  int32_T c2_i6;
  for (c2_i6 = 0; c2_i6 < 2; c2_i6++) {
    c2_z[c2_i6] = c2_x / c2_y[c2_i6];
  }
}

static void c2_diag(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    real_T c2_v[2], real_T c2_d[4])
{
  int32_T c2_j;
  int32_T c2_b_j;
  for (c2_j = 0; c2_j < 4; c2_j++) {
    c2_d[c2_j] = 0.0;
  }

  c2_j = 0;
  for (c2_b_j = 0; c2_b_j < 2; c2_b_j++) {
    c2_d[c2_j] = c2_v[c2_b_j];
    c2_j += 3;
  }
}

static void c2_eml_error(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance)
{
  int32_T c2_i7;
  static char_T c2_varargin_1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 's', 'q', 'r', 't', '_', 'd', 'o', 'm', 'a',
    'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c2_u[30];
  const mxArray *c2_y = NULL;
  for (c2_i7 = 0; c2_i7 < 30; c2_i7++) {
    c2_u[c2_i7] = c2_varargin_1[c2_i7];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 30), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static void c2_damp(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    real_T c2_a[49], real_T c2_h, real_T c2_wn[7], real_T c2_z[7])
{
  creal_T c2_r[7];
  int32_T c2_i;
  real_T c2_x[7];
  int8_T c2_idx[7];
  int32_T c2_b_k;
  boolean_T c2_y;
  int8_T c2_idx0[7];
  int32_T c2_i2;
  int32_T c2_j;
  int32_T c2_pEnd;
  int32_T c2_p;
  int32_T c2_q;
  int32_T c2_qEnd;
  int32_T c2_kEnd;
  creal_T c2_b_r[7];
  boolean_T c2_exitg1;
  boolean_T c2_b0;
  real_T c2_b_a;
  real_T c2_b;
  c2_eig(chartInstance, c2_a, c2_r);
  for (c2_i = 0; c2_i < 7; c2_i++) {
    c2_x[c2_i] = c2_r[c2_i].re;
    c2_idx[c2_i] = (int8_T)(c2_i + 1);
  }

  for (c2_b_k = 0; c2_b_k < 5; c2_b_k += 2) {
    c2_i = c2_b_k + 1;
    if ((c2_x[c2_b_k] <= c2_x[c2_i]) || muDoubleScalarIsNaN(c2_x[c2_i])) {
      c2_y = TRUE;
    } else {
      c2_y = FALSE;
    }

    if (!c2_y) {
      c2_idx[c2_b_k] = (int8_T)(c2_b_k + 2);
      c2_idx[c2_b_k + 1] = (int8_T)(c2_b_k + 1);
    }
  }

  for (c2_i = 0; c2_i < 7; c2_i++) {
    c2_idx0[c2_i] = 1;
  }

  c2_i = 2;
  while (c2_i < 7) {
    c2_i2 = c2_i << 1;
    c2_j = 1;
    for (c2_pEnd = 1 + c2_i; c2_pEnd < 8; c2_pEnd = c2_qEnd + c2_i) {
      c2_p = c2_j;
      c2_q = c2_pEnd;
      c2_qEnd = c2_j + c2_i2;
      if (c2_qEnd > 8) {
        c2_qEnd = 8;
      }

      c2_b_k = 1;
      c2_kEnd = c2_qEnd - c2_j;
      while (c2_b_k <= c2_kEnd) {
        sf_mex_lw_bounds_check(c2_p, 1, 7);
        sf_mex_lw_bounds_check(c2_q, 1, 7);
        if ((c2_x[c2_idx[c2_p - 1] - 1] <= c2_x[c2_idx[c2_q - 1] - 1]) ||
            muDoubleScalarIsNaN(c2_x[c2_idx[c2_q - 1] - 1])) {
          c2_y = TRUE;
        } else {
          c2_y = FALSE;
        }

        if (c2_y) {
          c2_idx0[sf_mex_lw_bounds_check(c2_b_k, 1, 7) - 1] = c2_idx[c2_p - 1];
          c2_p++;
          if (c2_p == c2_pEnd) {
            while (c2_q < c2_qEnd) {
              c2_b_k++;
              c2_idx0[sf_mex_lw_bounds_check(c2_b_k, 1, 7) - 1] = c2_idx[c2_q -
                1];
              c2_q++;
            }
          }
        } else {
          c2_idx0[sf_mex_lw_bounds_check(c2_b_k, 1, 7) - 1] = c2_idx[c2_q - 1];
          c2_q++;
          if (c2_q == c2_qEnd) {
            while (c2_p < c2_pEnd) {
              c2_b_k++;
              c2_idx0[sf_mex_lw_bounds_check(c2_b_k, 1, 7) - 1] = c2_idx[c2_p -
                1];
              c2_p++;
            }
          }
        }

        c2_b_k++;
      }

      for (c2_b_k = 1; c2_b_k <= c2_kEnd; c2_b_k++) {
        c2_idx[sf_mex_lw_bounds_check((c2_j + c2_b_k) - 1, 1, 7) - 1] =
          c2_idx0[c2_b_k - 1];
      }

      c2_j = c2_qEnd;
    }

    c2_i = c2_i2;
  }

  for (c2_i = 0; c2_i < 7; c2_i++) {
    c2_b_r[c2_i] = c2_r[c2_idx[c2_i] - 1];
  }

  for (c2_i = 0; c2_i < 7; c2_i++) {
    c2_r[c2_i] = c2_b_r[c2_i];
  }

  c2_y = FALSE;
  c2_b_k = 0;
  c2_exitg1 = 0U;
  while ((c2_exitg1 == 0U) && (c2_b_k < 7)) {
    if (((c2_r[c2_b_k].re == 0.0) && (c2_r[c2_b_k].im == 0.0)) ||
        (muDoubleScalarIsNaN(c2_r[c2_b_k].re) || muDoubleScalarIsNaN(c2_r[c2_b_k]
          .im))) {
      c2_b0 = TRUE;
    } else {
      c2_b0 = FALSE;
    }

    if (!c2_b0) {
      c2_y = TRUE;
      c2_exitg1 = 1U;
    } else {
      c2_b_k++;
    }
  }

  if (c2_y == 0) {
    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_z[c2_i] = 0.0;
      c2_wn[c2_i] = 0.0;
    }
  } else {
    for (c2_b_k = 0; c2_b_k < 7; c2_b_k++) {
      c2_b_a = c2_r[c2_b_k].re;
      c2_b = c2_r[c2_b_k].im;
      if ((c2_r[c2_b_k].im == 0.0) && muDoubleScalarIsNaN(c2_r[c2_b_k].re)) {
      } else if ((muDoubleScalarAbs(c2_r[c2_b_k].re) > 8.9884656743115785E+307) ||
                 (muDoubleScalarAbs(c2_r[c2_b_k].im) > 8.9884656743115785E+307))
      {
        c2_b_a = muDoubleScalarAbs(c2_r[c2_b_k].re / 2.0);
        c2_b = muDoubleScalarAbs(c2_r[c2_b_k].im / 2.0);
        if (c2_b_a < c2_b) {
          c2_b_a /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_b_a * c2_b_a + 1.0);
        } else if (c2_b_a > c2_b) {
          c2_b /= c2_b_a;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_b_a * 1.4142135623730951;
          }
        }

        c2_b_a = muDoubleScalarLog(c2_b) + 0.69314718055994529;
        c2_b = muDoubleScalarAtan2(c2_r[c2_b_k].im, c2_r[c2_b_k].re);
      } else {
        c2_b_a = muDoubleScalarAbs(c2_r[c2_b_k].re);
        c2_b = muDoubleScalarAbs(c2_r[c2_b_k].im);
        if (c2_b_a < c2_b) {
          c2_b_a /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_b_a * c2_b_a + 1.0);
        } else if (c2_b_a > c2_b) {
          c2_b /= c2_b_a;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_b_a * 1.4142135623730951;
          }
        }

        c2_b_a = muDoubleScalarLog(c2_b);
        c2_b = muDoubleScalarAtan2(c2_r[c2_b_k].im, c2_r[c2_b_k].re);
      }

      c2_r[c2_b_k].re = c2_b_a;
      c2_r[c2_b_k].im = c2_b;
    }

    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_b_a = c2_r[c2_i].im;
      if (c2_r[c2_i].im == 0.0) {
        c2_r[c2_i].re /= c2_h;
        c2_r[c2_i].im = 0.0;
      } else if (c2_r[c2_i].re == 0.0) {
        c2_r[c2_i].re = 0.0;
        c2_r[c2_i].im = c2_b_a / c2_h;
      } else {
        c2_r[c2_i].re /= c2_h;
        c2_r[c2_i].im = c2_b_a / c2_h;
      }
    }

    for (c2_b_k = 0; c2_b_k < 7; c2_b_k++) {
      c2_b_a = muDoubleScalarAbs(c2_r[c2_b_k].re);
      c2_b = muDoubleScalarAbs(c2_r[c2_b_k].im);
      if (c2_b_a < c2_b) {
        c2_b_a /= c2_b;
        c2_b *= muDoubleScalarSqrt(c2_b_a * c2_b_a + 1.0);
      } else if (c2_b_a > c2_b) {
        c2_b /= c2_b_a;
        c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_b_a;
      } else {
        if (!muDoubleScalarIsNaN(c2_b)) {
          c2_b = c2_b_a * 1.4142135623730951;
        }
      }

      c2_wn[c2_b_k] = c2_b;
    }

    for (c2_i = 0; c2_i < 7; c2_i++) {
      c2_z[c2_i] = -c2_r[c2_i].re / c2_wn[c2_i];
    }
  }
}

static void c2_eig(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   real_T c2_A[49], creal_T c2_V[7])
{
  int32_T c2_ii;
  creal_T c2_b_A[49];
  real_T c2_info;
  real_T c2_anrm;
  boolean_T c2_exitg7;
  real_T c2_a;
  real_T c2_b;
  creal_T c2_beta1[7];
  boolean_T c2_ilascl;
  real_T c2_anrmto;
  boolean_T c2_notdone;
  real_T c2_cfrom1;
  real_T c2_cto1;
  real_T c2_mul;
  int32_T c2_ilo;
  int32_T c2_ihi;
  int32_T c2_exitg2;
  int32_T c2_i;
  int32_T c2_j;
  boolean_T c2_exitg5;
  int32_T c2_nzcount;
  int32_T c2_jrow;
  boolean_T c2_exitg6;
  boolean_T c2_c_A;
  boolean_T c2_guard2 = FALSE;
  creal_T c2_atmp;
  int32_T c2_exitg1;
  boolean_T c2_exitg3;
  boolean_T c2_exitg4;
  boolean_T c2_guard1 = FALSE;
  int32_T c2_jrowm1;
  creal_T c2_s;
  static creal_T c2_dc0 = { 0.0, 0.0 };

  creal_T c2_d_A;
  creal_T c2_e_A;
  for (c2_ii = 0; c2_ii < 49; c2_ii++) {
    c2_b_A[c2_ii].re = c2_A[c2_ii];
    c2_b_A[c2_ii].im = 0.0;
  }

  c2_info = 0.0;
  c2_anrm = 0.0;
  c2_ii = 0;
  c2_exitg7 = 0U;
  while ((c2_exitg7 == 0U) && (c2_ii < 49)) {
    c2_a = muDoubleScalarAbs(c2_b_A[c2_ii].re);
    if (c2_a > 0.0) {
      c2_b = 0.0 / c2_a;
      c2_a *= muDoubleScalarSqrt(c2_b * c2_b + 1.0);
    } else {
      c2_a *= 1.4142135623730951;
    }

    if (muDoubleScalarIsNaN(c2_a)) {
      c2_anrm = rtNaN;
      c2_exitg7 = 1U;
    } else {
      if (c2_a > c2_anrm) {
        c2_anrm = c2_a;
      }

      c2_ii++;
    }
  }

  if (!((!muDoubleScalarIsInf(c2_anrm)) && (!muDoubleScalarIsNaN(c2_anrm)))) {
    for (c2_ii = 0; c2_ii < 7; c2_ii++) {
      c2_V[c2_ii].re = rtNaN;
      c2_V[c2_ii].im = 0.0;
      c2_beta1[c2_ii].re = rtNaN;
      c2_beta1[c2_ii].im = 0.0;
    }
  } else {
    c2_ilascl = FALSE;
    c2_anrmto = c2_anrm;
    if ((c2_anrm > 0.0) && (c2_anrm < 6.7178761075670888E-139)) {
      c2_anrmto = 6.7178761075670888E-139;
      c2_ilascl = TRUE;
    } else {
      if (c2_anrm > 1.4885657073574029E+138) {
        c2_anrmto = 1.4885657073574029E+138;
        c2_ilascl = TRUE;
      }
    }

    if (c2_ilascl) {
      c2_a = c2_anrm;
      c2_b = c2_anrmto;
      c2_notdone = TRUE;
      while (c2_notdone) {
        c2_cfrom1 = c2_a * 2.0041683600089728E-292;
        c2_cto1 = c2_b / 4.9896007738368E+291;
        if ((c2_cfrom1 > c2_b) && (c2_b != 0.0)) {
          c2_mul = 2.0041683600089728E-292;
          c2_a = c2_cfrom1;
        } else if (c2_cto1 > c2_a) {
          c2_mul = 4.9896007738368E+291;
          c2_b = c2_cto1;
        } else {
          c2_mul = c2_b / c2_a;
          c2_notdone = FALSE;
        }

        for (c2_ii = 0; c2_ii < 49; c2_ii++) {
          c2_b_A[c2_ii].re *= c2_mul;
          c2_b_A[c2_ii].im *= c2_mul;
        }
      }
    }

    c2_ilo = 1;
    c2_ihi = 7;
    do {
      c2_exitg2 = 0U;
      c2_i = 0;
      c2_j = 0;
      c2_notdone = FALSE;
      c2_ii = c2_ihi;
      c2_exitg5 = 0U;
      while ((c2_exitg5 == 0U) && (c2_ii > 0)) {
        c2_nzcount = 0;
        c2_i = c2_ii;
        c2_j = c2_ihi;
        c2_jrow = 1;
        c2_exitg6 = 0U;
        while ((c2_exitg6 == 0U) && (c2_jrow <= c2_ihi)) {
          c2_c_A = ((c2_b_A[(sf_mex_lw_bounds_check(c2_ii, 1, 7) + 7 * (c2_jrow
            - 1)) - 1].re != 0.0) || (c2_b_A[(sf_mex_lw_bounds_check(c2_ii, 1, 7)
                      + 7 * (c2_jrow - 1)) - 1].im != 0.0));
          c2_guard2 = FALSE;
          if (c2_c_A || (c2_ii == c2_jrow)) {
            if (c2_nzcount == 0) {
              c2_j = c2_jrow;
              c2_nzcount = 1;
              c2_guard2 = TRUE;
            } else {
              c2_nzcount = 2;
              c2_exitg6 = 1U;
            }
          } else {
            c2_guard2 = TRUE;
          }

          if (c2_guard2 == TRUE) {
            c2_jrow++;
          }
        }

        if (c2_nzcount < 2) {
          c2_notdone = TRUE;
          c2_exitg5 = 1U;
        } else {
          c2_ii--;
        }
      }

      if (!c2_notdone) {
        c2_exitg2 = 2U;
      } else {
        if (c2_i != c2_ihi) {
          for (c2_ii = 0; c2_ii < 7; c2_ii++) {
            c2_atmp = c2_b_A[(sf_mex_lw_bounds_check(c2_i, 1, 7) + 7 * c2_ii) -
              1];
            c2_b_A[(c2_i + 7 * c2_ii) - 1] = c2_b_A[(sf_mex_lw_bounds_check
              (c2_ihi, 1, 7) + 7 * c2_ii) - 1];
            c2_b_A[(c2_ihi + 7 * c2_ii) - 1] = c2_atmp;
          }
        }

        if (c2_j != c2_ihi) {
          for (c2_ii = 0; c2_ii + 1 <= c2_ihi; c2_ii++) {
            c2_atmp = c2_b_A[c2_ii + 7 * (sf_mex_lw_bounds_check(c2_j, 1, 7) - 1)];
            c2_b_A[c2_ii + 7 * (c2_j - 1)] = c2_b_A[c2_ii + 7 * (c2_ihi - 1)];
            c2_b_A[c2_ii + 7 * (c2_ihi - 1)] = c2_atmp;
          }
        }

        sf_mex_lw_bounds_check(c2_ihi, 1, 7);
        c2_ihi--;
        if (c2_ihi == 1) {
          c2_exitg2 = 1U;
        }
      }
    } while (c2_exitg2 == 0U);

    if (c2_exitg2 == 1U) {
    } else {
      do {
        c2_exitg1 = 0U;
        c2_i = 0;
        c2_j = 0;
        c2_notdone = FALSE;
        c2_jrow = c2_ilo;
        c2_exitg3 = 0U;
        while ((c2_exitg3 == 0U) && (c2_jrow <= c2_ihi)) {
          c2_nzcount = 0;
          c2_i = c2_ihi;
          c2_j = c2_jrow;
          c2_ii = c2_ilo;
          c2_exitg4 = 0U;
          while ((c2_exitg4 == 0U) && (c2_ii <= c2_ihi)) {
            c2_c_A = ((c2_b_A[(c2_ii + 7 * (sf_mex_lw_bounds_check(c2_jrow, 1, 7)
              - 1)) - 1].re != 0.0) || (c2_b_A[(c2_ii + 7 *
                        (sf_mex_lw_bounds_check(c2_jrow, 1, 7) - 1)) - 1].im !=
                       0.0));
            c2_guard1 = FALSE;
            if (c2_c_A || (c2_ii == c2_jrow)) {
              if (c2_nzcount == 0) {
                c2_i = c2_ii;
                c2_nzcount = 1;
                c2_guard1 = TRUE;
              } else {
                c2_nzcount = 2;
                c2_exitg4 = 1U;
              }
            } else {
              c2_guard1 = TRUE;
            }

            if (c2_guard1 == TRUE) {
              c2_ii++;
            }
          }

          if (c2_nzcount < 2) {
            c2_notdone = TRUE;
            c2_exitg3 = 1U;
          } else {
            c2_jrow++;
          }
        }

        if (!c2_notdone) {
          c2_exitg1 = 1U;
        } else {
          if (c2_i != c2_ilo) {
            for (c2_ii = c2_ilo - 1; c2_ii + 1 < 8; c2_ii++) {
              c2_atmp = c2_b_A[(sf_mex_lw_bounds_check(c2_i, 1, 7) + 7 * c2_ii)
                - 1];
              c2_b_A[(c2_i + 7 * c2_ii) - 1] = c2_b_A[(sf_mex_lw_bounds_check
                (c2_ilo, 1, 7) + 7 * c2_ii) - 1];
              c2_b_A[(c2_ilo + 7 * c2_ii) - 1] = c2_atmp;
            }
          }

          if (c2_j != c2_ilo) {
            for (c2_ii = 0; c2_ii + 1 <= c2_ihi; c2_ii++) {
              c2_atmp = c2_b_A[c2_ii + 7 * (sf_mex_lw_bounds_check(c2_j, 1, 7) -
                1)];
              c2_b_A[c2_ii + 7 * (c2_j - 1)] = c2_b_A[c2_ii + 7 *
                (sf_mex_lw_bounds_check(c2_ilo, 1, 7) - 1)];
              c2_b_A[c2_ii + 7 * (c2_ilo - 1)] = c2_atmp;
            }
          }

          sf_mex_lw_bounds_check(c2_ilo, 1, 7);
          c2_ilo++;
          if (c2_ilo == c2_ihi) {
            c2_exitg1 = 1U;
          }
        }
      } while (c2_exitg1 == 0U);
    }

    if (!(c2_ihi < c2_ilo + 2)) {
      c2_ii = c2_ilo - 1;
      while (c2_ii + 1 < c2_ihi - 1) {
        c2_nzcount = c2_ii + 1;
        c2_jrow = c2_ihi - 1;
        while (c2_jrow + 1 > c2_nzcount + 1) {
          c2_jrowm1 = c2_jrow - 1;
          c2_eml_matlab_zlartg(chartInstance, c2_b_A[c2_jrowm1 + 7 * c2_ii],
                               c2_b_A[c2_jrow + 7 * c2_ii], &c2_a, &c2_s,
                               &c2_atmp);
          c2_b_A[c2_jrowm1 + 7 * c2_ii] = c2_atmp;
          c2_b_A[c2_jrow + 7 * c2_ii] = c2_dc0;
          for (c2_j = c2_nzcount; c2_j + 1 <= c2_ihi; c2_j++) {
            c2_atmp.re = c2_a * c2_b_A[c2_jrowm1 + 7 * c2_j].re;
            c2_atmp.im = c2_a * c2_b_A[c2_jrowm1 + 7 * c2_j].im;
            c2_b = c2_s.re * c2_b_A[c2_jrow + 7 * c2_j].re - c2_s.im *
              c2_b_A[c2_jrow + 7 * c2_j].im;
            c2_cfrom1 = c2_s.re * c2_b_A[c2_jrow + 7 * c2_j].im + c2_s.im *
              c2_b_A[c2_jrow + 7 * c2_j].re;
            c2_d_A = c2_b_A[c2_jrowm1 + 7 * c2_j];
            c2_e_A = c2_b_A[c2_jrowm1 + 7 * c2_j];
            c2_b_A[c2_jrow + 7 * c2_j].re = c2_a * c2_b_A[c2_jrow + 7 * c2_j].re
              - (c2_s.re * c2_b_A[c2_jrowm1 + 7 * c2_j].re + c2_s.im *
                 c2_b_A[c2_jrowm1 + 7 * c2_j].im);
            c2_b_A[c2_jrow + 7 * c2_j].im = c2_a * c2_b_A[c2_jrow + 7 * c2_j].im
              - (c2_s.re * c2_d_A.im - c2_s.im * c2_e_A.re);
            c2_b_A[c2_jrowm1 + 7 * c2_j].re = c2_atmp.re + c2_b;
            c2_b_A[c2_jrowm1 + 7 * c2_j].im = c2_atmp.im + c2_cfrom1;
          }

          c2_s.re = -c2_s.re;
          c2_s.im = -c2_s.im;
          for (c2_i = c2_ilo - 1; c2_i + 1 <= c2_ihi; c2_i++) {
            c2_atmp.re = c2_a * c2_b_A[c2_i + 7 * c2_jrow].re;
            c2_atmp.im = c2_a * c2_b_A[c2_i + 7 * c2_jrow].im;
            c2_b = c2_s.re * c2_b_A[c2_i + 7 * c2_jrowm1].re - c2_s.im *
              c2_b_A[c2_i + 7 * c2_jrowm1].im;
            c2_cfrom1 = c2_s.re * c2_b_A[c2_i + 7 * c2_jrowm1].im + c2_s.im *
              c2_b_A[c2_i + 7 * c2_jrowm1].re;
            c2_d_A = c2_b_A[c2_i + 7 * c2_jrow];
            c2_e_A = c2_b_A[c2_i + 7 * c2_jrow];
            c2_b_A[c2_i + 7 * c2_jrowm1].re = c2_a * c2_b_A[c2_i + 7 * c2_jrowm1]
              .re - (c2_s.re * c2_b_A[c2_i + 7 * c2_jrow].re + c2_s.im *
                     c2_b_A[c2_i + 7 * c2_jrow].im);
            c2_b_A[c2_i + 7 * c2_jrowm1].im = c2_a * c2_b_A[c2_i + 7 * c2_jrowm1]
              .im - (c2_s.re * c2_d_A.im - c2_s.im * c2_e_A.re);
            c2_b_A[c2_i + 7 * c2_jrow].re = c2_atmp.re + c2_b;
            c2_b_A[c2_i + 7 * c2_jrow].im = c2_atmp.im + c2_cfrom1;
          }

          c2_jrow = c2_jrowm1;
        }

        c2_ii = c2_nzcount;
      }
    }

    c2_eml_matlab_zhgeqz(chartInstance, c2_b_A, c2_ilo, c2_ihi, &c2_info, c2_V,
                         c2_beta1);
    if ((!(c2_info != 0.0)) && c2_ilascl) {
      c2_notdone = TRUE;
      while (c2_notdone) {
        c2_cfrom1 = c2_anrmto * 2.0041683600089728E-292;
        c2_cto1 = c2_anrm / 4.9896007738368E+291;
        if ((c2_cfrom1 > c2_anrm) && (c2_anrm != 0.0)) {
          c2_mul = 2.0041683600089728E-292;
          c2_anrmto = c2_cfrom1;
        } else if (c2_cto1 > c2_anrmto) {
          c2_mul = 4.9896007738368E+291;
          c2_anrm = c2_cto1;
        } else {
          c2_mul = c2_anrm / c2_anrmto;
          c2_notdone = FALSE;
        }

        for (c2_ii = 0; c2_ii < 7; c2_ii++) {
          c2_V[c2_ii].re *= c2_mul;
          c2_V[c2_ii].im *= c2_mul;
        }
      }
    }
  }

  for (c2_ii = 0; c2_ii < 7; c2_ii++) {
    c2_cto1 = c2_V[c2_ii].re;
    c2_mul = c2_V[c2_ii].im;
    if (c2_beta1[c2_ii].im == 0.0) {
      if (c2_V[c2_ii].im == 0.0) {
        c2_V[c2_ii].re /= c2_beta1[c2_ii].re;
        c2_V[c2_ii].im = 0.0;
      } else if (c2_V[c2_ii].re == 0.0) {
        c2_V[c2_ii].re = 0.0;
        c2_V[c2_ii].im = c2_mul / c2_beta1[c2_ii].re;
      } else {
        c2_V[c2_ii].re /= c2_beta1[c2_ii].re;
        c2_V[c2_ii].im = c2_mul / c2_beta1[c2_ii].re;
      }
    } else if (c2_beta1[c2_ii].re == 0.0) {
      if (c2_V[c2_ii].re == 0.0) {
        c2_V[c2_ii].re = c2_V[c2_ii].im / c2_beta1[c2_ii].im;
        c2_V[c2_ii].im = 0.0;
      } else if (c2_V[c2_ii].im == 0.0) {
        c2_V[c2_ii].re = 0.0;
        c2_V[c2_ii].im = -(c2_cto1 / c2_beta1[c2_ii].im);
      } else {
        c2_V[c2_ii].re = c2_V[c2_ii].im / c2_beta1[c2_ii].im;
        c2_V[c2_ii].im = -(c2_cto1 / c2_beta1[c2_ii].im);
      }
    } else {
      c2_cfrom1 = muDoubleScalarAbs(c2_beta1[c2_ii].re);
      c2_a = muDoubleScalarAbs(c2_beta1[c2_ii].im);
      if (c2_cfrom1 > c2_a) {
        c2_a = c2_beta1[c2_ii].im / c2_beta1[c2_ii].re;
        c2_b = c2_beta1[c2_ii].re + c2_a * c2_beta1[c2_ii].im;
        c2_V[c2_ii].re = (c2_V[c2_ii].re + c2_a * c2_V[c2_ii].im) / c2_b;
        c2_V[c2_ii].im = (c2_mul - c2_a * c2_cto1) / c2_b;
      } else if (c2_a == c2_cfrom1) {
        c2_a = c2_beta1[c2_ii].re > 0.0 ? 0.5 : -0.5;
        c2_b = c2_beta1[c2_ii].im > 0.0 ? 0.5 : -0.5;
        c2_V[c2_ii].re = (c2_V[c2_ii].re * c2_a + c2_V[c2_ii].im * c2_b) /
          c2_cfrom1;
        c2_V[c2_ii].im = (c2_mul * c2_a - c2_cto1 * c2_b) / c2_cfrom1;
      } else {
        c2_a = c2_beta1[c2_ii].re / c2_beta1[c2_ii].im;
        c2_b = c2_beta1[c2_ii].im + c2_a * c2_beta1[c2_ii].re;
        c2_V[c2_ii].re = (c2_a * c2_V[c2_ii].re + c2_V[c2_ii].im) / c2_b;
        c2_V[c2_ii].im = (c2_a * c2_mul - c2_cto1) / c2_b;
      }
    }
  }

  if (c2_info < 0.0) {
    c2_b_eml_warning(chartInstance);
  } else {
    if (c2_info > 0.0) {
      c2_c_eml_warning(chartInstance);
    }
  }
}

static void c2_eml_matlab_zlartg(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_f, creal_T c2_g, real_T *c2_cs, creal_T *c2_sn,
  creal_T *c2_r)
{
  real_T c2_scale;
  real_T c2_dr;
  real_T c2_f2;
  real_T c2_fs_re;
  real_T c2_fs_im;
  real_T c2_gs_re;
  real_T c2_gs_im;
  int32_T c2_count;
  int32_T c2_rescaledir;
  boolean_T c2_guard1 = FALSE;
  real_T c2_b;
  real_T c2_g2s;
  c2_scale = muDoubleScalarAbs(c2_f.re);
  c2_dr = muDoubleScalarAbs(c2_f.im);
  if (c2_dr > c2_scale) {
    c2_scale = c2_dr;
  }

  c2_dr = muDoubleScalarAbs(c2_g.re);
  c2_f2 = muDoubleScalarAbs(c2_g.im);
  if (c2_f2 > c2_dr) {
    c2_dr = c2_f2;
  }

  if (c2_dr > c2_scale) {
    c2_scale = c2_dr;
  }

  c2_fs_re = c2_f.re;
  c2_fs_im = c2_f.im;
  c2_gs_re = c2_g.re;
  c2_gs_im = c2_g.im;
  c2_count = 0;
  c2_rescaledir = 0;
  c2_guard1 = FALSE;
  if (c2_scale >= 7.4428285367870146E+137) {
    do {
      c2_count++;
      c2_fs_re *= 1.3435752215134178E-138;
      c2_fs_im *= 1.3435752215134178E-138;
      c2_gs_re *= 1.3435752215134178E-138;
      c2_gs_im *= 1.3435752215134178E-138;
      c2_scale *= 1.3435752215134178E-138;
    } while (!(c2_scale < 7.4428285367870146E+137));

    c2_rescaledir = 1;
    c2_guard1 = TRUE;
  } else if (c2_scale <= 1.3435752215134178E-138) {
    if ((c2_g.re == 0.0) && (c2_g.im == 0.0)) {
      *c2_cs = 1.0;
      c2_sn->re = 0.0;
      c2_sn->im = 0.0;
      *c2_r = c2_f;
    } else {
      do {
        c2_count++;
        c2_fs_re *= 7.4428285367870146E+137;
        c2_fs_im *= 7.4428285367870146E+137;
        c2_gs_re *= 7.4428285367870146E+137;
        c2_gs_im *= 7.4428285367870146E+137;
        c2_scale *= 7.4428285367870146E+137;
      } while (!(c2_scale > 1.3435752215134178E-138));

      c2_rescaledir = -1;
      c2_guard1 = TRUE;
    }
  } else {
    c2_guard1 = TRUE;
  }

  if (c2_guard1 == TRUE) {
    c2_f2 = c2_fs_re * c2_fs_re + c2_fs_im * c2_fs_im;
    c2_scale = c2_gs_re * c2_gs_re + c2_gs_im * c2_gs_im;
    c2_dr = c2_scale;
    if (1.0 > c2_scale) {
      c2_dr = 1.0;
    }

    if (c2_f2 <= c2_dr * 2.0041683600089728E-292) {
      if ((c2_f.re == 0.0) && (c2_f.im == 0.0)) {
        *c2_cs = 0.0;
        c2_f2 = muDoubleScalarAbs(c2_g.re);
        c2_b = muDoubleScalarAbs(c2_g.im);
        if (c2_f2 < c2_b) {
          c2_f2 /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
        } else if (c2_f2 > c2_b) {
          c2_b /= c2_f2;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_f2 * 1.4142135623730951;
          }
        }

        c2_r->re = c2_b;
        c2_r->im = 0.0;
        c2_f2 = muDoubleScalarAbs(c2_gs_re);
        c2_b = muDoubleScalarAbs(c2_gs_im);
        if (c2_f2 < c2_b) {
          c2_f2 /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
        } else if (c2_f2 > c2_b) {
          c2_b /= c2_f2;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_f2 * 1.4142135623730951;
          }
        }

        c2_sn->re = c2_gs_re / c2_b;
        c2_sn->im = -c2_gs_im / c2_b;
      } else {
        c2_f2 = muDoubleScalarAbs(c2_fs_re);
        c2_b = muDoubleScalarAbs(c2_fs_im);
        if (c2_f2 < c2_b) {
          c2_f2 /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
        } else if (c2_f2 > c2_b) {
          c2_b /= c2_f2;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_f2 * 1.4142135623730951;
          }
        }

        if (c2_scale < 0.0) {
          c2_eml_error(chartInstance);
        }

        c2_g2s = muDoubleScalarSqrt(c2_scale);
        *c2_cs = c2_b / c2_g2s;
        c2_dr = muDoubleScalarAbs(c2_f.re);
        c2_f2 = muDoubleScalarAbs(c2_f.im);
        if (c2_f2 > c2_dr) {
          c2_dr = c2_f2;
        }

        if (c2_dr > 1.0) {
          c2_f2 = muDoubleScalarAbs(c2_f.re);
          c2_b = muDoubleScalarAbs(c2_f.im);
          if (c2_f2 < c2_b) {
            c2_f2 /= c2_b;
            c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
          } else if (c2_f2 > c2_b) {
            c2_b /= c2_f2;
            c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
          } else {
            if (!muDoubleScalarIsNaN(c2_b)) {
              c2_b = c2_f2 * 1.4142135623730951;
            }
          }

          c2_fs_re = c2_f.re / c2_b;
          c2_fs_im = c2_f.im / c2_b;
        } else {
          c2_dr = 7.4428285367870146E+137 * c2_f.re;
          c2_scale = 7.4428285367870146E+137 * c2_f.im;
          c2_f2 = muDoubleScalarAbs(c2_dr);
          c2_b = muDoubleScalarAbs(c2_scale);
          if (c2_f2 < c2_b) {
            c2_f2 /= c2_b;
            c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
          } else if (c2_f2 > c2_b) {
            c2_b /= c2_f2;
            c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
          } else {
            if (!muDoubleScalarIsNaN(c2_b)) {
              c2_b = c2_f2 * 1.4142135623730951;
            }
          }

          c2_fs_re = c2_dr / c2_b;
          c2_fs_im = c2_scale / c2_b;
        }

        c2_gs_re /= c2_g2s;
        c2_gs_im = -c2_gs_im / c2_g2s;
        c2_sn->re = c2_fs_re * c2_gs_re - c2_fs_im * c2_gs_im;
        c2_sn->im = c2_fs_re * c2_gs_im + c2_fs_im * c2_gs_re;
        c2_r->re = *c2_cs * c2_f.re + (c2_sn->re * c2_g.re - c2_sn->im * c2_g.im);
        c2_r->im = *c2_cs * c2_f.im + (c2_sn->re * c2_g.im + c2_sn->im * c2_g.re);
      }
    } else {
      c2_dr = 1.0 + c2_scale / c2_f2;
      if (c2_dr < 0.0) {
        c2_eml_error(chartInstance);
      }

      c2_b = muDoubleScalarSqrt(c2_dr);
      c2_r->re = c2_b * c2_fs_re;
      c2_r->im = c2_b * c2_fs_im;
      *c2_cs = 1.0 / c2_b;
      c2_b = c2_f2 + c2_scale;
      c2_f2 = c2_r->re / c2_b;
      c2_dr = c2_r->im / c2_b;
      c2_sn->re = c2_f2 * c2_gs_re - c2_dr * -c2_gs_im;
      c2_sn->im = c2_f2 * -c2_gs_im + c2_dr * c2_gs_re;
      if (c2_rescaledir > 0) {
        for (c2_rescaledir = 1; c2_rescaledir <= c2_count; c2_rescaledir++) {
          c2_r->re *= 7.4428285367870146E+137;
          c2_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c2_rescaledir < 0) {
          for (c2_rescaledir = 1; c2_rescaledir <= c2_count; c2_rescaledir++) {
            c2_r->re *= 1.3435752215134178E-138;
            c2_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c2_eml_matlab_zhgeqz(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_A[49], int32_T c2_ilo, int32_T c2_ihi, real_T
  *c2_info, creal_T c2_alpha1[7], creal_T c2_beta1[7])
{
  static creal_T c2_dc1 = { 0.0, 0.0 };

  int32_T c2_jm1;
  creal_T c2_b_A[49];
  real_T c2_eshift_re;
  real_T c2_eshift_im;
  static creal_T c2_dc2 = { 0.0, 0.0 };

  creal_T c2_ctemp;
  real_T c2_rho_re;
  real_T c2_rho_im;
  real_T c2_anorm;
  real_T c2_a;
  real_T c2_atol;
  real_T c2_ascale;
  boolean_T c2_failed;
  int32_T c2_j;
  boolean_T c2_guard1 = FALSE;
  boolean_T c2_guard2 = FALSE;
  int32_T c2_ifirst;
  int32_T c2_istart;
  int32_T c2_ilast;
  int32_T c2_ilastm1;
  int32_T c2_ifrstm;
  int32_T c2_iiter;
  int32_T c2_maxit;
  boolean_T c2_goto60;
  boolean_T c2_goto70;
  boolean_T c2_goto90;
  int32_T c2_jiter;
  int32_T c2_exitg1;
  boolean_T c2_exitg3;
  boolean_T c2_ilazro;
  boolean_T c2_b_guard1 = FALSE;
  creal_T c2_c_A;
  creal_T c2_t1;
  creal_T c2_d;
  creal_T c2_sigma1;
  real_T c2_sigma2_re;
  real_T c2_sigma2_im;
  real_T c2_b;
  int32_T c2_jp1;
  boolean_T c2_exitg2;
  creal_T c2_d_A;
  int32_T c2_i;
  c2_dc1.re = rtNaN;
  for (c2_jm1 = 0; c2_jm1 < 49; c2_jm1++) {
    c2_b_A[c2_jm1] = c2_A[c2_jm1];
  }

  for (c2_jm1 = 0; c2_jm1 < 7; c2_jm1++) {
    c2_alpha1[c2_jm1].re = 0.0;
    c2_alpha1[c2_jm1].im = 0.0;
    c2_beta1[c2_jm1].re = 1.0;
    c2_beta1[c2_jm1].im = 0.0;
  }

  c2_eshift_re = 0.0;
  c2_eshift_im = 0.0;
  c2_ctemp = c2_dc2;
  c2_rho_re = 0.0;
  c2_rho_im = 0.0;
  c2_anorm = c2_eml_matlab_zlanhs(chartInstance, c2_A, c2_ilo, c2_ihi);
  c2_a = 2.2204460492503131E-16 * c2_anorm;
  c2_atol = 2.2250738585072014E-308;
  if (c2_a > 2.2250738585072014E-308) {
    c2_atol = c2_a;
  }

  c2_a = 2.2250738585072014E-308;
  if (c2_anorm > 2.2250738585072014E-308) {
    c2_a = c2_anorm;
  }

  c2_ascale = 1.0 / c2_a;
  c2_failed = TRUE;
  for (c2_j = c2_ihi + 1; c2_j < 8; c2_j++) {
    c2_alpha1[sf_mex_lw_bounds_check(c2_j, 1, 7) - 1] = c2_A
      [(sf_mex_lw_bounds_check(c2_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c2_j, 1,
          7) - 1)) - 1];
  }

  c2_guard1 = FALSE;
  c2_guard2 = FALSE;
  if (c2_ihi >= c2_ilo) {
    c2_ifirst = c2_ilo;
    c2_istart = c2_ilo;
    c2_ilast = c2_ihi;
    c2_ilastm1 = c2_ihi - 2;
    c2_ifrstm = c2_ilo;
    c2_iiter = 0;
    c2_maxit = 30 * ((c2_ihi - c2_ilo) + 1);
    c2_goto60 = FALSE;
    c2_goto70 = FALSE;
    c2_goto90 = FALSE;
    c2_jiter = 1;
    do {
      c2_exitg1 = 0U;
      if (c2_jiter <= c2_maxit) {
        if (c2_ilast == c2_ilo) {
          c2_goto60 = TRUE;
        } else {
          sf_mex_lw_bounds_check(c2_ilastm1 + 1, 1, 7);
          sf_mex_lw_bounds_check(c2_ilast, 1, 7);
          if (muDoubleScalarAbs(c2_b_A[(c2_ilast + 7 * c2_ilastm1) - 1].re) +
              muDoubleScalarAbs(c2_b_A[(c2_ilast + 7 * c2_ilastm1) - 1].im) <=
              c2_atol) {
            c2_b_A[(c2_ilast + 7 * c2_ilastm1) - 1] = c2_dc2;
            c2_goto60 = TRUE;
          } else {
            c2_j = c2_ilastm1;
            c2_exitg3 = 0U;
            while ((c2_exitg3 == 0U) && (c2_j + 1 >= c2_ilo)) {
              c2_jm1 = c2_j - 1;
              if (c2_j + 1 == c2_ilo) {
                c2_ilazro = TRUE;
              } else {
                sf_mex_lw_bounds_check(c2_jm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c2_j + 1, 1, 7);
                if (muDoubleScalarAbs(c2_b_A[c2_j + 7 * c2_jm1].re) +
                    muDoubleScalarAbs(c2_b_A[c2_j + 7 * c2_jm1].im) <= c2_atol)
                {
                  c2_b_A[c2_j + 7 * c2_jm1] = c2_dc2;
                  c2_ilazro = TRUE;
                } else {
                  c2_ilazro = FALSE;
                }
              }

              if (c2_ilazro) {
                c2_ifirst = c2_j + 1;
                c2_goto70 = TRUE;
                c2_exitg3 = 1U;
              } else {
                c2_j = c2_jm1;
              }
            }
          }
        }

        if (c2_goto60 || c2_goto70) {
          c2_ilazro = TRUE;
        } else {
          c2_ilazro = FALSE;
        }

        if (!c2_ilazro) {
          for (c2_jm1 = 0; c2_jm1 < 7; c2_jm1++) {
            c2_alpha1[c2_jm1] = c2_dc1;
            c2_beta1[c2_jm1] = c2_dc1;
          }

          *c2_info = -1.0;
          c2_exitg1 = 1U;
        } else {
          c2_b_guard1 = FALSE;
          if (c2_goto60) {
            c2_goto60 = FALSE;
            c2_alpha1[sf_mex_lw_bounds_check(c2_ilast, 1, 7) - 1] = c2_b_A
              [(sf_mex_lw_bounds_check(c2_ilast, 1, 7) + 7 *
                (sf_mex_lw_bounds_check(c2_ilast, 1, 7) - 1)) - 1];
            c2_ilast = c2_ilastm1 + 1;
            c2_ilastm1--;
            if (c2_ilast < c2_ilo) {
              c2_failed = FALSE;
              c2_guard2 = TRUE;
              c2_exitg1 = 1U;
            } else {
              c2_iiter = 0;
              c2_eshift_re = 0.0;
              c2_eshift_im = 0.0;
              if (c2_ifrstm > c2_ilast) {
                c2_ifrstm = c2_ilo;
              }

              c2_b_guard1 = TRUE;
            }
          } else {
            if (c2_goto70) {
              c2_goto70 = FALSE;
              c2_iiter++;
              c2_ifrstm = c2_ifirst;
              if (c2_iiter - c2_div_s32_floor(chartInstance, c2_iiter, 10) * 10
                  != 0) {
                sf_mex_lw_bounds_check(c2_ilastm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c2_ilast, 1, 7);
                c2_c_A.re = -(c2_b_A[(c2_ilast + 7 * (c2_ilast - 1)) - 1].re -
                              c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].re);
                c2_c_A.im = -(c2_b_A[(c2_ilast + 7 * (c2_ilast - 1)) - 1].im -
                              c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].im);
                c2_t1 = c2_eml_div(chartInstance, c2_c_A, 2.0);
                c2_d.re = (c2_t1.re * c2_t1.re - c2_t1.im * c2_t1.im) +
                  (c2_b_A[c2_ilastm1 + 7 * (c2_ilast - 1)].re * c2_b_A[(c2_ilast
                    + 7 * c2_ilastm1) - 1].re - c2_b_A[c2_ilastm1 + 7 *
                   (c2_ilast - 1)].im * c2_b_A[(c2_ilast + 7 * c2_ilastm1) - 1].
                   im);
                c2_d.im = (c2_t1.re * c2_t1.im + c2_t1.im * c2_t1.re) +
                  (c2_b_A[c2_ilastm1 + 7 * (c2_ilast - 1)].re * c2_b_A[(c2_ilast
                    + 7 * c2_ilastm1) - 1].im + c2_b_A[c2_ilastm1 + 7 *
                   (c2_ilast - 1)].im * c2_b_A[(c2_ilast + 7 * c2_ilastm1) - 1].
                   re);
                c2_sqrt(chartInstance, &c2_d);
                c2_sigma1.re = c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].re -
                  (c2_t1.re - c2_d.re);
                c2_sigma1.im = c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].im -
                  (c2_t1.im - c2_d.im);
                c2_sigma2_re = c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].re -
                  (c2_t1.re + c2_d.re);
                c2_sigma2_im = c2_b_A[c2_ilastm1 + 7 * c2_ilastm1].im -
                  (c2_t1.im + c2_d.im);
                c2_a = muDoubleScalarAbs(c2_sigma1.re - c2_b_A[(c2_ilast + 7 *
                  (c2_ilast - 1)) - 1].re);
                c2_anorm = muDoubleScalarAbs(c2_sigma1.im - c2_b_A[(c2_ilast + 7
                  * (c2_ilast - 1)) - 1].im);
                if (c2_a < c2_anorm) {
                  c2_a /= c2_anorm;
                  c2_anorm *= muDoubleScalarSqrt(c2_a * c2_a + 1.0);
                } else if (c2_a > c2_anorm) {
                  c2_anorm /= c2_a;
                  c2_anorm = muDoubleScalarSqrt(c2_anorm * c2_anorm + 1.0) *
                    c2_a;
                } else {
                  if (!muDoubleScalarIsNaN(c2_anorm)) {
                    c2_anorm = c2_a * 1.4142135623730951;
                  }
                }

                c2_a = muDoubleScalarAbs(c2_sigma2_re - c2_b_A[(c2_ilast + 7 *
                  (c2_ilast - 1)) - 1].re);
                c2_b = muDoubleScalarAbs(c2_sigma2_im - c2_b_A[(c2_ilast + 7 *
                  (c2_ilast - 1)) - 1].im);
                if (c2_a < c2_b) {
                  c2_a /= c2_b;
                  c2_b *= muDoubleScalarSqrt(c2_a * c2_a + 1.0);
                } else if (c2_a > c2_b) {
                  c2_b /= c2_a;
                  c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_a;
                } else {
                  if (!muDoubleScalarIsNaN(c2_b)) {
                    c2_b = c2_a * 1.4142135623730951;
                  }
                }

                if (c2_anorm <= c2_b) {
                  c2_sigma2_re = c2_sigma1.re;
                  c2_sigma2_im = c2_sigma1.im;
                  c2_rho_re = c2_t1.re - c2_d.re;
                  c2_rho_im = c2_t1.im - c2_d.im;
                } else {
                  c2_rho_re = c2_t1.re + c2_d.re;
                  c2_rho_im = c2_t1.im + c2_d.im;
                }
              } else {
                c2_eshift_re += c2_b_A[(sf_mex_lw_bounds_check(c2_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c2_ilastm1 + 1, 1, 7) - 1)) - 1].
                  re;
                c2_eshift_im += c2_b_A[(sf_mex_lw_bounds_check(c2_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c2_ilastm1 + 1, 1, 7) - 1)) - 1].
                  im;
                c2_sigma2_re = c2_eshift_re;
                c2_sigma2_im = c2_eshift_im;
              }

              c2_j = c2_ilastm1;
              c2_jp1 = c2_ilastm1 + 1;
              c2_exitg2 = 0U;
              while ((c2_exitg2 == 0U) && (c2_j + 1 > c2_ifirst)) {
                c2_jm1 = c2_j - 1;
                c2_istart = c2_j + 1;
                c2_ctemp.re = c2_b_A[c2_j + 7 * c2_j].re - c2_sigma2_re;
                c2_ctemp.im = c2_b_A[c2_j + 7 * c2_j].im - c2_sigma2_im;
                c2_anorm = c2_ascale * (muDoubleScalarAbs(c2_ctemp.re) +
                  muDoubleScalarAbs(c2_ctemp.im));
                sf_mex_lw_bounds_check(c2_jp1 + 1, 1, 7);
                c2_a = c2_ascale * (muDoubleScalarAbs(c2_b_A[c2_jp1 + 7 * c2_j].
                  re) + muDoubleScalarAbs(c2_b_A[c2_jp1 + 7 * c2_j].im));
                c2_b = c2_anorm;
                if (c2_a > c2_anorm) {
                  c2_b = c2_a;
                }

                if ((c2_b < 1.0) && (c2_b != 0.0)) {
                  c2_anorm /= c2_b;
                  c2_a /= c2_b;
                }

                sf_mex_lw_bounds_check(c2_jm1 + 1, 1, 7);
                if ((muDoubleScalarAbs(c2_b_A[c2_j + 7 * c2_jm1].re) +
                     muDoubleScalarAbs(c2_b_A[c2_j + 7 * c2_jm1].im)) * c2_a <=
                    c2_anorm * c2_atol) {
                  c2_goto90 = TRUE;
                  c2_exitg2 = 1U;
                } else {
                  c2_jp1 = c2_j;
                  c2_j = c2_jm1;
                }
              }

              if (!c2_goto90) {
                c2_istart = c2_ifirst;
                if (c2_ifirst == c2_ilastm1 + 1) {
                  c2_ctemp.re = c2_rho_re;
                  c2_ctemp.im = c2_rho_im;
                } else {
                  c2_ctemp.re = c2_b_A[(sf_mex_lw_bounds_check(c2_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c2_ifirst,
                    1, 7) - 1)) - 1].re - c2_sigma2_re;
                  c2_ctemp.im = c2_b_A[(sf_mex_lw_bounds_check(c2_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c2_ifirst,
                    1, 7) - 1)) - 1].im - c2_sigma2_im;
                }

                c2_goto90 = TRUE;
              }
            }

            if (c2_goto90) {
              c2_goto90 = FALSE;
              sf_mex_lw_bounds_check(c2_istart, 1, 7);
              sf_mex_lw_bounds_check(c2_istart + 1, 1, 7);
              c2_b_eml_matlab_zlartg(chartInstance, c2_ctemp, c2_b_A[c2_istart +
                7 * (c2_istart - 1)], &c2_a, &c2_sigma1);
              c2_j = c2_istart - 1;
              c2_jm1 = c2_istart - 2;
              while (c2_j + 1 < c2_ilast) {
                c2_jp1 = c2_j + 1;
                if (c2_j + 1 > c2_istart) {
                  c2_c_A = c2_b_A[(sf_mex_lw_bounds_check(c2_j + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c2_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c2_d_A = c2_b_A[(sf_mex_lw_bounds_check(c2_jp1 + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c2_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c2_eml_matlab_zlartg(chartInstance, c2_c_A, c2_d_A, &c2_a,
                                       &c2_sigma1, &c2_t1);
                  c2_b_A[c2_j + 7 * c2_jm1] = c2_t1;
                  c2_b_A[c2_jp1 + 7 * c2_jm1] = c2_dc2;
                }

                for (c2_jm1 = c2_j; c2_jm1 + 1 <= c2_ilast; c2_jm1++) {
                  sf_mex_lw_bounds_check(c2_jm1 + 1, 1, 7);
                  c2_t1.re = c2_a * c2_b_A[c2_j + 7 * c2_jm1].re;
                  c2_t1.im = c2_a * c2_b_A[c2_j + 7 * c2_jm1].im;
                  sf_mex_lw_bounds_check(c2_jp1 + 1, 1, 7);
                  c2_d.re = c2_sigma1.re * c2_b_A[c2_jp1 + 7 * c2_jm1].re -
                    c2_sigma1.im * c2_b_A[c2_jp1 + 7 * c2_jm1].im;
                  c2_d.im = c2_sigma1.re * c2_b_A[c2_jp1 + 7 * c2_jm1].im +
                    c2_sigma1.im * c2_b_A[c2_jp1 + 7 * c2_jm1].re;
                  c2_c_A = c2_b_A[c2_j + 7 * c2_jm1];
                  c2_d_A = c2_b_A[c2_j + 7 * c2_jm1];
                  c2_b_A[c2_jp1 + 7 * c2_jm1].re = c2_a * c2_b_A[c2_jp1 + 7 *
                    c2_jm1].re - (c2_sigma1.re * c2_b_A[c2_j + 7 * c2_jm1].re +
                                  c2_sigma1.im * c2_b_A[c2_j + 7 * c2_jm1].im);
                  c2_b_A[c2_jp1 + 7 * c2_jm1].im = c2_a * c2_b_A[c2_jp1 + 7 *
                    c2_jm1].im - (c2_sigma1.re * c2_c_A.im - c2_sigma1.im *
                                  c2_d_A.re);
                  c2_b_A[c2_j + 7 * c2_jm1].re = c2_t1.re + c2_d.re;
                  c2_b_A[c2_j + 7 * c2_jm1].im = c2_t1.im + c2_d.im;
                }

                c2_sigma1.re = -c2_sigma1.re;
                c2_sigma1.im = -c2_sigma1.im;
                c2_jm1 = c2_jp1 + 2;
                if (c2_ilast < c2_jm1) {
                  c2_jm1 = c2_ilast;
                }

                for (c2_i = c2_ifrstm - 1; c2_i + 1 <= c2_jm1; c2_i++) {
                  sf_mex_lw_bounds_check(c2_jp1 + 1, 1, 7);
                  sf_mex_lw_bounds_check(c2_i + 1, 1, 7);
                  c2_t1.re = c2_a * c2_b_A[c2_i + 7 * c2_jp1].re;
                  c2_t1.im = c2_a * c2_b_A[c2_i + 7 * c2_jp1].im;
                  c2_d.re = c2_sigma1.re * c2_b_A[c2_i + 7 * c2_j].re -
                    c2_sigma1.im * c2_b_A[c2_i + 7 * c2_j].im;
                  c2_d.im = c2_sigma1.re * c2_b_A[c2_i + 7 * c2_j].im +
                    c2_sigma1.im * c2_b_A[c2_i + 7 * c2_j].re;
                  c2_c_A = c2_b_A[c2_i + 7 * c2_jp1];
                  c2_d_A = c2_b_A[c2_i + 7 * c2_jp1];
                  c2_b_A[c2_i + 7 * c2_j].re = c2_a * c2_b_A[c2_i + 7 * c2_j].re
                    - (c2_sigma1.re * c2_b_A[c2_i + 7 * c2_jp1].re +
                       c2_sigma1.im * c2_b_A[c2_i + 7 * c2_jp1].im);
                  c2_b_A[c2_i + 7 * c2_j].im = c2_a * c2_b_A[c2_i + 7 * c2_j].im
                    - (c2_sigma1.re * c2_c_A.im - c2_sigma1.im * c2_d_A.re);
                  c2_b_A[c2_i + 7 * c2_jp1].re = c2_t1.re + c2_d.re;
                  c2_b_A[c2_i + 7 * c2_jp1].im = c2_t1.im + c2_d.im;
                }

                c2_jm1 = c2_j;
                c2_j = c2_jp1;
              }
            }

            c2_b_guard1 = TRUE;
          }

          if (c2_b_guard1 == TRUE) {
            c2_jiter++;
          }
        }
      } else {
        c2_guard2 = TRUE;
        c2_exitg1 = 1U;
      }
    } while (c2_exitg1 == 0U);
  } else {
    c2_guard1 = TRUE;
  }

  if (c2_guard2 == TRUE) {
    if (c2_failed) {
      *c2_info = (real_T)c2_ilast;
      for (c2_jm1 = 1; c2_jm1 <= c2_ilast; c2_jm1++) {
        c2_alpha1[sf_mex_lw_bounds_check(c2_jm1, 1, 7) - 1] = c2_dc1;
        c2_beta1[c2_jm1 - 1] = c2_dc1;
      }
    } else {
      c2_guard1 = TRUE;
    }
  }

  if (c2_guard1 == TRUE) {
    for (c2_j = 1; c2_j <= c2_ilo - 1; c2_j++) {
      c2_alpha1[sf_mex_lw_bounds_check(c2_j, 1, 7) - 1] = c2_b_A
        [(sf_mex_lw_bounds_check(c2_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c2_j,
            1, 7) - 1)) - 1];
    }

    *c2_info = 0.0;
  }
}

static real_T c2_eml_matlab_zlanhs(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_A[49], int32_T c2_ilo, int32_T c2_ihi)
{
  real_T c2_f;
  real_T c2_scale;
  real_T c2_sumsq;
  boolean_T c2_firstNonZero;
  int32_T c2_j;
  int32_T c2_c;
  int32_T c2_i;
  real_T c2_temp1;
  real_T c2_temp2;
  c2_f = 0.0;
  if (!(c2_ilo > c2_ihi)) {
    c2_scale = 0.0;
    c2_sumsq = 0.0;
    c2_firstNonZero = TRUE;
    for (c2_j = c2_ilo; c2_j <= c2_ihi; c2_j++) {
      c2_c = c2_j + 1;
      if (c2_ihi < c2_c) {
        c2_c = c2_ihi;
      }

      for (c2_i = c2_ilo; c2_i <= c2_c; c2_i++) {
        sf_mex_lw_bounds_check(c2_i, 1, 7);
        sf_mex_lw_bounds_check(c2_j, 1, 7);
        if (c2_A[(c2_i + 7 * (c2_j - 1)) - 1].re != 0.0) {
          c2_temp1 = muDoubleScalarAbs(c2_A[(c2_i + 7 * (c2_j - 1)) - 1].re);
          if (c2_firstNonZero) {
            c2_sumsq = 1.0;
            c2_scale = c2_temp1;
            c2_firstNonZero = FALSE;
          } else if (c2_scale < c2_temp1) {
            c2_temp2 = c2_scale / c2_temp1;
            c2_sumsq = 1.0 + c2_sumsq * c2_temp2 * c2_temp2;
            c2_scale = c2_temp1;
          } else {
            c2_temp2 = c2_temp1 / c2_scale;
            c2_sumsq += c2_temp2 * c2_temp2;
          }
        }

        if (c2_A[(c2_i + 7 * (c2_j - 1)) - 1].im != 0.0) {
          c2_temp1 = muDoubleScalarAbs(c2_A[(c2_i + 7 * (c2_j - 1)) - 1].im);
          if (c2_firstNonZero) {
            c2_sumsq = 1.0;
            c2_scale = c2_temp1;
            c2_firstNonZero = FALSE;
          } else if (c2_scale < c2_temp1) {
            c2_temp2 = c2_scale / c2_temp1;
            c2_sumsq = 1.0 + c2_sumsq * c2_temp2 * c2_temp2;
            c2_scale = c2_temp1;
          } else {
            c2_temp2 = c2_temp1 / c2_scale;
            c2_sumsq += c2_temp2 * c2_temp2;
          }
        }
      }
    }

    if (c2_sumsq < 0.0) {
      c2_eml_error(chartInstance);
    }

    c2_f = c2_scale * muDoubleScalarSqrt(c2_sumsq);
  }

  return c2_f;
}

static creal_T c2_eml_div(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  creal_T c2_x, real_T c2_y)
{
  creal_T c2_z;
  if (c2_x.im == 0.0) {
    c2_z.re = c2_x.re / c2_y;
    c2_z.im = 0.0;
  } else if (c2_x.re == 0.0) {
    c2_z.re = 0.0;
    c2_z.im = c2_x.im / c2_y;
  } else {
    c2_z.re = c2_x.re / c2_y;
    c2_z.im = c2_x.im / c2_y;
  }

  return c2_z;
}

static void c2_b_eml_matlab_zlartg(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_f, creal_T c2_g, real_T *c2_cs, creal_T *c2_sn)
{
  real_T c2_scale;
  real_T c2_b;
  real_T c2_f2;
  real_T c2_fs_re;
  real_T c2_fs_im;
  real_T c2_gs_re;
  real_T c2_gs_im;
  boolean_T c2_guard1 = FALSE;
  real_T c2_b_b;
  real_T c2_g2s;
  c2_scale = muDoubleScalarAbs(c2_f.re);
  c2_b = muDoubleScalarAbs(c2_f.im);
  if (c2_b > c2_scale) {
    c2_scale = c2_b;
  }

  c2_b = muDoubleScalarAbs(c2_g.re);
  c2_f2 = muDoubleScalarAbs(c2_g.im);
  if (c2_f2 > c2_b) {
    c2_b = c2_f2;
  }

  if (c2_b > c2_scale) {
    c2_scale = c2_b;
  }

  c2_fs_re = c2_f.re;
  c2_fs_im = c2_f.im;
  c2_gs_re = c2_g.re;
  c2_gs_im = c2_g.im;
  c2_guard1 = FALSE;
  if (c2_scale >= 7.4428285367870146E+137) {
    do {
      c2_fs_re *= 1.3435752215134178E-138;
      c2_fs_im *= 1.3435752215134178E-138;
      c2_gs_re *= 1.3435752215134178E-138;
      c2_gs_im *= 1.3435752215134178E-138;
      c2_scale *= 1.3435752215134178E-138;
    } while (!(c2_scale < 7.4428285367870146E+137));

    c2_guard1 = TRUE;
  } else if (c2_scale <= 1.3435752215134178E-138) {
    if ((c2_g.re == 0.0) && (c2_g.im == 0.0)) {
      *c2_cs = 1.0;
      c2_sn->re = 0.0;
      c2_sn->im = 0.0;
    } else {
      do {
        c2_fs_re *= 7.4428285367870146E+137;
        c2_fs_im *= 7.4428285367870146E+137;
        c2_gs_re *= 7.4428285367870146E+137;
        c2_gs_im *= 7.4428285367870146E+137;
        c2_scale *= 7.4428285367870146E+137;
      } while (!(c2_scale > 1.3435752215134178E-138));

      c2_guard1 = TRUE;
    }
  } else {
    c2_guard1 = TRUE;
  }

  if (c2_guard1 == TRUE) {
    c2_f2 = c2_fs_re * c2_fs_re + c2_fs_im * c2_fs_im;
    c2_scale = c2_gs_re * c2_gs_re + c2_gs_im * c2_gs_im;
    c2_b = c2_scale;
    if (1.0 > c2_scale) {
      c2_b = 1.0;
    }

    if (c2_f2 <= c2_b * 2.0041683600089728E-292) {
      if ((c2_f.re == 0.0) && (c2_f.im == 0.0)) {
        *c2_cs = 0.0;
        c2_f2 = muDoubleScalarAbs(c2_gs_re);
        c2_b_b = muDoubleScalarAbs(c2_gs_im);
        if (c2_f2 < c2_b_b) {
          c2_f2 /= c2_b_b;
          c2_b_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
        } else if (c2_f2 > c2_b_b) {
          c2_b_b /= c2_f2;
          c2_b_b = muDoubleScalarSqrt(c2_b_b * c2_b_b + 1.0) * c2_f2;
        } else {
          if (!muDoubleScalarIsNaN(c2_b_b)) {
            c2_b_b = c2_f2 * 1.4142135623730951;
          }
        }

        c2_sn->re = c2_gs_re / c2_b_b;
        c2_sn->im = -c2_gs_im / c2_b_b;
      } else {
        c2_f2 = muDoubleScalarAbs(c2_fs_re);
        c2_b = muDoubleScalarAbs(c2_fs_im);
        if (c2_f2 < c2_b) {
          c2_f2 /= c2_b;
          c2_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
        } else if (c2_f2 > c2_b) {
          c2_b /= c2_f2;
          c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_f2;
        } else {
          if (!muDoubleScalarIsNaN(c2_b)) {
            c2_b = c2_f2 * 1.4142135623730951;
          }
        }

        if (c2_scale < 0.0) {
          c2_eml_error(chartInstance);
        }

        c2_g2s = muDoubleScalarSqrt(c2_scale);
        *c2_cs = c2_b / c2_g2s;
        c2_b = muDoubleScalarAbs(c2_f.re);
        c2_f2 = muDoubleScalarAbs(c2_f.im);
        if (c2_f2 > c2_b) {
          c2_b = c2_f2;
        }

        if (c2_b > 1.0) {
          c2_f2 = muDoubleScalarAbs(c2_f.re);
          c2_b_b = muDoubleScalarAbs(c2_f.im);
          if (c2_f2 < c2_b_b) {
            c2_f2 /= c2_b_b;
            c2_b_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
          } else if (c2_f2 > c2_b_b) {
            c2_b_b /= c2_f2;
            c2_b_b = muDoubleScalarSqrt(c2_b_b * c2_b_b + 1.0) * c2_f2;
          } else {
            if (!muDoubleScalarIsNaN(c2_b_b)) {
              c2_b_b = c2_f2 * 1.4142135623730951;
            }
          }

          c2_fs_re = c2_f.re / c2_b_b;
          c2_fs_im = c2_f.im / c2_b_b;
        } else {
          c2_b = 7.4428285367870146E+137 * c2_f.re;
          c2_scale = 7.4428285367870146E+137 * c2_f.im;
          c2_f2 = muDoubleScalarAbs(c2_b);
          c2_b_b = muDoubleScalarAbs(c2_scale);
          if (c2_f2 < c2_b_b) {
            c2_f2 /= c2_b_b;
            c2_b_b *= muDoubleScalarSqrt(c2_f2 * c2_f2 + 1.0);
          } else if (c2_f2 > c2_b_b) {
            c2_b_b /= c2_f2;
            c2_b_b = muDoubleScalarSqrt(c2_b_b * c2_b_b + 1.0) * c2_f2;
          } else {
            if (!muDoubleScalarIsNaN(c2_b_b)) {
              c2_b_b = c2_f2 * 1.4142135623730951;
            }
          }

          c2_fs_re = c2_b / c2_b_b;
          c2_fs_im = c2_scale / c2_b_b;
        }

        c2_gs_re /= c2_g2s;
        c2_gs_im = -c2_gs_im / c2_g2s;
        c2_sn->re = c2_fs_re * c2_gs_re - c2_fs_im * c2_gs_im;
        c2_sn->im = c2_fs_re * c2_gs_im + c2_fs_im * c2_gs_re;
      }
    } else {
      c2_b = 1.0 + c2_scale / c2_f2;
      if (c2_b < 0.0) {
        c2_eml_error(chartInstance);
      }

      c2_b = muDoubleScalarSqrt(c2_b);
      *c2_cs = 1.0 / c2_b;
      c2_b_b = c2_f2 + c2_scale;
      c2_fs_re = c2_b * c2_fs_re / c2_b_b;
      c2_fs_im = c2_b * c2_fs_im / c2_b_b;
      c2_sn->re = c2_fs_re * c2_gs_re - c2_fs_im * -c2_gs_im;
      c2_sn->im = c2_fs_re * -c2_gs_im + c2_fs_im * c2_gs_re;
    }
  }
}

static void c2_b_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  int32_T c2_i8;
  static char_T c2_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c2_u[26];
  const mxArray *c2_y = NULL;
  for (c2_i8 = 0; c2_i8 < 26; c2_i8++) {
    c2_u[c2_i8] = c2_varargin_1[c2_i8];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 26), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static void c2_c_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  int32_T c2_i9;
  static char_T c2_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c2_u[34];
  const mxArray *c2_y = NULL;
  for (c2_i9 = 0; c2_i9 < 34; c2_i9++) {
    c2_u[c2_i9] = c2_varargin_1[c2_i9];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 34), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static void c2_d_eye(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T c2_I[49])
{
  int32_T c2_i;
  int32_T c2_b_i;
  for (c2_i = 0; c2_i < 49; c2_i++) {
    c2_I[c2_i] = 0.0;
  }

  c2_i = 0;
  for (c2_b_i = 0; c2_b_i < 7; c2_b_i++) {
    c2_I[c2_i] = 1.0;
    c2_i += 8;
  }
}

static void c2_mldivide(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  creal_T c2_A[49], real_T c2_B[7], creal_T c2_Y[7])
{
  int32_T c2_i10;
  creal_T c2_b_A[49];
  int8_T c2_ipiv[7];
  int32_T c2_info;
  int32_T c2_j;
  int32_T c2_jj;
  int32_T c2_jp1j;
  int32_T c2_iy;
  int32_T c2_ix;
  real_T c2_smax;
  int32_T c2_b_k;
  real_T c2_s;
  int32_T c2_jrow;
  creal_T c2_temp;
  int32_T c2_jA;
  for (c2_i10 = 0; c2_i10 < 49; c2_i10++) {
    c2_b_A[c2_i10] = c2_A[c2_i10];
  }

  for (c2_i10 = 0; c2_i10 < 7; c2_i10++) {
    c2_ipiv[c2_i10] = (int8_T)(1 + c2_i10);
  }

  c2_info = 0;
  for (c2_j = 0; c2_j < 6; c2_j++) {
    c2_jj = c2_j << 3;
    c2_jp1j = c2_jj + 2;
    c2_iy = 1;
    c2_ix = c2_jj;
    c2_smax = muDoubleScalarAbs(c2_b_A[c2_jj].re) + muDoubleScalarAbs
      (c2_b_A[c2_jj].im);
    for (c2_b_k = 2; c2_b_k <= 7 - c2_j; c2_b_k++) {
      c2_ix++;
      sf_mex_lw_bounds_check(c2_ix + 1, 1, 49);
      c2_s = muDoubleScalarAbs(c2_b_A[c2_ix].re) + muDoubleScalarAbs
        (c2_b_A[c2_ix].im);
      if (c2_s > c2_smax) {
        c2_iy = c2_b_k;
        c2_smax = c2_s;
      }
    }

    if ((c2_b_A[(c2_jj + c2_iy) - 1].re != 0.0) || (c2_b_A[(c2_jj + c2_iy) - 1].
         im != 0.0)) {
      if (c2_iy - 1 != 0) {
        c2_ipiv[c2_j] = (int8_T)(c2_j + c2_iy);
        c2_jrow = 1 + c2_j;
        c2_iy = (c2_jrow + c2_iy) - 1;
        for (c2_b_k = 0; c2_b_k < 7; c2_b_k++) {
          c2_temp = c2_b_A[sf_mex_lw_bounds_check(c2_jrow, 1, 49) - 1];
          c2_b_A[c2_jrow - 1] = c2_b_A[sf_mex_lw_bounds_check(c2_iy, 1, 49) - 1];
          c2_b_A[c2_iy - 1] = c2_temp;
          c2_jrow += 7;
          c2_iy += 7;
        }
      }

      c2_i10 = (c2_jp1j - c2_j) + 5;
      for (c2_jrow = c2_jp1j; c2_jrow <= c2_i10; c2_jrow++) {
        c2_b_A[c2_jrow - 1] = c2_b_eml_div(chartInstance, c2_b_A[c2_jrow - 1],
          c2_b_A[c2_jj]);
      }
    } else {
      c2_info = c2_j + 1;
    }

    c2_jA = c2_jj + 8;
    c2_iy = c2_jj + 7;
    for (c2_jrow = 1; c2_jrow <= 6 - c2_j; c2_jrow++) {
      sf_mex_lw_bounds_check(c2_iy + 1, 1, 49);
      if ((c2_b_A[c2_iy].re != 0.0) || (c2_b_A[c2_iy].im != 0.0)) {
        c2_temp.re = -c2_b_A[c2_iy].re - c2_b_A[c2_iy].im * 0.0;
        c2_temp.im = c2_b_A[c2_iy].re * 0.0 + -c2_b_A[c2_iy].im;
        c2_ix = c2_jp1j;
        c2_i10 = (c2_jA - c2_j) + 6;
        for (c2_b_k = 1 + c2_jA; c2_b_k <= c2_i10; c2_b_k++) {
          c2_smax = c2_b_A[sf_mex_lw_bounds_check(c2_ix, 1, 49) - 1].re *
            c2_temp.re - c2_b_A[sf_mex_lw_bounds_check(c2_ix, 1, 49) - 1].im *
            c2_temp.im;
          c2_s = c2_b_A[sf_mex_lw_bounds_check(c2_ix, 1, 49) - 1].re *
            c2_temp.im + c2_b_A[sf_mex_lw_bounds_check(c2_ix, 1, 49) - 1].im *
            c2_temp.re;
          c2_b_A[sf_mex_lw_bounds_check(c2_b_k, 1, 49) - 1].re =
            c2_b_A[sf_mex_lw_bounds_check(c2_b_k, 1, 49) - 1].re + c2_smax;
          c2_b_A[sf_mex_lw_bounds_check(c2_b_k, 1, 49) - 1].im =
            c2_b_A[sf_mex_lw_bounds_check(c2_b_k, 1, 49) - 1].im + c2_s;
          c2_ix++;
        }
      }

      c2_iy += 7;
      c2_jA += 7;
    }
  }

  if ((c2_info == 0) && (!((c2_b_A[48].re != 0.0) || (c2_b_A[48].im != 0.0)))) {
    c2_info = 7;
  }

  if (c2_info > 0) {
    c2_d_eml_warning(chartInstance);
  }

  for (c2_i10 = 0; c2_i10 < 7; c2_i10++) {
    c2_Y[c2_i10].re = c2_B[c2_i10];
    c2_Y[c2_i10].im = 0.0;
  }

  for (c2_jrow = 0; c2_jrow < 7; c2_jrow++) {
    if (c2_ipiv[c2_jrow] != c2_jrow + 1) {
      c2_temp = c2_Y[c2_jrow];
      c2_Y[c2_jrow] = c2_Y[sf_mex_lw_bounds_check((int32_T)c2_ipiv[c2_jrow], 1,
        7) - 1];
      c2_Y[sf_mex_lw_bounds_check((int32_T)c2_ipiv[c2_jrow], 1, 7) - 1] =
        c2_temp;
    }
  }

  for (c2_b_k = 0; c2_b_k < 7; c2_b_k++) {
    c2_iy = 7 * c2_b_k;
    if ((c2_Y[c2_b_k].re != 0.0) || (c2_Y[c2_b_k].im != 0.0)) {
      for (c2_jrow = c2_b_k + 2; c2_jrow < 8; c2_jrow++) {
        c2_smax = c2_Y[c2_b_k].re * c2_b_A[(c2_jrow + c2_iy) - 1].im +
          c2_Y[c2_b_k].im * c2_b_A[(c2_jrow + c2_iy) - 1].re;
        c2_Y[c2_jrow - 1].re -= c2_Y[c2_b_k].re * c2_b_A[(c2_jrow + c2_iy) - 1].
          re - c2_Y[c2_b_k].im * c2_b_A[(c2_jrow + c2_iy) - 1].im;
        c2_Y[c2_jrow - 1].im -= c2_smax;
      }
    }
  }

  for (c2_b_k = 6; c2_b_k > -1; c2_b_k += -1) {
    c2_iy = 7 * c2_b_k;
    if ((c2_Y[c2_b_k].re != 0.0) || (c2_Y[c2_b_k].im != 0.0)) {
      c2_Y[c2_b_k] = c2_b_eml_div(chartInstance, c2_Y[c2_b_k], c2_b_A[c2_b_k +
        c2_iy]);
      for (c2_jrow = 0; c2_jrow + 1 <= c2_b_k; c2_jrow++) {
        c2_smax = c2_Y[c2_b_k].re * c2_b_A[c2_jrow + c2_iy].im + c2_Y[c2_b_k].im
          * c2_b_A[c2_jrow + c2_iy].re;
        c2_Y[c2_jrow].re -= c2_Y[c2_b_k].re * c2_b_A[c2_jrow + c2_iy].re -
          c2_Y[c2_b_k].im * c2_b_A[c2_jrow + c2_iy].im;
        c2_Y[c2_jrow].im -= c2_smax;
      }
    }
  }
}

static creal_T c2_b_eml_div(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, creal_T c2_x, creal_T c2_y)
{
  creal_T c2_z;
  real_T c2_brm;
  real_T c2_bim;
  real_T c2_d;
  if (c2_y.im == 0.0) {
    if (c2_x.im == 0.0) {
      c2_z.re = c2_x.re / c2_y.re;
      c2_z.im = 0.0;
    } else if (c2_x.re == 0.0) {
      c2_z.re = 0.0;
      c2_z.im = c2_x.im / c2_y.re;
    } else {
      c2_z.re = c2_x.re / c2_y.re;
      c2_z.im = c2_x.im / c2_y.re;
    }
  } else if (c2_y.re == 0.0) {
    if (c2_x.re == 0.0) {
      c2_z.re = c2_x.im / c2_y.im;
      c2_z.im = 0.0;
    } else if (c2_x.im == 0.0) {
      c2_z.re = 0.0;
      c2_z.im = -(c2_x.re / c2_y.im);
    } else {
      c2_z.re = c2_x.im / c2_y.im;
      c2_z.im = -(c2_x.re / c2_y.im);
    }
  } else {
    c2_brm = muDoubleScalarAbs(c2_y.re);
    c2_bim = muDoubleScalarAbs(c2_y.im);
    if (c2_brm > c2_bim) {
      c2_bim = c2_y.im / c2_y.re;
      c2_d = c2_y.re + c2_bim * c2_y.im;
      c2_z.re = (c2_x.re + c2_bim * c2_x.im) / c2_d;
      c2_z.im = (c2_x.im - c2_bim * c2_x.re) / c2_d;
    } else if (c2_bim == c2_brm) {
      c2_bim = c2_y.re > 0.0 ? 0.5 : -0.5;
      c2_d = c2_y.im > 0.0 ? 0.5 : -0.5;
      c2_z.re = (c2_x.re * c2_bim + c2_x.im * c2_d) / c2_brm;
      c2_z.im = (c2_x.im * c2_bim - c2_x.re * c2_d) / c2_brm;
    } else {
      c2_bim = c2_y.re / c2_y.im;
      c2_d = c2_y.im + c2_bim * c2_y.re;
      c2_z.re = (c2_bim * c2_x.re + c2_x.im) / c2_d;
      c2_z.im = (c2_bim * c2_x.im - c2_x.re) / c2_d;
    }
  }

  return c2_z;
}

static void c2_d_eml_warning(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  int32_T c2_i11;
  static char_T c2_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c2_u[27];
  const mxArray *c2_y = NULL;
  for (c2_i11 = 0; c2_i11 < 27; c2_i11++) {
    c2_u[c2_i11] = c2_varargin_1[c2_i11];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static void c2_abs(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   creal_T c2_x[2], real_T c2_y[2])
{
  int32_T c2_b_k;
  real_T c2_a;
  real_T c2_b;
  for (c2_b_k = 0; c2_b_k < 2; c2_b_k++) {
    c2_a = muDoubleScalarAbs(c2_x[c2_b_k].re);
    c2_b = muDoubleScalarAbs(c2_x[c2_b_k].im);
    if (c2_a < c2_b) {
      c2_a /= c2_b;
      c2_b *= muDoubleScalarSqrt(c2_a * c2_a + 1.0);
    } else if (c2_a > c2_b) {
      c2_b /= c2_a;
      c2_b = muDoubleScalarSqrt(c2_b * c2_b + 1.0) * c2_a;
    } else {
      if (!muDoubleScalarIsNaN(c2_b)) {
        c2_b = c2_a * 1.4142135623730951;
      }
    }

    c2_y[c2_b_k] = c2_b;
  }
}

static void c2_b_eml_error(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
  int32_T c2_i12;
  static char_T c2_varargin_1[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'l', 'o', 'g', '1', '0', '_', 'd', 'o', 'm',
    'a', 'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c2_u[31];
  const mxArray *c2_y = NULL;
  for (c2_i12 = 0; c2_i12 < 31; c2_i12++) {
    c2_u[c2_i12] = c2_varargin_1[c2_i12];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 31), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c2_y));
}

static real_T c2_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_Fs1, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_Fs1), &c2_thisId);
  sf_mex_destroy(&c2_Fs1);
  return c2_y;
}

static real_T c2_b_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d1;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d1, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d1;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_c_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_Wn, const char_T *c2_identifier, real_T
  c2_y[7])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_Wn), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_Wn);
}

static void c2_d_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7])
{
  real_T c2_dv20[7];
  int32_T c2_i13;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv20, 1, 0, 0U, 1, 0U, 1, 7);
  for (c2_i13 = 0; c2_i13 < 7; c2_i13++) {
    c2_y[c2_i13] = c2_dv20[c2_i13];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_e_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_ABK, const char_T *c2_identifier, real_T
  c2_y[70])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_ABK), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_ABK);
}

static void c2_f_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[70])
{
  real_T c2_dv21[70];
  int32_T c2_i14;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_ABK_not_empty = FALSE;
  } else {
    chartInstance->c2_ABK_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv21, 1, 0, 0U, 1, 0U, 2, 7,
                  10);
    for (c2_i14 = 0; c2_i14 < 70; c2_i14++) {
      c2_y[c2_i14] = c2_dv21[c2_i14];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_g_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_CD, const char_T *c2_identifier, real_T
  c2_y[16])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_CD), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_CD);
}

static void c2_h_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[16])
{
  real_T c2_dv22[16];
  int32_T c2_i15;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_CD_not_empty = FALSE;
  } else {
    chartInstance->c2_CD_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv22, 1, 0, 0U, 1, 0U, 2, 2,
                  8);
    for (c2_i15 = 0; c2_i15 < 16; c2_i15++) {
      c2_y[c2_i15] = c2_dv22[c2_i15];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_i_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_LK, const char_T *c2_identifier, real_T
  c2_y[15000])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_LK), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_LK);
}

static void c2_j_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[15000])
{
  real_T c2_dv23[15000];
  int32_T c2_i16;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_LK_not_empty = FALSE;
  } else {
    chartInstance->c2_LK_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv23, 1, 0, 0U, 1, 0U, 2,
                  100, 150);
    for (c2_i16 = 0; c2_i16 < 15000; c2_i16++) {
      c2_y[c2_i16] = c2_dv23[c2_i16];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_k_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Pabk, const char_T *c2_identifier, real_T
  c2_y[100])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_Pabk), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_Pabk);
}

static void c2_l_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[100])
{
  real_T c2_dv24[100];
  int32_T c2_i17;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_Pabk_not_empty = FALSE;
  } else {
    chartInstance->c2_Pabk_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv24, 1, 0, 0U, 1, 0U, 2, 10,
                  10);
    for (c2_i17 = 0; c2_i17 < 100; c2_i17++) {
      c2_y[c2_i17] = c2_dv24[c2_i17];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_m_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Pcd, const char_T *c2_identifier, real_T
  c2_y[64])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_Pcd), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_Pcd);
}

static void c2_n_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[64])
{
  real_T c2_dv25[64];
  int32_T c2_i18;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_Pcd_not_empty = FALSE;
  } else {
    chartInstance->c2_Pcd_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv25, 1, 0, 0U, 1, 0U, 2, 8,
                  8);
    for (c2_i18 = 0; c2_i18 < 64; c2_i18++) {
      c2_y[c2_i18] = c2_dv25[c2_i18];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_o_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Plk, const char_T *c2_identifier, real_T
  c2_y[23104])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_Plk), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_Plk);
}

static void c2_p_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[23104])
{
  static real_T c2_dv26[23104];
  int32_T c2_i19;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_Plk_not_empty = FALSE;
  } else {
    chartInstance->c2_Plk_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv26, 1, 0, 0U, 1, 0U, 2,
                  152, 152);
    for (c2_i19 = 0; c2_i19 < 23104; c2_i19++) {
      c2_y[c2_i19] = c2_dv26[c2_i19];
    }
  }

  sf_mex_destroy(&c2_u);
}

static real_T c2_q_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_U1, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_U1), &c2_thisId);
  sf_mex_destroy(&c2_b_U1);
  return c2_y;
}

static real_T c2_r_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d2;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_U1_not_empty = FALSE;
  } else {
    chartInstance->c2_U1_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d2, 1, 0, 0U, 0, 0U, 0);
    c2_y = c2_d2;
  }

  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_s_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_VARX, const char_T *c2_identifier, real_T
  c2_y[304])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_VARX), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_VARX);
}

static void c2_t_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[304])
{
  real_T c2_dv27[304];
  int32_T c2_i20;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_VARX_not_empty = FALSE;
  } else {
    chartInstance->c2_VARX_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv27, 1, 0, 0U, 1, 0U, 2, 2,
                  152);
    for (c2_i20 = 0; c2_i20 < 304; c2_i20++) {
      c2_y[c2_i20] = c2_dv27[c2_i20];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_u_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_X, const char_T *c2_identifier, real_T
  c2_y[7])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_v_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_X), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_X);
}

static void c2_v_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7])
{
  real_T c2_dv28[7];
  int32_T c2_i21;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_X_not_empty = FALSE;
  } else {
    chartInstance->c2_X_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv28, 1, 0, 0U, 1, 0U, 1, 7);
    for (c2_i21 = 0; c2_i21 < 7; c2_i21++) {
      c2_y[c2_i21] = c2_dv28[c2_i21];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_w_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_X1, const char_T *c2_identifier, real_T
  c2_y[7])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_x_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_X1), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_X1);
}

static void c2_x_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[7])
{
  real_T c2_dv29[7];
  int32_T c2_i22;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_X1_not_empty = FALSE;
  } else {
    chartInstance->c2_X1_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv29, 1, 0, 0U, 1, 0U, 1, 7);
    for (c2_i22 = 0; c2_i22 < 7; c2_i22++) {
      c2_y[c2_i22] = c2_dv29[c2_i22];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_y_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Y1, const char_T *c2_identifier, real_T
  c2_y[2])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_ab_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_Y1), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_Y1);
}

static void c2_ab_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[2])
{
  real_T c2_dv30[2];
  int32_T c2_i23;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_Y1_not_empty = FALSE;
  } else {
    chartInstance->c2_Y1_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv30, 1, 0, 0U, 1, 0U, 1, 2);
    for (c2_i23 = 0; c2_i23 < 2; c2_i23++) {
      c2_y[c2_i23] = c2_dv30[c2_i23];
    }
  }

  sf_mex_destroy(&c2_u);
}

static void c2_bb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_Zp, const char_T *c2_identifier, real_T
  c2_y[150])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_cb_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_Zp), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_Zp);
}

static void c2_cb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[150])
{
  real_T c2_dv31[150];
  int32_T c2_i24;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_Zp_not_empty = FALSE;
  } else {
    chartInstance->c2_Zp_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv31, 1, 0, 0U, 1, 0U, 1,
                  150);
    for (c2_i24 = 0; c2_i24 < 150; c2_i24++) {
      c2_y[c2_i24] = c2_dv31[c2_i24];
    }
  }

  sf_mex_destroy(&c2_u);
}

static real_T c2_db_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_k, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_eb_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_k), &c2_thisId);
  sf_mex_destroy(&c2_b_k);
  return c2_y;
}

static real_T c2_eb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d3;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_k_not_empty = FALSE;
  } else {
    chartInstance->c2_k_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d3, 1, 0, 0U, 0, 0U, 0);
    c2_y = c2_d3;
  }

  sf_mex_destroy(&c2_u);
  return c2_y;
}

static real_T c2_fb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_saw, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_gb_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_saw), &c2_thisId);
  sf_mex_destroy(&c2_b_saw);
  return c2_y;
}

static real_T c2_gb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d4;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_saw_not_empty = FALSE;
  } else {
    chartInstance->c2_saw_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d4, 1, 0, 0U, 0, 0U, 0);
    c2_y = c2_d4;
  }

  sf_mex_destroy(&c2_u);
  return c2_y;
}

static real_T c2_hb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_start, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_ib_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_start),
    &c2_thisId);
  sf_mex_destroy(&c2_b_start);
  return c2_y;
}

static real_T c2_ib_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d5;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_start_not_empty = FALSE;
  } else {
    chartInstance->c2_start_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d5, 1, 0, 0U, 0, 0U, 0);
    c2_y = c2_d5;
  }

  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_jb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_w, const char_T *c2_identifier, real_T
  c2_y[300])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_kb_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_w), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_w);
}

static void c2_kb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[300])
{
  real_T c2_dv32[300];
  int32_T c2_i25;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_w_not_empty = FALSE;
  } else {
    chartInstance->c2_w_not_empty = TRUE;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv32, 1, 0, 0U, 1, 0U, 2, 1,
                  300);
    for (c2_i25 = 0; c2_i25 < 300; c2_i25++) {
      c2_y[c2_i25] = c2_dv32[c2_i25];
    }
  }

  sf_mex_destroy(&c2_u);
}

static uint8_T c2_lb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_sim1_rlti_varx_gust2, const
  char_T *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_mb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_sim1_rlti_varx_gust2), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_sim1_rlti_varx_gust2);
  return c2_y;
}

static uint8_T c2_mb_emlrt_marshallIn(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_rls_ew_reg(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_z[152], real_T c2_y[2], real_T c2_theta[304], real_T c2_P[23104],
  real_T c2_lambda, real_T c2_reg, real_T *c2_b_saw)
{
  int32_T c2_i26;
  real_T c2_eta[152];
  real_T c2_b_y;
  int32_T c2_m;
  int32_T c2_n;
  static real_T c2_c_y[23104];
  int32_T c2_b_k;
  real_T c2_alpha1;
  int32_T c2_lda;
  int32_T c2_ldb;
  real_T c2_beta1;
  int32_T c2_ldc;
  char_T c2_TRANSA;
  char_T c2_TRANSB;
  static real_T c2_d_y[23104];
  real_T c2_a[152];
  real_T c2_e_y[152];
  real_T c2_f_y[2];
  real_T c2_g_y[304];
  real_T c2_h_y[304];
  for (c2_i26 = 0; c2_i26 < 152; c2_i26++) {
    c2_eta[c2_i26] = 0.0;
  }

  if (*c2_b_saw > 151.5) {
    *c2_b_saw = 1.0;
  } else {
    (*c2_b_saw)++;
  }

  c2_b_y = (1.0 - muDoubleScalarPower(c2_lambda, 152.0)) * c2_reg;
  if (c2_b_y < 0.0) {
    c2_eml_error(chartInstance);
  }

  c2_eta[sf_mex_lw_bounds_check((int32_T)*c2_b_saw, 1, 152) - 1] =
    muDoubleScalarSqrt(c2_b_y);
  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_c_y[c2_n + c2_i26] = c2_eta[c2_n] * c2_eta[c2_m];
    }

    c2_i26 += 152;
  }

  c2_m = 152;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 152;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 152;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 23104; c2_i26++) {
    c2_d_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_P[0],
          &c2_lda, &c2_c_y[0], &c2_ldb, &c2_beta1, &c2_d_y[0], &c2_ldc);
  c2_m = 152;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 152;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 152;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 23104; c2_i26++) {
    c2_c_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_d_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_c_y[0], &c2_ldc);
  c2_m = 1;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 1;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 1;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 152; c2_i26++) {
    c2_a[c2_i26] = c2_eta[c2_i26];
    c2_e_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_a[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_e_y[0], &c2_ldc);
  c2_b_y = 0.0;
  for (c2_b_k = 0; c2_b_k < 152; c2_b_k++) {
    c2_b_y += c2_e_y[c2_b_k] * c2_eta[c2_b_k];
  }

  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_P[c2_n + c2_i26] -= c2_c_y[c2_n + c2_i26] / (1.0 + c2_b_y);
    }

    c2_i26 += 152;
  }

  c2_b_y = 1.0 / c2_lambda;
  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_c_y[c2_n + c2_i26] = c2_z[c2_n] * c2_z[c2_m];
    }

    c2_i26 += 152;
  }

  c2_m = 152;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 152;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 152;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 23104; c2_i26++) {
    c2_d_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_P[0],
          &c2_lda, &c2_c_y[0], &c2_ldb, &c2_beta1, &c2_d_y[0], &c2_ldc);
  c2_m = 152;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 152;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 152;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 23104; c2_i26++) {
    c2_c_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_d_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_c_y[0], &c2_ldc);
  c2_m = 1;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 1;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 1;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 152; c2_i26++) {
    c2_a[c2_i26] = c2_z[c2_i26];
    c2_e_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_a[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_e_y[0], &c2_ldc);
  c2_alpha1 = 0.0;
  for (c2_b_k = 0; c2_b_k < 152; c2_b_k++) {
    c2_alpha1 += c2_e_y[c2_b_k] * c2_z[c2_b_k];
  }

  c2_alpha1 += c2_lambda;
  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_P[c2_n + c2_i26] = c2_b_y * (c2_P[c2_n + c2_i26] - c2_c_y[c2_n + c2_i26]
        / c2_alpha1);
    }

    c2_i26 += 152;
  }

  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    c2_n = 0;
    for (c2_b_k = 0; c2_b_k < 152; c2_b_k++) {
      c2_c_y[c2_b_k + c2_i26] = 0.5 * (c2_P[c2_b_k + c2_i26] + c2_P[c2_n + c2_m]);
      c2_n += 152;
    }

    c2_i26 += 152;
  }

  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_P[c2_n + c2_i26] = c2_c_y[c2_n + c2_i26];
    }

    c2_i26 += 152;
  }

  for (c2_i26 = 0; c2_i26 < 2; c2_i26++) {
    c2_alpha1 = 0.0;
    c2_m = 0;
    for (c2_n = 0; c2_n < 152; c2_n++) {
      c2_alpha1 += c2_theta[c2_m + c2_i26] * c2_z[c2_n];
      c2_m += 2;
    }

    c2_f_y[c2_i26] = c2_y[c2_i26] - c2_alpha1;
  }

  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 2; c2_n++) {
      c2_g_y[c2_n + c2_i26] = c2_f_y[c2_n] * c2_z[c2_m];
    }

    c2_i26 += 2;
  }

  c2_m = 2;
  c2_n = 152;
  c2_b_k = 152;
  c2_alpha1 = 1.0;
  c2_lda = 2;
  c2_ldb = 152;
  c2_beta1 = 0.0;
  c2_ldc = 2;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i26 = 0; c2_i26 < 304; c2_i26++) {
    c2_h_y[c2_i26] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_g_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_h_y[0], &c2_ldc);
  c2_i26 = 0;
  for (c2_m = 0; c2_m < 152; c2_m++) {
    for (c2_n = 0; c2_n < 2; c2_n++) {
      c2_theta[c2_n + c2_i26] += c2_h_y[c2_n + c2_i26];
    }

    c2_i26 += 2;
  }
}

static void c2_rls_ew(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                      real_T c2_z[8], real_T c2_y[2], real_T c2_theta[16],
                      real_T c2_P[64], real_T c2_lambda)
{
  real_T c2_b_y;
  int32_T c2_i27;
  int32_T c2_m;
  int32_T c2_n;
  real_T c2_c_y[64];
  int32_T c2_b_k;
  real_T c2_alpha1;
  int32_T c2_lda;
  int32_T c2_ldb;
  real_T c2_beta1;
  int32_T c2_ldc;
  char_T c2_TRANSA;
  char_T c2_TRANSB;
  real_T c2_d_y[64];
  real_T c2_e_y[8];
  real_T c2_f_y[2];
  real_T c2_g_y[16];
  real_T c2_b_theta[16];
  c2_b_y = 1.0 / c2_lambda;
  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    for (c2_n = 0; c2_n < 8; c2_n++) {
      c2_c_y[c2_n + c2_i27] = c2_z[c2_n] * c2_z[c2_m];
    }

    c2_i27 += 8;
  }

  c2_m = 8;
  c2_n = 8;
  c2_b_k = 8;
  c2_alpha1 = 1.0;
  c2_lda = 8;
  c2_ldb = 8;
  c2_beta1 = 0.0;
  c2_ldc = 8;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i27 = 0; c2_i27 < 64; c2_i27++) {
    c2_d_y[c2_i27] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_P[0],
          &c2_lda, &c2_c_y[0], &c2_ldb, &c2_beta1, &c2_d_y[0], &c2_ldc);
  c2_m = 8;
  c2_n = 8;
  c2_b_k = 8;
  c2_alpha1 = 1.0;
  c2_lda = 8;
  c2_ldb = 8;
  c2_beta1 = 0.0;
  c2_ldc = 8;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i27 = 0; c2_i27 < 64; c2_i27++) {
    c2_c_y[c2_i27] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_d_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_c_y[0], &c2_ldc);
  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    c2_e_y[c2_m] = 0.0;
    for (c2_n = 0; c2_n < 8; c2_n++) {
      c2_alpha1 = c2_e_y[c2_m] + c2_z[c2_n] * c2_P[c2_n + c2_i27];
      c2_e_y[c2_m] = c2_alpha1;
    }

    c2_i27 += 8;
  }

  c2_alpha1 = 0.0;
  for (c2_b_k = 0; c2_b_k < 8; c2_b_k++) {
    c2_alpha1 += c2_e_y[c2_b_k] * c2_z[c2_b_k];
  }

  c2_alpha1 += c2_lambda;
  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    for (c2_n = 0; c2_n < 8; c2_n++) {
      c2_P[c2_n + c2_i27] = c2_b_y * (c2_P[c2_n + c2_i27] - c2_c_y[c2_n + c2_i27]
        / c2_alpha1);
    }

    c2_i27 += 8;
  }

  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    c2_n = 0;
    for (c2_b_k = 0; c2_b_k < 8; c2_b_k++) {
      c2_c_y[c2_b_k + c2_i27] = 0.5 * (c2_P[c2_b_k + c2_i27] + c2_P[c2_n + c2_m]);
      c2_n += 8;
    }

    c2_i27 += 8;
  }

  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    for (c2_n = 0; c2_n < 8; c2_n++) {
      c2_P[c2_n + c2_i27] = c2_c_y[c2_n + c2_i27];
    }

    c2_i27 += 8;
  }

  for (c2_i27 = 0; c2_i27 < 2; c2_i27++) {
    c2_alpha1 = 0.0;
    c2_m = 0;
    for (c2_n = 0; c2_n < 8; c2_n++) {
      c2_alpha1 += c2_theta[c2_m + c2_i27] * c2_z[c2_n];
      c2_m += 2;
    }

    c2_f_y[c2_i27] = c2_y[c2_i27] - c2_alpha1;
  }

  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    for (c2_n = 0; c2_n < 2; c2_n++) {
      c2_g_y[c2_n + c2_i27] = c2_f_y[c2_n] * c2_z[c2_m];
    }

    c2_i27 += 2;
  }

  for (c2_i27 = 0; c2_i27 < 2; c2_i27++) {
    c2_m = 0;
    c2_n = 0;
    for (c2_b_k = 0; c2_b_k < 8; c2_b_k++) {
      c2_alpha1 = 0.0;
      c2_lda = 0;
      for (c2_ldb = 0; c2_ldb < 8; c2_ldb++) {
        c2_alpha1 += c2_g_y[c2_lda + c2_i27] * c2_P[c2_ldb + c2_n];
        c2_lda += 2;
      }

      c2_b_theta[c2_m + c2_i27] = c2_theta[c2_m + c2_i27] + c2_alpha1;
      c2_m += 2;
      c2_n += 8;
    }
  }

  c2_i27 = 0;
  for (c2_m = 0; c2_m < 8; c2_m++) {
    for (c2_n = 0; c2_n < 2; c2_n++) {
      c2_theta[c2_n + c2_i27] = c2_b_theta[c2_n + c2_i27];
    }

    c2_i27 += 2;
  }
}

static void c2_b_rls_ew(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
  real_T c2_z[10], real_T c2_y[7], real_T c2_theta[70], real_T c2_P[100], real_T
  c2_lambda)
{
  real_T c2_b_y;
  int32_T c2_i28;
  int32_T c2_m;
  int32_T c2_n;
  real_T c2_c_y[100];
  int32_T c2_b_k;
  real_T c2_alpha1;
  int32_T c2_lda;
  int32_T c2_ldb;
  real_T c2_beta1;
  int32_T c2_ldc;
  char_T c2_TRANSA;
  char_T c2_TRANSB;
  real_T c2_d_y[100];
  real_T c2_e_y[10];
  real_T c2_f_y[7];
  real_T c2_g_y[70];
  real_T c2_h_y[70];
  c2_b_y = 1.0 / c2_lambda;
  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    for (c2_n = 0; c2_n < 10; c2_n++) {
      c2_c_y[c2_n + c2_i28] = c2_z[c2_n] * c2_z[c2_m];
    }

    c2_i28 += 10;
  }

  c2_m = 10;
  c2_n = 10;
  c2_b_k = 10;
  c2_alpha1 = 1.0;
  c2_lda = 10;
  c2_ldb = 10;
  c2_beta1 = 0.0;
  c2_ldc = 10;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i28 = 0; c2_i28 < 100; c2_i28++) {
    c2_d_y[c2_i28] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_P[0],
          &c2_lda, &c2_c_y[0], &c2_ldb, &c2_beta1, &c2_d_y[0], &c2_ldc);
  c2_m = 10;
  c2_n = 10;
  c2_b_k = 10;
  c2_alpha1 = 1.0;
  c2_lda = 10;
  c2_ldb = 10;
  c2_beta1 = 0.0;
  c2_ldc = 10;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i28 = 0; c2_i28 < 100; c2_i28++) {
    c2_c_y[c2_i28] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_d_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_c_y[0], &c2_ldc);
  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    c2_e_y[c2_m] = 0.0;
    for (c2_n = 0; c2_n < 10; c2_n++) {
      c2_alpha1 = c2_e_y[c2_m] + c2_z[c2_n] * c2_P[c2_n + c2_i28];
      c2_e_y[c2_m] = c2_alpha1;
    }

    c2_i28 += 10;
  }

  c2_alpha1 = 0.0;
  for (c2_b_k = 0; c2_b_k < 10; c2_b_k++) {
    c2_alpha1 += c2_e_y[c2_b_k] * c2_z[c2_b_k];
  }

  c2_alpha1 += c2_lambda;
  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    for (c2_n = 0; c2_n < 10; c2_n++) {
      c2_P[c2_n + c2_i28] = c2_b_y * (c2_P[c2_n + c2_i28] - c2_c_y[c2_n + c2_i28]
        / c2_alpha1);
    }

    c2_i28 += 10;
  }

  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    c2_n = 0;
    for (c2_b_k = 0; c2_b_k < 10; c2_b_k++) {
      c2_c_y[c2_b_k + c2_i28] = 0.5 * (c2_P[c2_b_k + c2_i28] + c2_P[c2_n + c2_m]);
      c2_n += 10;
    }

    c2_i28 += 10;
  }

  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    for (c2_n = 0; c2_n < 10; c2_n++) {
      c2_P[c2_n + c2_i28] = c2_c_y[c2_n + c2_i28];
    }

    c2_i28 += 10;
  }

  for (c2_i28 = 0; c2_i28 < 7; c2_i28++) {
    c2_alpha1 = 0.0;
    c2_m = 0;
    for (c2_n = 0; c2_n < 10; c2_n++) {
      c2_alpha1 += c2_theta[c2_m + c2_i28] * c2_z[c2_n];
      c2_m += 7;
    }

    c2_f_y[c2_i28] = c2_y[c2_i28] - c2_alpha1;
  }

  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    for (c2_n = 0; c2_n < 7; c2_n++) {
      c2_g_y[c2_n + c2_i28] = c2_f_y[c2_n] * c2_z[c2_m];
    }

    c2_i28 += 7;
  }

  c2_m = 7;
  c2_n = 10;
  c2_b_k = 10;
  c2_alpha1 = 1.0;
  c2_lda = 7;
  c2_ldb = 10;
  c2_beta1 = 0.0;
  c2_ldc = 7;
  c2_TRANSA = 'N';
  c2_TRANSB = 'N';
  for (c2_i28 = 0; c2_i28 < 70; c2_i28++) {
    c2_h_y[c2_i28] = 0.0;
  }

  dgemm32(&c2_TRANSA, &c2_TRANSB, &c2_m, &c2_n, &c2_b_k, &c2_alpha1, &c2_g_y[0],
          &c2_lda, &c2_P[0], &c2_ldb, &c2_beta1, &c2_h_y[0], &c2_ldc);
  c2_i28 = 0;
  for (c2_m = 0; c2_m < 10; c2_m++) {
    for (c2_n = 0; c2_n < 7; c2_n++) {
      c2_theta[c2_n + c2_i28] += c2_h_y[c2_n + c2_i28];
    }

    c2_i28 += 7;
  }
}

static void c2_sqrt(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                    creal_T *c2_x)
{
  real_T c2_absxi;
  real_T c2_a;
  real_T c2_absxr;
  if (c2_x->im == 0.0) {
    if (c2_x->re < 0.0) {
      c2_absxi = 0.0;
      c2_a = muDoubleScalarSqrt(muDoubleScalarAbs(c2_x->re));
    } else {
      c2_absxi = muDoubleScalarSqrt(c2_x->re);
      c2_a = 0.0;
    }
  } else if (c2_x->re == 0.0) {
    if (c2_x->im < 0.0) {
      c2_absxi = muDoubleScalarSqrt(-c2_x->im / 2.0);
      c2_a = -c2_absxi;
    } else {
      c2_absxi = muDoubleScalarSqrt(c2_x->im / 2.0);
      c2_a = c2_absxi;
    }
  } else if (muDoubleScalarIsNaN(c2_x->re) || muDoubleScalarIsNaN(c2_x->im)) {
    c2_absxi = rtNaN;
    c2_a = rtNaN;
  } else if (muDoubleScalarIsInf(c2_x->im)) {
    c2_absxi = rtInf;
    c2_a = c2_x->im;
  } else if (muDoubleScalarIsInf(c2_x->re)) {
    if (c2_x->re < 0.0) {
      c2_absxi = 0.0;
      c2_a = rtInf;
    } else {
      c2_absxi = rtInf;
      c2_a = 0.0;
    }
  } else {
    c2_absxr = muDoubleScalarAbs(c2_x->re);
    c2_absxi = muDoubleScalarAbs(c2_x->im);
    if ((c2_absxr > 4.4942328371557893E+307) || (c2_absxi >
         4.4942328371557893E+307)) {
      c2_absxr *= 0.5;
      c2_absxi *= 0.5;
      if (c2_absxr < c2_absxi) {
        c2_a = c2_absxr / c2_absxi;
        c2_absxi *= muDoubleScalarSqrt(c2_a * c2_a + 1.0);
      } else if (c2_absxr > c2_absxi) {
        c2_absxi /= c2_absxr;
        c2_absxi = muDoubleScalarSqrt(c2_absxi * c2_absxi + 1.0) * c2_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c2_absxi)) {
          c2_absxi = c2_absxr * 1.4142135623730951;
        }
      }

      if (c2_absxi > c2_absxr) {
        c2_absxi = muDoubleScalarSqrt(c2_absxi) * muDoubleScalarSqrt(1.0 +
          c2_absxr / c2_absxi);
      } else {
        c2_absxi = muDoubleScalarSqrt(c2_absxi) * 1.4142135623730951;
      }
    } else {
      if (c2_absxr < c2_absxi) {
        c2_a = c2_absxr / c2_absxi;
        c2_absxi *= muDoubleScalarSqrt(c2_a * c2_a + 1.0);
      } else if (c2_absxr > c2_absxi) {
        c2_absxi /= c2_absxr;
        c2_absxi = muDoubleScalarSqrt(c2_absxi * c2_absxi + 1.0) * c2_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c2_absxi)) {
          c2_absxi = c2_absxr * 1.4142135623730951;
        }
      }

      c2_absxi = muDoubleScalarSqrt((c2_absxi + c2_absxr) * 0.5);
    }

    if (c2_x->re > 0.0) {
      c2_a = 0.5 * (c2_x->im / c2_absxi);
    } else {
      if (c2_x->im < 0.0) {
        c2_a = -c2_absxi;
      } else {
        c2_a = c2_absxi;
      }

      c2_absxi = 0.5 * (c2_x->im / c2_a);
    }
  }

  c2_x->re = c2_absxi;
  c2_x->im = c2_a;
}

static void c2_exp(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                   creal_T *c2_x)
{
  real_T c2_r;
  creal_T c2_b_x;
  creal_T c2_c_x;
  c2_r = muDoubleScalarExp(c2_x->re / 2.0);
  c2_b_x = *c2_x;
  c2_c_x = *c2_x;
  c2_x->re = c2_r * (c2_r * muDoubleScalarCos(c2_b_x.im));
  c2_x->im = c2_r * (c2_r * muDoubleScalarSin(c2_c_x.im));
}

static void c2_log10(SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance,
                     real_T *c2_x)
{
  if (*c2_x < 0.0) {
    c2_b_eml_error(chartInstance);
  }

  *c2_x = muDoubleScalarLog10(*c2_x);
}

static int32_T c2_div_s32_floor(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance, int32_T c2_numerator, int32_T c2_denominator)
{
  int32_T c2_quotient;
  uint32_T c2_absNumerator;
  uint32_T c2_absDenominator;
  int32_T c2_quotientNeedsNegation;
  uint32_T c2_tempAbsQuotient;
  if (c2_denominator == 0) {
    c2_quotient = c2_numerator >= 0 ? MAX_int32_T : MIN_int32_T;
    sf_mex_dividebyzero_error();
  } else {
    c2_absNumerator = (uint32_T)(c2_numerator >= 0 ? c2_numerator :
      -c2_numerator);
    c2_absDenominator = (uint32_T)(c2_denominator >= 0 ? c2_denominator :
      -c2_denominator);
    c2_quotientNeedsNegation = (c2_numerator < 0 != c2_denominator < 0);
    c2_tempAbsQuotient = c2_absNumerator / c2_absDenominator;
    if ((uint32_T)c2_quotientNeedsNegation) {
      c2_absNumerator %= c2_absDenominator;
      if (c2_absNumerator > (uint32_T)0) {
        c2_tempAbsQuotient++;
      }
    }

    c2_quotient = (uint32_T)c2_quotientNeedsNegation ? -(int32_T)
      c2_tempAbsQuotient : (int32_T)c2_tempAbsQuotient;
  }

  return c2_quotient;
}

static void init_dsm_address_info(SFc2_sim1_rlti_varx_gust2InstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
void sf_c2_sim1_rlti_varx_gust2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4094331313U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3399383051U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1142657676U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4282830528U);
}

mxArray *sf_c2_sim1_rlti_varx_gust2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("kq27u6EyuNq9Now62ZJQsF");
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

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

static const mxArray *sf_get_sim_state_info_c2_sim1_rlti_varx_gust2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x10'type','srcId','name','auxInfo'{{M[1],M[5],T\"Fs1\",},{M[1],M[6],T\"Fs2\",},{M[1],M[9],T\"Wn\",},{M[1],M[4],T\"Ws\",},{M[1],M[10],T\"Zn\",},{M[4],M[0],T\"ABK\",S'l','i','p'{{M1x2[1209 1212],M[0],}}},{M[4],M[0],T\"CD\",S'l','i','p'{{M1x2[1213 1215],M[0],}}},{M[4],M[0],T\"LK\",S'l','i','p'{{M1x2[1197 1199],M[0],}}},{M[4],M[0],T\"Pabk\",S'l','i','p'{{M1x2[1204 1208],M[0],}}},{M[4],M[0],T\"Pcd\",S'l','i','p'{{M1x2[1200 1203],M[0],}}}}",
    "100 S1x10'type','srcId','name','auxInfo'{{M[4],M[0],T\"Plk\",S'l','i','p'{{M1x2[1188 1191],M[0],}}},{M[4],M[0],T\"U1\",S'l','i','p'{{M1x2[1216 1218],M[0],}}},{M[4],M[0],T\"VARX\",S'l','i','p'{{M1x2[1192 1196],M[0],}}},{M[4],M[0],T\"X\",S'l','i','p'{{M1x2[1222 1223],M[0],}}},{M[4],M[0],T\"X1\",S'l','i','p'{{M1x2[1224 1226],M[0],}}},{M[4],M[0],T\"Y1\",S'l','i','p'{{M1x2[1219 1221],M[0],}}},{M[4],M[0],T\"Zp\",S'l','i','p'{{M1x2[1227 1229],M[0],}}},{M[4],M[0],T\"k\",S'l','i','p'{{M1x2[1232 1233],M[0],}}},{M[4],M[0],T\"saw\",S'l','i','p'{{M1x2[1240 1243],M[0],}}},{M[4],M[0],T\"start\",S'l','i','p'{{M1x2[1234 1239],M[0],}}}}",
    "100 S1x2'type','srcId','name','auxInfo'{{M[4],M[0],T\"w\",S'l','i','p'{{M1x2[1230 1231],M[0],}}},{M[8],M[0],T\"is_active_c2_sim1_rlti_varx_gust2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 22, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_sim1_rlti_varx_gust2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void sf_opaque_initialize_c2_sim1_rlti_varx_gust2(void *chartInstanceVar)
{
  initialize_params_c2_sim1_rlti_varx_gust2
    ((SFc2_sim1_rlti_varx_gust2InstanceStruct*) chartInstanceVar);
  initialize_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c2_sim1_rlti_varx_gust2(void *chartInstanceVar)
{
  enable_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c2_sim1_rlti_varx_gust2(void *chartInstanceVar)
{
  disable_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c2_sim1_rlti_varx_gust2(void *chartInstanceVar)
{
  sf_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c2_sim1_rlti_varx_gust2
  (SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_sim1_rlti_varx_gust2
    ((SFc2_sim1_rlti_varx_gust2InstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_sim1_rlti_varx_gust2();/* state var info */
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

extern void sf_internal_set_sim_state_c2_sim1_rlti_varx_gust2(SimStruct* S,
  const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_sim1_rlti_varx_gust2();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_sim1_rlti_varx_gust2(SimStruct*
  S)
{
  return sf_internal_get_sim_state_c2_sim1_rlti_varx_gust2(S);
}

static void sf_opaque_set_sim_state_c2_sim1_rlti_varx_gust2(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c2_sim1_rlti_varx_gust2(S, st);
}

static void sf_opaque_terminate_c2_sim1_rlti_varx_gust2(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_sim1_rlti_varx_gust2InstanceStruct*) chartInstanceVar
      )->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
    }

    finalize_c2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
      chartInstanceVar);
    free((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }

  unload_sim1_rlti_varx_gust2_optimization_info();
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_sim1_rlti_varx_gust2((SFc2_sim1_rlti_varx_gust2InstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_sim1_rlti_varx_gust2(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_sim1_rlti_varx_gust2
      ((SFc2_sim1_rlti_varx_gust2InstanceStruct*)(((ChartInfoStruct *)
         ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_sim1_rlti_varx_gust2(SimStruct *S)
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
    mxArray *infoStruct = load_sim1_rlti_varx_gust2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,infoStruct,2,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,infoStruct,2,4);
      sf_mark_chart_reusable_outputs(S,infoStruct,2,5);
    }

    sf_set_rtw_dwork_info(S,infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(934508857U));
  ssSetChecksum1(S,(1561996025U));
  ssSetChecksum2(S,(4007073007U));
  ssSetChecksum3(S,(2989630701U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
}

static void mdlRTW_c2_sim1_rlti_varx_gust2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_sim1_rlti_varx_gust2(SimStruct *S)
{
  SFc2_sim1_rlti_varx_gust2InstanceStruct *chartInstance;
  chartInstance = (SFc2_sim1_rlti_varx_gust2InstanceStruct *)malloc(sizeof
    (SFc2_sim1_rlti_varx_gust2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_sim1_rlti_varx_gust2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_sim1_rlti_varx_gust2;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c2_sim1_rlti_varx_gust2;
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

void c2_sim1_rlti_varx_gust2_method_dispatcher(SimStruct *S, int_T method, void *
  data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_sim1_rlti_varx_gust2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_sim1_rlti_varx_gust2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_sim1_rlti_varx_gust2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_sim1_rlti_varx_gust2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
