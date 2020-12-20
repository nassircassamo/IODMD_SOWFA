/* Include files */

#include "blascompat32.h"
#include "sim1_rlti_varmax_gust2_sfun.h"
#include "c5_sim1_rlti_varmax_gust2.h"
#include "mwmathutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void initialize_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void initialize_params_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void enable_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void disable_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void set_sim_state_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance, const mxArray
   *c5_st);
static void finalize_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void sf_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void c5_chartstep_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void initSimStructsc5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber);
static void c5_info_helper(c5_ResolvedFunctionInfo c5_info[106]);
static void c5_b_info_helper(c5_ResolvedFunctionInfo c5_info[106]);
static real_T c5_mrdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_A, real_T c5_B);
static real_T c5_rdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_x, real_T c5_y);
static void c5_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   real_T c5_I[63504]);
static void c5_b_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[100]);
static void c5_c_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[64]);
static void c5_logspace(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  real_T c5_d1, real_T c5_d2, real_T c5_y[300]);
static void c5_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static void c5_b_rdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_x, real_T c5_y[2], real_T c5_z[2]);
static void c5_diag(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    real_T c5_v[2], real_T c5_d[4]);
static void c5_eml_error(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static void c5_damp(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    real_T c5_a[49], real_T c5_h, real_T c5_wn[7], real_T c5_z[7]);
static void c5_eig(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   real_T c5_A[49], creal_T c5_V[7]);
static void c5_eml_matlab_zlartg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn,
  creal_T *c5_r);
static void c5_eml_matlab_zhgeqz(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_A[49], int32_T c5_ilo, int32_T c5_ihi, real_T
  *c5_info, creal_T c5_alpha1[7], creal_T c5_beta1[7]);
static real_T c5_eml_matlab_zlanhs(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_A[49], int32_T c5_ilo, int32_T c5_ihi);
static creal_T c5_eml_div(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_x, real_T c5_y);
static void c5_b_eml_matlab_zlartg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn);
static void c5_b_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static void c5_c_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static void c5_d_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[49]);
static void c5_mldivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  creal_T c5_A[49], real_T c5_B[7], creal_T c5_Y[7]);
static creal_T c5_b_eml_div(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_x, creal_T c5_y);
static void c5_d_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static void c5_abs(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   creal_T c5_x[2], real_T c5_y[2]);
static void c5_b_eml_error(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);
static real_T c5_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_Fs1, const char_T *c5_identifier);
static real_T c5_b_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_c_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_Wn, const char_T *c5_identifier, real_T
  c5_y[7]);
static void c5_d_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7]);
static void c5_e_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_ABK, const char_T *c5_identifier, real_T
  c5_y[70]);
static void c5_f_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[70]);
static void c5_g_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_CD, const char_T *c5_identifier, real_T
  c5_y[16]);
static void c5_h_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[16]);
static void c5_i_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_E1, const char_T *c5_identifier, real_T
  c5_y[2]);
static void c5_j_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[2]);
static void c5_k_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_LK, const char_T *c5_identifier, real_T
  c5_y[25000]);
static void c5_l_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[25000]);
static void c5_m_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Pabk, const char_T *c5_identifier, real_T
  c5_y[100]);
static void c5_n_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[100]);
static void c5_o_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Pcd, const char_T *c5_identifier, real_T
  c5_y[64]);
static void c5_p_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[64]);
static void c5_q_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Plk, const char_T *c5_identifier, real_T
  c5_y[63504]);
static void c5_r_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[63504]);
static real_T c5_s_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_U1, const char_T *c5_identifier);
static real_T c5_t_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_u_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_VARX, const char_T *c5_identifier, real_T
  c5_y[504]);
static void c5_v_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[504]);
static void c5_w_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_X, const char_T *c5_identifier, real_T
  c5_y[7]);
static void c5_x_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7]);
static void c5_y_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_X1, const char_T *c5_identifier, real_T
  c5_y[7]);
static void c5_ab_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7]);
static void c5_bb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Y1, const char_T *c5_identifier, real_T
  c5_y[2]);
static void c5_cb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[2]);
static void c5_db_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Zp, const char_T *c5_identifier, real_T
  c5_y[250]);
static void c5_eb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[250]);
static real_T c5_fb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_k, const char_T *c5_identifier);
static real_T c5_gb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static real_T c5_hb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_saw, const char_T *c5_identifier);
static real_T c5_ib_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static real_T c5_jb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_start, const char_T *c5_identifier);
static real_T c5_kb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_lb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_w, const char_T *c5_identifier, real_T
  c5_y[300]);
static void c5_mb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[300]);
static uint8_T c5_nb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct *
  chartInstance, const mxArray *c5_b_is_active_c5_sim1_rlti_varmax_gust2, const
  char_T *c5_identifier);
static uint8_T c5_ob_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct *
  chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_rls_ew_reg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_z[252], real_T c5_y[2], real_T c5_theta[504], real_T
  c5_P[63504], real_T c5_lambda, real_T c5_reg, real_T *c5_b_saw);
static void c5_rls_ew(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                      real_T c5_z[8], real_T c5_y[2], real_T c5_theta[16],
                      real_T c5_P[64], real_T c5_lambda);
static void c5_b_rls_ew(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  real_T c5_z[10], real_T c5_y[7], real_T c5_theta[70], real_T c5_P[100], real_T
  c5_lambda);
static void c5_sqrt(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    creal_T *c5_x);
static void c5_exp(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   creal_T *c5_x);
static void c5_log10(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T *c5_x);
static int32_T c5_div_s32_floor(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, int32_T c5_numerator, int32_T c5_denominator);
static void init_dsm_address_info(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c5_Plk_not_empty = FALSE;
  chartInstance->c5_VARX_not_empty = FALSE;
  chartInstance->c5_LK_not_empty = FALSE;
  chartInstance->c5_Pcd_not_empty = FALSE;
  chartInstance->c5_Pabk_not_empty = FALSE;
  chartInstance->c5_ABK_not_empty = FALSE;
  chartInstance->c5_CD_not_empty = FALSE;
  chartInstance->c5_U1_not_empty = FALSE;
  chartInstance->c5_Y1_not_empty = FALSE;
  chartInstance->c5_X_not_empty = FALSE;
  chartInstance->c5_X1_not_empty = FALSE;
  chartInstance->c5_E1_not_empty = FALSE;
  chartInstance->c5_Zp_not_empty = FALSE;
  chartInstance->c5_w_not_empty = FALSE;
  chartInstance->c5_k_not_empty = FALSE;
  chartInstance->c5_start_not_empty = FALSE;
  chartInstance->c5_saw_not_empty = FALSE;
  chartInstance->c5_is_active_c5_sim1_rlti_varmax_gust2 = 0U;
}

static void initialize_params_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  real_T c5_dv0[700];
  int32_T c5_i0;
  real_T c5_d0;
  real_T c5_dv1[4];
  sf_set_error_prefix_string(
    "Error evaluating data 'W' in the parent workspace.\n");
  sf_mex_import_named("W", sf_mex_get_sfun_param(chartInstance->S, 2, 0), c5_dv0,
                      0, 0, 0U, 1, 0U, 2, 7, 100);
  for (c5_i0 = 0; c5_i0 < 700; c5_i0++) {
    chartInstance->c5_W[c5_i0] = c5_dv0[c5_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Du' in the parent workspace.\n");
  sf_mex_import_named("Du", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c5_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c5_Du = c5_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Dy' in the parent workspace.\n");
  sf_mex_import_named("Dy", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      c5_dv1, 0, 0, 0U, 1, 0U, 2, 2, 2);
  for (c5_i0 = 0; c5_i0 < 4; c5_i0++) {
    chartInstance->c5_Dy[c5_i0] = c5_dv1[c5_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static const mxArray *get_sim_state_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  const mxArray *c5_st;
  const mxArray *c5_y = NULL;
  real_T c5_u;
  const mxArray *c5_b_y = NULL;
  real_T c5_b_u;
  const mxArray *c5_c_y = NULL;
  int32_T c5_i1;
  real_T c5_c_u[7];
  const mxArray *c5_d_y = NULL;
  real_T c5_d_u;
  const mxArray *c5_e_y = NULL;
  real_T c5_e_u[7];
  const mxArray *c5_f_y = NULL;
  real_T c5_f_u[70];
  const mxArray *c5_g_y = NULL;
  real_T c5_g_u[16];
  const mxArray *c5_h_y = NULL;
  real_T c5_h_u[2];
  const mxArray *c5_i_y = NULL;
  static real_T c5_i_u[25000];
  const mxArray *c5_j_y = NULL;
  real_T c5_j_u[100];
  const mxArray *c5_k_y = NULL;
  real_T c5_k_u[64];
  const mxArray *c5_l_y = NULL;
  static real_T c5_l_u[63504];
  const mxArray *c5_m_y = NULL;
  real_T c5_m_u;
  const mxArray *c5_n_y = NULL;
  real_T c5_n_u[504];
  const mxArray *c5_o_y = NULL;
  real_T c5_o_u[7];
  const mxArray *c5_p_y = NULL;
  real_T c5_p_u[7];
  const mxArray *c5_q_y = NULL;
  real_T c5_q_u[2];
  const mxArray *c5_r_y = NULL;
  real_T c5_r_u[250];
  const mxArray *c5_s_y = NULL;
  real_T c5_s_u;
  const mxArray *c5_t_y = NULL;
  real_T c5_t_u;
  const mxArray *c5_u_y = NULL;
  real_T c5_u_u;
  const mxArray *c5_v_y = NULL;
  real_T c5_v_u[300];
  const mxArray *c5_w_y = NULL;
  uint8_T c5_w_u;
  const mxArray *c5_x_y = NULL;
  real_T *c5_Fs1;
  real_T *c5_Fs2;
  real_T *c5_Ws;
  real_T (*c5_Zn)[7];
  real_T (*c5_Wn)[7];
  c5_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_st = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_createcellarray(23), FALSE);
  c5_u = *c5_Fs1;
  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 0, c5_b_y);
  c5_b_u = *c5_Fs2;
  c5_c_y = NULL;
  sf_mex_assign(&c5_c_y, sf_mex_create("y", &c5_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 1, c5_c_y);
  for (c5_i1 = 0; c5_i1 < 7; c5_i1++) {
    c5_c_u[c5_i1] = (*c5_Wn)[c5_i1];
  }

  c5_d_y = NULL;
  sf_mex_assign(&c5_d_y, sf_mex_create("y", c5_c_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c5_y, 2, c5_d_y);
  c5_d_u = *c5_Ws;
  c5_e_y = NULL;
  sf_mex_assign(&c5_e_y, sf_mex_create("y", &c5_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 3, c5_e_y);
  for (c5_i1 = 0; c5_i1 < 7; c5_i1++) {
    c5_e_u[c5_i1] = (*c5_Zn)[c5_i1];
  }

  c5_f_y = NULL;
  sf_mex_assign(&c5_f_y, sf_mex_create("y", c5_e_u, 0, 0U, 1U, 0U, 1, 7), FALSE);
  sf_mex_setcell(c5_y, 4, c5_f_y);
  for (c5_i1 = 0; c5_i1 < 70; c5_i1++) {
    c5_f_u[c5_i1] = chartInstance->c5_ABK[c5_i1];
  }

  c5_g_y = NULL;
  if (!chartInstance->c5_ABK_not_empty) {
    sf_mex_assign(&c5_g_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_g_y, sf_mex_create("y", c5_f_u, 0, 0U, 1U, 0U, 2, 7, 10),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 5, c5_g_y);
  for (c5_i1 = 0; c5_i1 < 16; c5_i1++) {
    c5_g_u[c5_i1] = chartInstance->c5_CD[c5_i1];
  }

  c5_h_y = NULL;
  if (!chartInstance->c5_CD_not_empty) {
    sf_mex_assign(&c5_h_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_h_y, sf_mex_create("y", c5_g_u, 0, 0U, 1U, 0U, 2, 2, 8),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 6, c5_h_y);
  for (c5_i1 = 0; c5_i1 < 2; c5_i1++) {
    c5_h_u[c5_i1] = chartInstance->c5_E1[c5_i1];
  }

  c5_i_y = NULL;
  if (!chartInstance->c5_E1_not_empty) {
    sf_mex_assign(&c5_i_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_i_y, sf_mex_create("y", c5_h_u, 0, 0U, 1U, 0U, 1, 2),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 7, c5_i_y);
  for (c5_i1 = 0; c5_i1 < 25000; c5_i1++) {
    c5_i_u[c5_i1] = chartInstance->c5_LK[c5_i1];
  }

  c5_j_y = NULL;
  if (!chartInstance->c5_LK_not_empty) {
    sf_mex_assign(&c5_j_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_j_y, sf_mex_create("y", c5_i_u, 0, 0U, 1U, 0U, 2, 100, 250),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 8, c5_j_y);
  for (c5_i1 = 0; c5_i1 < 100; c5_i1++) {
    c5_j_u[c5_i1] = chartInstance->c5_Pabk[c5_i1];
  }

  c5_k_y = NULL;
  if (!chartInstance->c5_Pabk_not_empty) {
    sf_mex_assign(&c5_k_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_k_y, sf_mex_create("y", c5_j_u, 0, 0U, 1U, 0U, 2, 10, 10),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 9, c5_k_y);
  for (c5_i1 = 0; c5_i1 < 64; c5_i1++) {
    c5_k_u[c5_i1] = chartInstance->c5_Pcd[c5_i1];
  }

  c5_l_y = NULL;
  if (!chartInstance->c5_Pcd_not_empty) {
    sf_mex_assign(&c5_l_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_l_y, sf_mex_create("y", c5_k_u, 0, 0U, 1U, 0U, 2, 8, 8),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 10, c5_l_y);
  for (c5_i1 = 0; c5_i1 < 63504; c5_i1++) {
    c5_l_u[c5_i1] = chartInstance->c5_Plk[c5_i1];
  }

  c5_m_y = NULL;
  if (!chartInstance->c5_Plk_not_empty) {
    sf_mex_assign(&c5_m_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_m_y, sf_mex_create("y", c5_l_u, 0, 0U, 1U, 0U, 2, 252, 252),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 11, c5_m_y);
  c5_m_u = chartInstance->c5_U1;
  c5_n_y = NULL;
  if (!chartInstance->c5_U1_not_empty) {
    sf_mex_assign(&c5_n_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_n_y, sf_mex_create("y", &c5_m_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c5_y, 12, c5_n_y);
  for (c5_i1 = 0; c5_i1 < 504; c5_i1++) {
    c5_n_u[c5_i1] = chartInstance->c5_VARX[c5_i1];
  }

  c5_o_y = NULL;
  if (!chartInstance->c5_VARX_not_empty) {
    sf_mex_assign(&c5_o_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_o_y, sf_mex_create("y", c5_n_u, 0, 0U, 1U, 0U, 2, 2, 252),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 13, c5_o_y);
  for (c5_i1 = 0; c5_i1 < 7; c5_i1++) {
    c5_o_u[c5_i1] = chartInstance->c5_X[c5_i1];
  }

  c5_p_y = NULL;
  if (!chartInstance->c5_X_not_empty) {
    sf_mex_assign(&c5_p_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_p_y, sf_mex_create("y", c5_o_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 14, c5_p_y);
  for (c5_i1 = 0; c5_i1 < 7; c5_i1++) {
    c5_p_u[c5_i1] = chartInstance->c5_X1[c5_i1];
  }

  c5_q_y = NULL;
  if (!chartInstance->c5_X1_not_empty) {
    sf_mex_assign(&c5_q_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_q_y, sf_mex_create("y", c5_p_u, 0, 0U, 1U, 0U, 1, 7),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 15, c5_q_y);
  for (c5_i1 = 0; c5_i1 < 2; c5_i1++) {
    c5_q_u[c5_i1] = chartInstance->c5_Y1[c5_i1];
  }

  c5_r_y = NULL;
  if (!chartInstance->c5_Y1_not_empty) {
    sf_mex_assign(&c5_r_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_r_y, sf_mex_create("y", c5_q_u, 0, 0U, 1U, 0U, 1, 2),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 16, c5_r_y);
  for (c5_i1 = 0; c5_i1 < 250; c5_i1++) {
    c5_r_u[c5_i1] = chartInstance->c5_Zp[c5_i1];
  }

  c5_s_y = NULL;
  if (!chartInstance->c5_Zp_not_empty) {
    sf_mex_assign(&c5_s_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_s_y, sf_mex_create("y", c5_r_u, 0, 0U, 1U, 0U, 1, 250),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 17, c5_s_y);
  c5_s_u = chartInstance->c5_k;
  c5_t_y = NULL;
  if (!chartInstance->c5_k_not_empty) {
    sf_mex_assign(&c5_t_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_t_y, sf_mex_create("y", &c5_s_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c5_y, 18, c5_t_y);
  c5_t_u = chartInstance->c5_saw;
  c5_u_y = NULL;
  if (!chartInstance->c5_saw_not_empty) {
    sf_mex_assign(&c5_u_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_u_y, sf_mex_create("y", &c5_t_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c5_y, 19, c5_u_y);
  c5_u_u = chartInstance->c5_start;
  c5_v_y = NULL;
  if (!chartInstance->c5_start_not_empty) {
    sf_mex_assign(&c5_v_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_v_y, sf_mex_create("y", &c5_u_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c5_y, 20, c5_v_y);
  for (c5_i1 = 0; c5_i1 < 300; c5_i1++) {
    c5_v_u[c5_i1] = chartInstance->c5_w[c5_i1];
  }

  c5_w_y = NULL;
  if (!chartInstance->c5_w_not_empty) {
    sf_mex_assign(&c5_w_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c5_w_y, sf_mex_create("y", c5_v_u, 0, 0U, 1U, 0U, 2, 1, 300),
                  FALSE);
  }

  sf_mex_setcell(c5_y, 21, c5_w_y);
  c5_w_u = chartInstance->c5_is_active_c5_sim1_rlti_varmax_gust2;
  c5_x_y = NULL;
  sf_mex_assign(&c5_x_y, sf_mex_create("y", &c5_w_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 22, c5_x_y);
  sf_mex_assign(&c5_st, c5_y, FALSE);
  return c5_st;
}

static void set_sim_state_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance, const mxArray
   *c5_st)
{
  const mxArray *c5_u;
  real_T c5_dv2[7];
  int32_T c5_i2;
  real_T c5_dv3[70];
  real_T c5_dv4[16];
  real_T c5_dv5[2];
  static real_T c5_dv6[25000];
  real_T c5_dv7[100];
  real_T c5_dv8[64];
  static real_T c5_dv9[63504];
  real_T c5_dv10[504];
  real_T c5_dv11[250];
  real_T c5_dv12[300];
  real_T *c5_Fs1;
  real_T *c5_Fs2;
  real_T *c5_Ws;
  real_T (*c5_Wn)[7];
  real_T (*c5_Zn)[7];
  c5_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_u = sf_mex_dup(c5_st);
  *c5_Fs1 = c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 0)),
    "Fs1");
  *c5_Fs2 = c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 1)),
    "Fs2");
  c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 2)), "Wn",
                        c5_dv2);
  for (c5_i2 = 0; c5_i2 < 7; c5_i2++) {
    (*c5_Wn)[c5_i2] = c5_dv2[c5_i2];
  }

  *c5_Ws = c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 3)),
    "Ws");
  c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 4)), "Zn",
                        c5_dv2);
  for (c5_i2 = 0; c5_i2 < 7; c5_i2++) {
    (*c5_Zn)[c5_i2] = c5_dv2[c5_i2];
  }

  c5_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 5)),
                        "ABK", c5_dv3);
  for (c5_i2 = 0; c5_i2 < 70; c5_i2++) {
    chartInstance->c5_ABK[c5_i2] = c5_dv3[c5_i2];
  }

  c5_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 6)), "CD",
                        c5_dv4);
  for (c5_i2 = 0; c5_i2 < 16; c5_i2++) {
    chartInstance->c5_CD[c5_i2] = c5_dv4[c5_i2];
  }

  c5_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 7)), "E1",
                        c5_dv5);
  for (c5_i2 = 0; c5_i2 < 2; c5_i2++) {
    chartInstance->c5_E1[c5_i2] = c5_dv5[c5_i2];
  }

  c5_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 8)), "LK",
                        c5_dv6);
  for (c5_i2 = 0; c5_i2 < 25000; c5_i2++) {
    chartInstance->c5_LK[c5_i2] = c5_dv6[c5_i2];
  }

  c5_m_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 9)),
                        "Pabk", c5_dv7);
  for (c5_i2 = 0; c5_i2 < 100; c5_i2++) {
    chartInstance->c5_Pabk[c5_i2] = c5_dv7[c5_i2];
  }

  c5_o_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 10)),
                        "Pcd", c5_dv8);
  for (c5_i2 = 0; c5_i2 < 64; c5_i2++) {
    chartInstance->c5_Pcd[c5_i2] = c5_dv8[c5_i2];
  }

  c5_q_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 11)),
                        "Plk", c5_dv9);
  for (c5_i2 = 0; c5_i2 < 63504; c5_i2++) {
    chartInstance->c5_Plk[c5_i2] = c5_dv9[c5_i2];
  }

  chartInstance->c5_U1 = c5_s_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 12)), "U1");
  c5_u_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 13)),
                        "VARX", c5_dv10);
  for (c5_i2 = 0; c5_i2 < 504; c5_i2++) {
    chartInstance->c5_VARX[c5_i2] = c5_dv10[c5_i2];
  }

  c5_w_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 14)), "X",
                        c5_dv2);
  for (c5_i2 = 0; c5_i2 < 7; c5_i2++) {
    chartInstance->c5_X[c5_i2] = c5_dv2[c5_i2];
  }

  c5_y_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 15)),
                        "X1", c5_dv2);
  for (c5_i2 = 0; c5_i2 < 7; c5_i2++) {
    chartInstance->c5_X1[c5_i2] = c5_dv2[c5_i2];
  }

  c5_bb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 16)),
    "Y1", c5_dv5);
  for (c5_i2 = 0; c5_i2 < 2; c5_i2++) {
    chartInstance->c5_Y1[c5_i2] = c5_dv5[c5_i2];
  }

  c5_db_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 17)),
    "Zp", c5_dv11);
  for (c5_i2 = 0; c5_i2 < 250; c5_i2++) {
    chartInstance->c5_Zp[c5_i2] = c5_dv11[c5_i2];
  }

  chartInstance->c5_k = c5_fb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 18)), "k");
  chartInstance->c5_saw = c5_hb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 19)), "saw");
  chartInstance->c5_start = c5_jb_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 20)), "start");
  c5_lb_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 21)),
    "w", c5_dv12);
  for (c5_i2 = 0; c5_i2 < 300; c5_i2++) {
    chartInstance->c5_w[c5_i2] = c5_dv12[c5_i2];
  }

  chartInstance->c5_is_active_c5_sim1_rlti_varmax_gust2 = c5_nb_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 22)),
     "is_active_c5_sim1_rlti_varmax_gust2");
  sf_mex_destroy(&c5_u);
  sf_mex_destroy(&c5_st);
}

static void finalize_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
}

static void sf_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  c5_chartstep_c5_sim1_rlti_varmax_gust2(chartInstance);
}

static void c5_chartstep_c5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
  real_T c5_ON;
  real_T c5_U;
  int32_T c5_i3;
  real_T c5_Fs[2];
  real_T c5_b_W[700];
  real_T c5_b_Du;
  real_T c5_b_Dy[4];
  real_T c5_h;
  real_T c5_b_saw;
  int32_T c5_ldc;
  int32_T c5_i;
  real_T c5_dv13[300];
  real_T c5_c_Dy[2];
  real_T c5_dv14[2];
  real_T c5_a[4];
  real_T c5_dv15[245];
  real_T c5_Z[252];
  real_T c5_b_VARX[504];
  static real_T c5_b_Plk[63504];
  int32_T c5_n;
  int32_T c5_b_k;
  int32_T c5_lda;
  int32_T c5_ldb;
  char_T c5_TRANSA;
  char_T c5_TRANSB;
  real_T c5_y[1750];
  real_T c5_dv16[8];
  real_T c5_dv17[10];
  real_T c5_A[49];
  real_T c5_b_a[7];
  real_T c5_C[14];
  real_T c5_dv18[49];
  real_T c5_Zn[7];
  real_T c5_Wn[7];
  creal_T c5_b_y;
  creal_T c5_c_y[49];
  real_T c5_c_a[7];
  creal_T c5_b[7];
  creal_T c5_b_C[14];
  creal_T c5_c_C[2];
  real_T *c5_Ws;
  real_T *c5_Fs1;
  real_T *c5_Fs2;
  real_T *c5_b_ON;
  real_T *c5_b_U;
  real_T *c5_RESET;
  real_T (*c5_b_Wn)[7];
  real_T (*c5_b_Zn)[7];
  real_T (*c5_Y)[2];
  c5_b_Zn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_b_Wn = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_Y = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c5_b_U = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c5_Fs2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_Fs1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_Ws = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_RESET = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c5_b_ON = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c5_ON = *c5_b_ON;
  c5_U = *c5_b_U;
  for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
    c5_Fs[c5_i3] = (*c5_Y)[c5_i3];
  }

  for (c5_i3 = 0; c5_i3 < 700; c5_i3++) {
    c5_b_W[c5_i3] = chartInstance->c5_W[c5_i3];
  }

  c5_b_Du = chartInstance->c5_Du;
  for (c5_i3 = 0; c5_i3 < 4; c5_i3++) {
    c5_b_Dy[c5_i3] = chartInstance->c5_Dy[c5_i3];
  }

  c5_h = 0.05;
  if ((!chartInstance->c5_Plk_not_empty) || (*c5_RESET != 0.0)) {
    c5_eye(chartInstance, chartInstance->c5_Plk);
    c5_b_saw = c5_mrdivide(chartInstance, 1.0, 100.0);
    for (c5_i3 = 0; c5_i3 < 63504; c5_i3++) {
      chartInstance->c5_Plk[c5_i3] *= c5_b_saw;
    }

    chartInstance->c5_Plk_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 504; c5_i3++) {
      chartInstance->c5_VARX[c5_i3] = 0.0;
    }

    chartInstance->c5_VARX_not_empty = TRUE;
    chartInstance->c5_saw = 0.0;
    chartInstance->c5_saw_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 25000; c5_i3++) {
      chartInstance->c5_LK[c5_i3] = 0.0;
    }

    chartInstance->c5_LK_not_empty = TRUE;
    c5_b_eye(chartInstance, chartInstance->c5_Pabk);
    c5_b_saw = c5_mrdivide(chartInstance, 1.0, 1.0E-6);
    for (c5_i3 = 0; c5_i3 < 100; c5_i3++) {
      chartInstance->c5_Pabk[c5_i3] *= c5_b_saw;
    }

    chartInstance->c5_Pabk_not_empty = TRUE;
    c5_i3 = 0;
    for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
      for (c5_i = 0; c5_i < 7; c5_i++) {
        chartInstance->c5_ABK[c5_i + c5_i3] = 0.0;
      }

      c5_i3 += 7;
    }

    for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
      chartInstance->c5_ABK[c5_i3 + 49] = 0.0;
    }

    c5_i3 = 0;
    for (c5_ldc = 0; c5_ldc < 2; c5_ldc++) {
      for (c5_i = 0; c5_i < 7; c5_i++) {
        chartInstance->c5_ABK[(c5_i + c5_i3) + 56] = 0.0;
      }

      c5_i3 += 7;
    }

    chartInstance->c5_ABK_not_empty = TRUE;
    c5_c_eye(chartInstance, chartInstance->c5_Pcd);
    c5_b_saw = c5_mrdivide(chartInstance, 1.0, 1.0E-6);
    for (c5_i3 = 0; c5_i3 < 64; c5_i3++) {
      chartInstance->c5_Pcd[c5_i3] *= c5_b_saw;
    }

    chartInstance->c5_Pcd_not_empty = TRUE;
    c5_i3 = 0;
    for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
      for (c5_i = 0; c5_i < 2; c5_i++) {
        chartInstance->c5_CD[c5_i + c5_i3] = 0.0;
      }

      c5_i3 += 2;
    }

    chartInstance->c5_CD_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      chartInstance->c5_CD[c5_i3 + 14] = 0.0;
      chartInstance->c5_Y1[c5_i3] = 0.0;
    }

    chartInstance->c5_Y1_not_empty = TRUE;
    chartInstance->c5_U1 = 0.0;
    chartInstance->c5_U1_not_empty = TRUE;
    chartInstance->c5_X1_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
      chartInstance->c5_X1[c5_i3] = 0.0;
      chartInstance->c5_X[c5_i3] = 0.0;
    }

    chartInstance->c5_X_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      chartInstance->c5_E1[c5_i3] = 0.0;
    }

    chartInstance->c5_E1_not_empty = TRUE;
    for (c5_i3 = 0; c5_i3 < 250; c5_i3++) {
      chartInstance->c5_Zp[c5_i3] = 0.0;
    }

    chartInstance->c5_Zp_not_empty = TRUE;
    c5_logspace(chartInstance, -1.0, 1.0, c5_dv13);
    for (c5_i3 = 0; c5_i3 < 300; c5_i3++) {
      chartInstance->c5_w[c5_i3] = c5_dv13[c5_i3];
    }

    chartInstance->c5_w_not_empty = TRUE;
    chartInstance->c5_k = 1.0;
    chartInstance->c5_k_not_empty = TRUE;
    chartInstance->c5_start = 1.0;
    chartInstance->c5_start_not_empty = TRUE;
  }

  if (c5_ON > 0.5) {
    c5_U *= c5_rdivide(chartInstance, 1.0, c5_b_Du);
    c5_i3 = 0;
    for (c5_ldc = 0; c5_ldc < 2; c5_ldc++) {
      c5_c_Dy[c5_ldc] = c5_b_Dy[c5_i3];
      c5_i3 += 3;
    }

    c5_b_rdivide(chartInstance, 1.0, c5_c_Dy, c5_dv14);
    c5_diag(chartInstance, c5_dv14, c5_a);
    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      c5_c_Dy[c5_i3] = 0.0;
      c5_ldc = 0;
      for (c5_i = 0; c5_i < 2; c5_i++) {
        c5_c_Dy[c5_i3] += c5_a[c5_ldc + c5_i3] * c5_Fs[c5_i];
        c5_ldc += 2;
      }
    }

    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      c5_Fs[c5_i3] = c5_c_Dy[c5_i3];
    }

    for (c5_i3 = 0; c5_i3 < 245; c5_i3++) {
      c5_dv15[c5_i3] = chartInstance->c5_Zp[c5_i3 + 5];
    }

    for (c5_i3 = 0; c5_i3 < 245; c5_i3++) {
      chartInstance->c5_Zp[c5_i3] = c5_dv15[c5_i3];
    }

    chartInstance->c5_Zp[245] = chartInstance->c5_U1;
    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      chartInstance->c5_Zp[c5_i3 + 246] = chartInstance->c5_Y1[c5_i3];
    }

    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      chartInstance->c5_Zp[c5_i3 + 248] = chartInstance->c5_E1[c5_i3];
    }

    for (c5_i3 = 0; c5_i3 < 250; c5_i3++) {
      c5_Z[c5_i3] = chartInstance->c5_Zp[c5_i3];
    }

    c5_Z[250] = c5_U;
    c5_Z[251] = 1.0;
    if ((chartInstance->c5_start < 0.5) || (chartInstance->c5_k >= 50.0)) {
      for (c5_i3 = 0; c5_i3 < 504; c5_i3++) {
        c5_b_VARX[c5_i3] = chartInstance->c5_VARX[c5_i3];
      }

      for (c5_i3 = 0; c5_i3 < 63504; c5_i3++) {
        c5_b_Plk[c5_i3] = chartInstance->c5_Plk[c5_i3];
      }

      c5_b_saw = chartInstance->c5_saw;
      c5_rls_ew_reg(chartInstance, c5_Z, c5_Fs, c5_b_VARX, c5_b_Plk, 0.999, 1.0,
                    &c5_b_saw);
      for (c5_i3 = 0; c5_i3 < 504; c5_i3++) {
        chartInstance->c5_VARX[c5_i3] = c5_b_VARX[c5_i3];
      }

      for (c5_i3 = 0; c5_i3 < 63504; c5_i3++) {
        chartInstance->c5_Plk[c5_i3] = c5_b_Plk[c5_i3];
      }

      chartInstance->c5_saw = c5_b_saw;
      for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
        c5_b_saw = 0.0;
        c5_ldc = 0;
        for (c5_i = 0; c5_i < 252; c5_i++) {
          c5_b_saw += chartInstance->c5_VARX[c5_ldc + c5_i3] * c5_Z[c5_i];
          c5_ldc += 2;
        }

        chartInstance->c5_E1[c5_i3] = c5_Fs[c5_i3] - c5_b_saw;
      }
    }

    if ((chartInstance->c5_start < 0.5) || (chartInstance->c5_k >= 150.0)) {
      for (c5_i3 = 0; c5_i3 < 25000; c5_i3++) {
        chartInstance->c5_LK[c5_i3] = 0.0;
      }

      for (c5_i = 0; c5_i < 50; c5_i++) {
        c5_h = 0.0;
        c5_n = 0;
        while (c5_n <= 49 - c5_i) {
          c5_h = (real_T)c5_n;
          c5_b_k = (c5_i << 1) - 1;
          c5_lda = (c5_i + c5_n) * 5;
          c5_ldb = c5_n * 5 - 1;
          for (c5_i3 = 0; c5_i3 < 5; c5_i3++) {
            for (c5_ldc = 0; c5_ldc < 2; c5_ldc++) {
              chartInstance->c5_LK[((c5_ldc + c5_b_k) + 100 *
                                    (sf_mex_lw_bounds_check((c5_i3 + c5_lda) + 1,
                1, 250) - 1)) + 1] = chartInstance->c5_VARX[c5_ldc + (((c5_i3 +
                c5_ldb) + 1) << 1)];
            }
          }

          c5_n++;
          sf_mex_listen_for_ctrl_c(chartInstance->S);
        }

        sf_mex_listen_for_ctrl_c(chartInstance->S);
      }

      c5_i = 7;
      c5_n = 250;
      c5_b_k = 100;
      c5_b_saw = 1.0;
      c5_lda = 7;
      c5_ldb = 100;
      c5_ON = 0.0;
      c5_ldc = 7;
      c5_TRANSA = 'N';
      c5_TRANSB = 'N';
      for (c5_i3 = 0; c5_i3 < 1750; c5_i3++) {
        c5_y[c5_i3] = 0.0;
      }

      dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_i, &c5_n, &c5_b_k, &c5_b_saw, &c5_b_W
              [0], &c5_lda, &chartInstance->c5_LK[0], &c5_ldb, &c5_ON, &c5_y[0],
              &c5_ldc);
      c5_i = 7;
      c5_n = 1;
      c5_b_k = 250;
      c5_b_saw = 1.0;
      c5_lda = 7;
      c5_ldb = 250;
      c5_ON = 0.0;
      c5_ldc = 7;
      c5_TRANSA = 'N';
      c5_TRANSB = 'N';
      for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
        chartInstance->c5_X[c5_i3] = 0.0;
      }

      dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_i, &c5_n, &c5_b_k, &c5_b_saw, &c5_y[0],
              &c5_lda, &chartInstance->c5_Zp[0], &c5_ldb, &c5_ON,
              &chartInstance->c5_X[0], &c5_ldc);
    }

    if ((chartInstance->c5_start < 0.5) || (chartInstance->c5_k >= 151.0)) {
      for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
        c5_dv16[c5_i3] = chartInstance->c5_X1[c5_i3];
      }

      c5_dv16[7] = chartInstance->c5_U1;
      c5_rls_ew(chartInstance, c5_dv16, chartInstance->c5_Y1,
                chartInstance->c5_CD, chartInstance->c5_Pcd, 0.999);
      for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
        c5_dv16[c5_i3] = chartInstance->c5_X1[c5_i3];
      }

      c5_dv16[7] = chartInstance->c5_U1;
      for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
        c5_b_saw = 0.0;
        c5_ldc = 0;
        for (c5_i = 0; c5_i < 8; c5_i++) {
          c5_b_saw += chartInstance->c5_CD[c5_ldc + c5_i3] * c5_dv16[c5_i];
          c5_ldc += 2;
        }

        c5_c_Dy[c5_i3] = chartInstance->c5_Y1[c5_i3] - c5_b_saw;
      }

      for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
        c5_dv17[c5_i3] = chartInstance->c5_X1[c5_i3];
      }

      c5_dv17[7] = chartInstance->c5_U1;
      for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
        c5_dv17[c5_i3 + 8] = c5_c_Dy[c5_i3];
      }

      c5_b_rls_ew(chartInstance, c5_dv17, chartInstance->c5_X,
                  chartInstance->c5_ABK, chartInstance->c5_Pabk, 0.999);
    }

    if ((chartInstance->c5_start < 0.5) || (chartInstance->c5_k >= 150.0)) {
      for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
        chartInstance->c5_X1[c5_i3] = chartInstance->c5_X[c5_i3];
      }
    }

    for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
      chartInstance->c5_Y1[c5_i3] = c5_Fs[c5_i3];
    }

    chartInstance->c5_U1 = c5_U;
  }

  c5_b_saw = c5_rdivide(chartInstance, 1.0, c5_b_Du);
  c5_i3 = 0;
  for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_A[c5_i + c5_i3] = chartInstance->c5_ABK[c5_i + c5_i3];
    }

    c5_b_a[c5_ldc] = chartInstance->c5_ABK[c5_ldc + 49];
    c5_i3 += 7;
  }

  for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
    for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
      c5_C[c5_i3 + (c5_ldc << 1)] = 0.0;
      for (c5_i = 0; c5_i < 2; c5_i++) {
        c5_C[c5_i3 + (c5_ldc << 1)] += c5_b_Dy[c5_i3 + (c5_i << 1)] *
          chartInstance->c5_CD[c5_i + (c5_ldc << 1)];
      }
    }

    c5_Fs[c5_i3] = 0.0;
    for (c5_ldc = 0; c5_ldc < 2; c5_ldc++) {
      c5_Fs[c5_i3] += c5_b_Dy[c5_i3 + (c5_ldc << 1)] * chartInstance->c5_CD[14 +
        c5_ldc];
    }
  }

  c5_U = c5_rdivide(chartInstance, 1.0, c5_b_Du);
  c5_i3 = 0;
  for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_dv18[c5_i + c5_i3] = chartInstance->c5_ABK[c5_i + c5_i3];
    }

    c5_i3 += 7;
  }

  c5_damp(chartInstance, c5_dv18, c5_h, c5_Wn, c5_Zn);
  c5_b_Du = chartInstance->c5_w[sf_mex_lw_bounds_check((int32_T)
    chartInstance->c5_k, 1, 300) - 1];
  c5_b_y.re = 3.1415926535897931 * (2.0 * (chartInstance->c5_w[(int32_T)
    chartInstance->c5_k - 1] * 0.0));
  c5_b_y.im = 3.1415926535897931 * (2.0 * (chartInstance->c5_w[(int32_T)
    chartInstance->c5_k - 1] * c5_h));
  c5_exp(chartInstance, &c5_b_y);
  c5_d_eye(chartInstance, c5_dv18);
  for (c5_i3 = 0; c5_i3 < 49; c5_i3++) {
    c5_c_y[c5_i3].re = c5_b_y.re * c5_dv18[c5_i3] - c5_A[c5_i3];
    c5_c_y[c5_i3].im = c5_b_y.im * c5_dv18[c5_i3];
  }

  for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
    c5_c_a[c5_i3] = c5_b_a[c5_i3] * c5_b_saw;
  }

  c5_mldivide(chartInstance, c5_c_y, c5_c_a, c5_b);
  c5_i3 = 0;
  for (c5_ldc = 0; c5_ldc < 7; c5_ldc++) {
    for (c5_i = 0; c5_i < 2; c5_i++) {
      c5_b_C[c5_i + c5_i3].re = c5_C[c5_i + c5_i3];
      c5_b_C[c5_i + c5_i3].im = 0.0;
    }

    c5_i3 += 2;
  }

  for (c5_i3 = 0; c5_i3 < 2; c5_i3++) {
    c5_ON = 0.0;
    c5_b_saw = 0.0;
    c5_ldc = 0;
    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_ON += c5_b_C[c5_ldc + c5_i3].re * c5_b[c5_i].re - 0.0 * c5_b[c5_i].im;
      c5_b_saw += c5_b_C[c5_ldc + c5_i3].re * c5_b[c5_i].im + 0.0 * c5_b[c5_i].
        re;
      c5_ldc += 2;
    }

    c5_c_C[c5_i3].re = c5_ON + c5_Fs[c5_i3] * c5_U;
    c5_c_C[c5_i3].im = c5_b_saw;
  }

  c5_abs(chartInstance, c5_c_C, c5_Fs);
  c5_log10(chartInstance, &c5_b_Du);
  if (c5_Fs[0] < 1.0E-5) {
    c5_ON = -100.0;
  } else {
    c5_b_saw = c5_Fs[0];
    c5_log10(chartInstance, &c5_b_saw);
    c5_ON = 20.0 * c5_b_saw;
  }

  if (c5_Fs[1] < 1.0E-5) {
    *c5_Fs2 = -100.0;
  } else {
    c5_b_saw = c5_Fs[1];
    c5_log10(chartInstance, &c5_b_saw);
    *c5_Fs2 = 20.0 * c5_b_saw;
  }

  if (chartInstance->c5_k >= 300.0) {
    chartInstance->c5_start = 0.0;
    chartInstance->c5_k = 1.0;
  } else {
    chartInstance->c5_k++;
  }

  *c5_Ws = c5_b_Du;
  *c5_Fs1 = c5_ON;
  for (c5_i3 = 0; c5_i3 < 7; c5_i3++) {
    (*c5_b_Wn)[c5_i3] = c5_Wn[c5_i3];
    (*c5_b_Zn)[c5_i3] = c5_Zn[c5_i3];
  }
}

static void initSimStructsc5_sim1_rlti_varmax_gust2
  (SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber)
{
}

const mxArray *sf_c5_sim1_rlti_varmax_gust2_get_eml_resolved_functions_info(void)
{
  const mxArray *c5_nameCaptureInfo;
  c5_ResolvedFunctionInfo c5_info[106];
  const mxArray *c5_m0 = NULL;
  int32_T c5_i4;
  c5_ResolvedFunctionInfo *c5_r0;
  c5_nameCaptureInfo = NULL;
  c5_info_helper(c5_info);
  c5_b_info_helper(c5_info);
  sf_mex_assign(&c5_m0, sf_mex_createstruct("nameCaptureInfo", 1, 106), FALSE);
  for (c5_i4 = 0; c5_i4 < 106; c5_i4++) {
    c5_r0 = &c5_info[c5_i4];
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c5_r0->context)), "context", "nameCaptureInfo",
                    c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c5_r0->name)), "name", "nameCaptureInfo", c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c5_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c5_r0->resolved)), "resolved", "nameCaptureInfo",
                    c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c5_i4);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c5_i4);
  }

  sf_mex_assign(&c5_nameCaptureInfo, c5_m0, FALSE);
  return c5_nameCaptureInfo;
}

static void c5_info_helper(c5_ResolvedFunctionInfo c5_info[106])
{
  c5_info[0].context = "";
  c5_info[0].name = "mtimes";
  c5_info[0].dominantType = "double";
  c5_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c5_info[0].fileTimeLo = 1289519692U;
  c5_info[0].fileTimeHi = 0U;
  c5_info[0].mFileTimeLo = 0U;
  c5_info[0].mFileTimeHi = 0U;
  c5_info[1].context = "";
  c5_info[1].name = "mrdivide";
  c5_info[1].dominantType = "double";
  c5_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c5_info[1].fileTimeLo = 1310137456U;
  c5_info[1].fileTimeHi = 0U;
  c5_info[1].mFileTimeLo = 1289519692U;
  c5_info[1].mFileTimeHi = 0U;
  c5_info[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c5_info[2].name = "rdivide";
  c5_info[2].dominantType = "double";
  c5_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c5_info[2].fileTimeLo = 1286818844U;
  c5_info[2].fileTimeHi = 0U;
  c5_info[2].mFileTimeLo = 0U;
  c5_info[2].mFileTimeHi = 0U;
  c5_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c5_info[3].name = "eml_div";
  c5_info[3].dominantType = "double";
  c5_info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c5_info[3].fileTimeLo = 1305318000U;
  c5_info[3].fileTimeHi = 0U;
  c5_info[3].mFileTimeLo = 0U;
  c5_info[3].mFileTimeHi = 0U;
  c5_info[4].context = "";
  c5_info[4].name = "eye";
  c5_info[4].dominantType = "double";
  c5_info[4].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m";
  c5_info[4].fileTimeLo = 1286818688U;
  c5_info[4].fileTimeHi = 0U;
  c5_info[4].mFileTimeLo = 0U;
  c5_info[4].mFileTimeHi = 0U;
  c5_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c5_info[5].name = "eml_assert_valid_size_arg";
  c5_info[5].dominantType = "double";
  c5_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c5_info[5].fileTimeLo = 1286818694U;
  c5_info[5].fileTimeHi = 0U;
  c5_info[5].mFileTimeLo = 0U;
  c5_info[5].mFileTimeHi = 0U;
  c5_info[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral";
  c5_info[6].name = "isinf";
  c5_info[6].dominantType = "double";
  c5_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m";
  c5_info[6].fileTimeLo = 1286818760U;
  c5_info[6].fileTimeHi = 0U;
  c5_info[6].mFileTimeLo = 0U;
  c5_info[6].mFileTimeHi = 0U;
  c5_info[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c5_info[7].name = "eml_index_class";
  c5_info[7].dominantType = "";
  c5_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c5_info[7].fileTimeLo = 1286818778U;
  c5_info[7].fileTimeHi = 0U;
  c5_info[7].mFileTimeLo = 0U;
  c5_info[7].mFileTimeHi = 0U;
  c5_info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c5_info[8].name = "intmax";
  c5_info[8].dominantType = "char";
  c5_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c5_info[8].fileTimeLo = 1286818756U;
  c5_info[8].fileTimeHi = 0U;
  c5_info[8].mFileTimeLo = 0U;
  c5_info[8].mFileTimeHi = 0U;
  c5_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c5_info[9].name = "eml_is_float_class";
  c5_info[9].dominantType = "char";
  c5_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c5_info[9].fileTimeLo = 1286818782U;
  c5_info[9].fileTimeHi = 0U;
  c5_info[9].mFileTimeLo = 0U;
  c5_info[9].mFileTimeHi = 0U;
  c5_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c5_info[10].name = "min";
  c5_info[10].dominantType = "double";
  c5_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c5_info[10].fileTimeLo = 1308747330U;
  c5_info[10].fileTimeHi = 0U;
  c5_info[10].mFileTimeLo = 0U;
  c5_info[10].mFileTimeHi = 0U;
  c5_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c5_info[11].name = "eml_min_or_max";
  c5_info[11].dominantType = "char";
  c5_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  c5_info[11].fileTimeLo = 1303146212U;
  c5_info[11].fileTimeHi = 0U;
  c5_info[11].mFileTimeLo = 0U;
  c5_info[11].mFileTimeHi = 0U;
  c5_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c5_info[12].name = "eml_scalar_eg";
  c5_info[12].dominantType = "double";
  c5_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c5_info[12].fileTimeLo = 1286818796U;
  c5_info[12].fileTimeHi = 0U;
  c5_info[12].mFileTimeLo = 0U;
  c5_info[12].mFileTimeHi = 0U;
  c5_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c5_info[13].name = "eml_scalexp_alloc";
  c5_info[13].dominantType = "double";
  c5_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c5_info[13].fileTimeLo = 1286818796U;
  c5_info[13].fileTimeHi = 0U;
  c5_info[13].mFileTimeLo = 0U;
  c5_info[13].mFileTimeHi = 0U;
  c5_info[14].context = "";
  c5_info[14].name = "logspace";
  c5_info[14].dominantType = "double";
  c5_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[14].fileTimeLo = 1305318000U;
  c5_info[14].fileTimeHi = 0U;
  c5_info[14].mFileTimeLo = 0U;
  c5_info[14].mFileTimeHi = 0U;
  c5_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[15].name = "eml_scalar_floor";
  c5_info[15].dominantType = "double";
  c5_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c5_info[15].fileTimeLo = 1286818726U;
  c5_info[15].fileTimeHi = 0U;
  c5_info[15].mFileTimeLo = 0U;
  c5_info[15].mFileTimeHi = 0U;
  c5_info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[16].name = "abs";
  c5_info[16].dominantType = "double";
  c5_info[16].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c5_info[16].fileTimeLo = 1286818694U;
  c5_info[16].fileTimeHi = 0U;
  c5_info[16].mFileTimeLo = 0U;
  c5_info[16].mFileTimeHi = 0U;
  c5_info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c5_info[17].name = "eml_scalar_abs";
  c5_info[17].dominantType = "double";
  c5_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c5_info[17].fileTimeLo = 1286818712U;
  c5_info[17].fileTimeHi = 0U;
  c5_info[17].mFileTimeLo = 0U;
  c5_info[17].mFileTimeHi = 0U;
  c5_info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[18].name = "eps";
  c5_info[18].dominantType = "char";
  c5_info[18].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c5_info[18].fileTimeLo = 1286818686U;
  c5_info[18].fileTimeHi = 0U;
  c5_info[18].mFileTimeLo = 0U;
  c5_info[18].mFileTimeHi = 0U;
  c5_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[19].name = "eml_warning";
  c5_info[19].dominantType = "char";
  c5_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c5_info[19].fileTimeLo = 1286818802U;
  c5_info[19].fileTimeHi = 0U;
  c5_info[19].mFileTimeLo = 0U;
  c5_info[19].mFileTimeHi = 0U;
  c5_info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/logspace.m";
  c5_info[20].name = "linspace";
  c5_info[20].dominantType = "double";
  c5_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c5_info[20].fileTimeLo = 1286818762U;
  c5_info[20].fileTimeHi = 0U;
  c5_info[20].mFileTimeLo = 0U;
  c5_info[20].mFileTimeHi = 0U;
  c5_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/linspace.m";
  c5_info[21].name = "realmax";
  c5_info[21].dominantType = "char";
  c5_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c5_info[21].fileTimeLo = 1286818766U;
  c5_info[21].fileTimeHi = 0U;
  c5_info[21].mFileTimeLo = 0U;
  c5_info[21].mFileTimeHi = 0U;
  c5_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m";
  c5_info[22].name = "mpower";
  c5_info[22].dominantType = "double";
  c5_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c5_info[22].fileTimeLo = 1286818842U;
  c5_info[22].fileTimeHi = 0U;
  c5_info[22].mFileTimeLo = 0U;
  c5_info[22].mFileTimeHi = 0U;
  c5_info[23].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c5_info[23].name = "power";
  c5_info[23].dominantType = "double";
  c5_info[23].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  c5_info[23].fileTimeLo = 1294067944U;
  c5_info[23].fileTimeHi = 0U;
  c5_info[23].mFileTimeLo = 0U;
  c5_info[23].mFileTimeHi = 0U;
  c5_info[24].context = "";
  c5_info[24].name = "diag";
  c5_info[24].dominantType = "double";
  c5_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c5_info[24].fileTimeLo = 1286818686U;
  c5_info[24].fileTimeHi = 0U;
  c5_info[24].mFileTimeLo = 0U;
  c5_info[24].mFileTimeHi = 0U;
  c5_info[25].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c5_info[25].name = "eml_index_plus";
  c5_info[25].dominantType = "int32";
  c5_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c5_info[25].fileTimeLo = 1286818778U;
  c5_info[25].fileTimeHi = 0U;
  c5_info[25].mFileTimeLo = 0U;
  c5_info[25].mFileTimeHi = 0U;
  c5_info[26].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c5_info[26].name = "eml_index_times";
  c5_info[26].dominantType = "int32";
  c5_info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c5_info[26].fileTimeLo = 1286818780U;
  c5_info[26].fileTimeHi = 0U;
  c5_info[26].mFileTimeLo = 0U;
  c5_info[26].mFileTimeHi = 0U;
  c5_info[27].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m";
  c5_info[27].name = "eml_index_minus";
  c5_info[27].dominantType = "int32";
  c5_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c5_info[27].fileTimeLo = 1286818778U;
  c5_info[27].fileTimeHi = 0U;
  c5_info[27].mFileTimeLo = 0U;
  c5_info[27].mFileTimeHi = 0U;
  c5_info[28].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c5_info[28].name = "eml_xgemm";
  c5_info[28].dominantType = "int32";
  c5_info[28].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c5_info[28].fileTimeLo = 1299076772U;
  c5_info[28].fileTimeHi = 0U;
  c5_info[28].mFileTimeLo = 0U;
  c5_info[28].mFileTimeHi = 0U;
  c5_info[29].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c5_info[29].name = "eml_blas_inline";
  c5_info[29].dominantType = "";
  c5_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c5_info[29].fileTimeLo = 1299076768U;
  c5_info[29].fileTimeHi = 0U;
  c5_info[29].mFileTimeLo = 0U;
  c5_info[29].mFileTimeHi = 0U;
  c5_info[30].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c5_info[30].name = "eml_refblas_xgemm";
  c5_info[30].dominantType = "int32";
  c5_info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c5_info[30].fileTimeLo = 1299076774U;
  c5_info[30].fileTimeHi = 0U;
  c5_info[30].mFileTimeLo = 0U;
  c5_info[30].mFileTimeHi = 0U;
  c5_info[31].context = "";
  c5_info[31].name = "sqrt";
  c5_info[31].dominantType = "double";
  c5_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[31].fileTimeLo = 1286818752U;
  c5_info[31].fileTimeHi = 0U;
  c5_info[31].mFileTimeLo = 0U;
  c5_info[31].mFileTimeHi = 0U;
  c5_info[32].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[32].name = "eml_error";
  c5_info[32].dominantType = "char";
  c5_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c5_info[32].fileTimeLo = 1305318000U;
  c5_info[32].fileTimeHi = 0U;
  c5_info[32].mFileTimeLo = 0U;
  c5_info[32].mFileTimeHi = 0U;
  c5_info[33].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[33].name = "eml_scalar_sqrt";
  c5_info[33].dominantType = "double";
  c5_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c5_info[33].fileTimeLo = 1286818738U;
  c5_info[33].fileTimeHi = 0U;
  c5_info[33].mFileTimeLo = 0U;
  c5_info[33].mFileTimeHi = 0U;
  c5_info[34].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c5_info[34].name = "eml_xdotu";
  c5_info[34].dominantType = "int32";
  c5_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c5_info[34].fileTimeLo = 1299076772U;
  c5_info[34].fileTimeHi = 0U;
  c5_info[34].mFileTimeLo = 0U;
  c5_info[34].mFileTimeHi = 0U;
  c5_info[35].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c5_info[35].name = "eml_xdot";
  c5_info[35].dominantType = "int32";
  c5_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c5_info[35].fileTimeLo = 1299076772U;
  c5_info[35].fileTimeHi = 0U;
  c5_info[35].mFileTimeLo = 0U;
  c5_info[35].mFileTimeHi = 0U;
  c5_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c5_info[36].name = "eml_refblas_xdot";
  c5_info[36].dominantType = "int32";
  c5_info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c5_info[36].fileTimeLo = 1299076772U;
  c5_info[36].fileTimeHi = 0U;
  c5_info[36].mFileTimeLo = 0U;
  c5_info[36].mFileTimeHi = 0U;
  c5_info[37].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c5_info[37].name = "eml_refblas_xdotx";
  c5_info[37].dominantType = "int32";
  c5_info[37].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c5_info[37].fileTimeLo = 1299076774U;
  c5_info[37].fileTimeHi = 0U;
  c5_info[37].mFileTimeLo = 0U;
  c5_info[37].mFileTimeHi = 0U;
  c5_info[38].context = "";
  c5_info[38].name = "norm";
  c5_info[38].dominantType = "double";
  c5_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m";
  c5_info[38].fileTimeLo = 1286818826U;
  c5_info[38].fileTimeHi = 0U;
  c5_info[38].mFileTimeLo = 0U;
  c5_info[38].mFileTimeHi = 0U;
  c5_info[39].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm";
  c5_info[39].name = "eml_xnrm2";
  c5_info[39].dominantType = "int32";
  c5_info[39].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m";
  c5_info[39].fileTimeLo = 1299076776U;
  c5_info[39].fileTimeHi = 0U;
  c5_info[39].mFileTimeLo = 0U;
  c5_info[39].mFileTimeHi = 0U;
  c5_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m";
  c5_info[40].name = "eml_refblas_xnrm2";
  c5_info[40].dominantType = "int32";
  c5_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c5_info[40].fileTimeLo = 1299076784U;
  c5_info[40].fileTimeHi = 0U;
  c5_info[40].mFileTimeLo = 0U;
  c5_info[40].mFileTimeHi = 0U;
  c5_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m";
  c5_info[41].name = "realmin";
  c5_info[41].dominantType = "char";
  c5_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c5_info[41].fileTimeLo = 1286818768U;
  c5_info[41].fileTimeHi = 0U;
  c5_info[41].mFileTimeLo = 0U;
  c5_info[41].mFileTimeHi = 0U;
  c5_info[42].context = "";
  c5_info[42].name = "colon";
  c5_info[42].dominantType = "double";
  c5_info[42].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c5_info[42].fileTimeLo = 1286818838U;
  c5_info[42].fileTimeHi = 0U;
  c5_info[42].mFileTimeLo = 0U;
  c5_info[42].mFileTimeHi = 0U;
  c5_info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c5_info[43].name = "isfinite";
  c5_info[43].dominantType = "double";
  c5_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c5_info[43].fileTimeLo = 1286818758U;
  c5_info[43].fileTimeHi = 0U;
  c5_info[43].mFileTimeLo = 0U;
  c5_info[43].mFileTimeHi = 0U;
  c5_info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m";
  c5_info[44].name = "isnan";
  c5_info[44].dominantType = "double";
  c5_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  c5_info[44].fileTimeLo = 1286818760U;
  c5_info[44].fileTimeHi = 0U;
  c5_info[44].mFileTimeLo = 0U;
  c5_info[44].mFileTimeHi = 0U;
  c5_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!is_flint_colon";
  c5_info[45].name = "floor";
  c5_info[45].dominantType = "double";
  c5_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c5_info[45].fileTimeLo = 1286818742U;
  c5_info[45].fileTimeHi = 0U;
  c5_info[45].mFileTimeLo = 0U;
  c5_info[45].mFileTimeHi = 0U;
  c5_info[46].context = "";
  c5_info[46].name = "eig";
  c5_info[46].dominantType = "double";
  c5_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c5_info[46].fileTimeLo = 1305318000U;
  c5_info[46].fileTimeHi = 0U;
  c5_info[46].mFileTimeLo = 0U;
  c5_info[46].mFileTimeHi = 0U;
  c5_info[47].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m";
  c5_info[47].name = "eml_xgeev";
  c5_info[47].dominantType = "double";
  c5_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c5_info[47].fileTimeLo = 1286818804U;
  c5_info[47].fileTimeHi = 0U;
  c5_info[47].mFileTimeLo = 0U;
  c5_info[47].mFileTimeHi = 0U;
  c5_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m";
  c5_info[48].name = "eml_lapack_xgeev";
  c5_info[48].dominantType = "double";
  c5_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c5_info[48].fileTimeLo = 1301328468U;
  c5_info[48].fileTimeHi = 0U;
  c5_info[48].mFileTimeLo = 0U;
  c5_info[48].mFileTimeHi = 0U;
  c5_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m";
  c5_info[49].name = "eml_matlab_zggev";
  c5_info[49].dominantType = "double";
  c5_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[49].fileTimeLo = 1286818818U;
  c5_info[49].fileTimeHi = 0U;
  c5_info[49].mFileTimeLo = 0U;
  c5_info[49].mFileTimeHi = 0U;
  c5_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[50].name = "eml_matlab_zlangeM";
  c5_info[50].dominantType = "double";
  c5_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c5_info[50].fileTimeLo = 1286818820U;
  c5_info[50].fileTimeHi = 0U;
  c5_info[50].mFileTimeLo = 0U;
  c5_info[50].mFileTimeHi = 0U;
  c5_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c5_info[51].name = "eml_dlapy2";
  c5_info[51].dominantType = "double";
  c5_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m";
  c5_info[51].fileTimeLo = 1286818698U;
  c5_info[51].fileTimeHi = 0U;
  c5_info[51].mFileTimeLo = 0U;
  c5_info[51].mFileTimeHi = 0U;
  c5_info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m";
  c5_info[52].name = "eml_guarded_nan";
  c5_info[52].dominantType = "char";
  c5_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m";
  c5_info[52].fileTimeLo = 1286818776U;
  c5_info[52].fileTimeHi = 0U;
  c5_info[52].mFileTimeLo = 0U;
  c5_info[52].mFileTimeHi = 0U;
  c5_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[53].name = "eml_matlab_zlascl";
  c5_info[53].dominantType = "double";
  c5_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m";
  c5_info[53].fileTimeLo = 1286818822U;
  c5_info[53].fileTimeHi = 0U;
  c5_info[53].mFileTimeLo = 0U;
  c5_info[53].mFileTimeHi = 0U;
  c5_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[54].name = "eml_matlab_zggbal";
  c5_info[54].dominantType = "double";
  c5_info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m";
  c5_info[54].fileTimeLo = 1286818818U;
  c5_info[54].fileTimeHi = 0U;
  c5_info[54].mFileTimeLo = 0U;
  c5_info[54].mFileTimeHi = 0U;
  c5_info[55].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[55].name = "eml_matlab_zgghrd";
  c5_info[55].dominantType = "int32";
  c5_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c5_info[55].fileTimeLo = 1286818820U;
  c5_info[55].fileTimeHi = 0U;
  c5_info[55].mFileTimeLo = 0U;
  c5_info[55].mFileTimeHi = 0U;
  c5_info[56].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c5_info[56].name = "eml_matlab_zlartg";
  c5_info[56].dominantType = "double";
  c5_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c5_info[56].fileTimeLo = 1286818822U;
  c5_info[56].fileTimeHi = 0U;
  c5_info[56].mFileTimeLo = 0U;
  c5_info[56].mFileTimeHi = 0U;
  c5_info[57].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m";
  c5_info[57].name = "fix";
  c5_info[57].dominantType = "double";
  c5_info[57].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c5_info[57].fileTimeLo = 1286818742U;
  c5_info[57].fileTimeHi = 0U;
  c5_info[57].mFileTimeLo = 0U;
  c5_info[57].mFileTimeHi = 0U;
  c5_info[58].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m";
  c5_info[58].name = "eml_scalar_fix";
  c5_info[58].dominantType = "double";
  c5_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m";
  c5_info[58].fileTimeLo = 1286818726U;
  c5_info[58].fileTimeHi = 0U;
  c5_info[58].mFileTimeLo = 0U;
  c5_info[58].mFileTimeHi = 0U;
  c5_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c5_info[59].name = "eml_zrot_rows";
  c5_info[59].dominantType = "int32";
  c5_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c5_info[59].fileTimeLo = 1286818826U;
  c5_info[59].fileTimeHi = 0U;
  c5_info[59].mFileTimeLo = 0U;
  c5_info[59].mFileTimeHi = 0U;
  c5_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m";
  c5_info[60].name = "eml_conjtimes";
  c5_info[60].dominantType = "double";
  c5_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_conjtimes.m";
  c5_info[60].fileTimeLo = 1286818696U;
  c5_info[60].fileTimeHi = 0U;
  c5_info[60].mFileTimeLo = 0U;
  c5_info[60].mFileTimeHi = 0U;
  c5_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m";
  c5_info[61].name = "eml_zrot_cols";
  c5_info[61].dominantType = "int32";
  c5_info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m";
  c5_info[61].fileTimeLo = 1286818826U;
  c5_info[61].fileTimeHi = 0U;
  c5_info[61].mFileTimeLo = 0U;
  c5_info[61].mFileTimeHi = 0U;
  c5_info[62].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m";
  c5_info[62].name = "eml_matlab_zhgeqz";
  c5_info[62].dominantType = "int32";
  c5_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c5_info[62].fileTimeLo = 1286818820U;
  c5_info[62].fileTimeHi = 0U;
  c5_info[62].mFileTimeLo = 0U;
  c5_info[62].mFileTimeHi = 0U;
  c5_info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c5_info[63].name = "eml_matlab_zlanhs";
  c5_info[63].dominantType = "int32";
  c5_info[63].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m";
  c5_info[63].fileTimeLo = 1286818820U;
  c5_info[63].fileTimeHi = 0U;
  c5_info[63].mFileTimeLo = 0U;
  c5_info[63].mFileTimeHi = 0U;
}

static void c5_b_info_helper(c5_ResolvedFunctionInfo c5_info[106])
{
  c5_info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m";
  c5_info[64].name = "mod";
  c5_info[64].dominantType = "int32";
  c5_info[64].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c5_info[64].fileTimeLo = 1286818744U;
  c5_info[64].fileTimeHi = 0U;
  c5_info[64].mFileTimeLo = 0U;
  c5_info[64].mFileTimeHi = 0U;
  c5_info[65].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m";
  c5_info[65].name = "eml_scalar_mod";
  c5_info[65].dominantType = "int32";
  c5_info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_mod.m";
  c5_info[65].fileTimeLo = 1286818730U;
  c5_info[65].fileTimeHi = 0U;
  c5_info[65].mFileTimeLo = 0U;
  c5_info[65].mFileTimeHi = 0U;
  c5_info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c5_info[66].name = "eml_guarded_inf";
  c5_info[66].dominantType = "char";
  c5_info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m";
  c5_info[66].fileTimeLo = 1286818776U;
  c5_info[66].fileTimeHi = 0U;
  c5_info[66].mFileTimeLo = 0U;
  c5_info[66].mFileTimeHi = 0U;
  c5_info[67].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt";
  c5_info[67].name = "eml_scalar_hypot";
  c5_info[67].dominantType = "double";
  c5_info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m";
  c5_info[67].fileTimeLo = 1286818726U;
  c5_info[67].fileTimeHi = 0U;
  c5_info[67].mFileTimeLo = 0U;
  c5_info[67].mFileTimeHi = 0U;
  c5_info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c5_info[68].name = "isequal";
  c5_info[68].dominantType = "double";
  c5_info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c5_info[68].fileTimeLo = 1286818758U;
  c5_info[68].fileTimeHi = 0U;
  c5_info[68].mFileTimeLo = 0U;
  c5_info[68].mFileTimeHi = 0U;
  c5_info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c5_info[69].name = "eml_isequal_core";
  c5_info[69].dominantType = "double";
  c5_info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c5_info[69].fileTimeLo = 1286818786U;
  c5_info[69].fileTimeHi = 0U;
  c5_info[69].mFileTimeLo = 0U;
  c5_info[69].mFileTimeHi = 0U;
  c5_info[70].context = "";
  c5_info[70].name = "sort";
  c5_info[70].dominantType = "double";
  c5_info[70].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c5_info[70].fileTimeLo = 1303146208U;
  c5_info[70].fileTimeHi = 0U;
  c5_info[70].mFileTimeLo = 0U;
  c5_info[70].mFileTimeHi = 0U;
  c5_info[71].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sort.m";
  c5_info[71].name = "eml_sort";
  c5_info[71].dominantType = "double";
  c5_info[71].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c5_info[71].fileTimeLo = 1305318002U;
  c5_info[71].fileTimeHi = 0U;
  c5_info[71].mFileTimeLo = 0U;
  c5_info[71].mFileTimeHi = 0U;
  c5_info[72].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c5_info[72].name = "eml_nonsingleton_dim";
  c5_info[72].dominantType = "double";
  c5_info[72].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_nonsingleton_dim.m";
  c5_info[72].fileTimeLo = 1286818788U;
  c5_info[72].fileTimeHi = 0U;
  c5_info[72].mFileTimeLo = 0U;
  c5_info[72].mFileTimeHi = 0U;
  c5_info[73].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c5_info[73].name = "eml_assert_valid_dim";
  c5_info[73].dominantType = "double";
  c5_info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_dim.m";
  c5_info[73].fileTimeLo = 1286818694U;
  c5_info[73].fileTimeHi = 0U;
  c5_info[73].mFileTimeLo = 0U;
  c5_info[73].mFileTimeHi = 0U;
  c5_info[74].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort.m";
  c5_info[74].name = "eml_sort_idx";
  c5_info[74].dominantType = "double";
  c5_info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c5_info[74].fileTimeLo = 1305318004U;
  c5_info[74].fileTimeHi = 0U;
  c5_info[74].mFileTimeLo = 0U;
  c5_info[74].mFileTimeHi = 0U;
  c5_info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c5_info[75].name = "eml_size_ispow2";
  c5_info[75].dominantType = "int32";
  c5_info[75].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c5_info[75].fileTimeLo = 1286818798U;
  c5_info[75].fileTimeHi = 0U;
  c5_info[75].mFileTimeLo = 0U;
  c5_info[75].mFileTimeHi = 0U;
  c5_info[76].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_size_ispow2.m";
  c5_info[76].name = "eml_unsigned_class";
  c5_info[76].dominantType = "char";
  c5_info[76].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c5_info[76].fileTimeLo = 1286818800U;
  c5_info[76].fileTimeHi = 0U;
  c5_info[76].mFileTimeLo = 0U;
  c5_info[76].mFileTimeHi = 0U;
  c5_info[77].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_idx.m";
  c5_info[77].name = "eml_sort_le";
  c5_info[77].dominantType = "int32";
  c5_info[77].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m";
  c5_info[77].fileTimeLo = 1292190510U;
  c5_info[77].fileTimeHi = 0U;
  c5_info[77].mFileTimeLo = 0U;
  c5_info[77].mFileTimeHi = 0U;
  c5_info[78].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_sort_le.m!eml_sort_ascending_le";
  c5_info[78].name = "eml_relop";
  c5_info[78].dominantType = "function_handle";
  c5_info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m";
  c5_info[78].fileTimeLo = 1292190510U;
  c5_info[78].fileTimeHi = 0U;
  c5_info[78].mFileTimeLo = 0U;
  c5_info[78].mFileTimeHi = 0U;
  c5_info[79].context = "";
  c5_info[79].name = "any";
  c5_info[79].dominantType = "double";
  c5_info[79].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c5_info[79].fileTimeLo = 1286818834U;
  c5_info[79].fileTimeHi = 0U;
  c5_info[79].mFileTimeLo = 0U;
  c5_info[79].mFileTimeHi = 0U;
  c5_info[80].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m";
  c5_info[80].name = "eml_all_or_any";
  c5_info[80].dominantType = "char";
  c5_info[80].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c5_info[80].fileTimeLo = 1286818694U;
  c5_info[80].fileTimeHi = 0U;
  c5_info[80].mFileTimeLo = 0U;
  c5_info[80].mFileTimeHi = 0U;
  c5_info[81].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_all_or_any.m";
  c5_info[81].name = "eml_const_nonsingleton_dim";
  c5_info[81].dominantType = "double";
  c5_info[81].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  c5_info[81].fileTimeLo = 1286818696U;
  c5_info[81].fileTimeHi = 0U;
  c5_info[81].mFileTimeLo = 0U;
  c5_info[81].mFileTimeHi = 0U;
  c5_info[82].context = "";
  c5_info[82].name = "log";
  c5_info[82].dominantType = "double";
  c5_info[82].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c5_info[82].fileTimeLo = 1286818742U;
  c5_info[82].fileTimeHi = 0U;
  c5_info[82].mFileTimeLo = 0U;
  c5_info[82].mFileTimeHi = 0U;
  c5_info[83].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m";
  c5_info[83].name = "eml_scalar_log";
  c5_info[83].dominantType = "double";
  c5_info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c5_info[83].fileTimeLo = 1286818728U;
  c5_info[83].fileTimeHi = 0U;
  c5_info[83].mFileTimeLo = 0U;
  c5_info[83].mFileTimeHi = 0U;
  c5_info[84].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m";
  c5_info[84].name = "eml_scalar_atan2";
  c5_info[84].dominantType = "double";
  c5_info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m";
  c5_info[84].fileTimeLo = 1286818720U;
  c5_info[84].fileTimeHi = 0U;
  c5_info[84].mFileTimeLo = 0U;
  c5_info[84].mFileTimeHi = 0U;
  c5_info[85].context = "";
  c5_info[85].name = "exp";
  c5_info[85].dominantType = "double";
  c5_info[85].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c5_info[85].fileTimeLo = 1286818740U;
  c5_info[85].fileTimeHi = 0U;
  c5_info[85].mFileTimeLo = 0U;
  c5_info[85].mFileTimeHi = 0U;
  c5_info[86].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m";
  c5_info[86].name = "eml_scalar_exp";
  c5_info[86].dominantType = "double";
  c5_info[86].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m";
  c5_info[86].fileTimeLo = 1301328464U;
  c5_info[86].fileTimeHi = 0U;
  c5_info[86].mFileTimeLo = 0U;
  c5_info[86].mFileTimeHi = 0U;
  c5_info[87].context = "";
  c5_info[87].name = "mldivide";
  c5_info[87].dominantType = "double";
  c5_info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c5_info[87].fileTimeLo = 1310137456U;
  c5_info[87].fileTimeHi = 0U;
  c5_info[87].mFileTimeLo = 1289519690U;
  c5_info[87].mFileTimeHi = 0U;
  c5_info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c5_info[88].name = "eml_lusolve";
  c5_info[88].dominantType = "double";
  c5_info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c5_info[88].fileTimeLo = 1305318000U;
  c5_info[88].fileTimeHi = 0U;
  c5_info[88].mFileTimeLo = 0U;
  c5_info[88].mFileTimeHi = 0U;
  c5_info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c5_info[89].name = "eml_xgetrf";
  c5_info[89].dominantType = "int32";
  c5_info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c5_info[89].fileTimeLo = 1286818806U;
  c5_info[89].fileTimeHi = 0U;
  c5_info[89].mFileTimeLo = 0U;
  c5_info[89].mFileTimeHi = 0U;
  c5_info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c5_info[90].name = "eml_lapack_xgetrf";
  c5_info[90].dominantType = "int32";
  c5_info[90].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c5_info[90].fileTimeLo = 1286818810U;
  c5_info[90].fileTimeHi = 0U;
  c5_info[90].mFileTimeLo = 0U;
  c5_info[90].mFileTimeHi = 0U;
  c5_info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c5_info[91].name = "eml_matlab_zgetrf";
  c5_info[91].dominantType = "int32";
  c5_info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c5_info[91].fileTimeLo = 1302688994U;
  c5_info[91].fileTimeHi = 0U;
  c5_info[91].mFileTimeLo = 0U;
  c5_info[91].mFileTimeHi = 0U;
  c5_info[92].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c5_info[92].name = "intmin";
  c5_info[92].dominantType = "char";
  c5_info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c5_info[92].fileTimeLo = 1286818756U;
  c5_info[92].fileTimeHi = 0U;
  c5_info[92].mFileTimeLo = 0U;
  c5_info[92].mFileTimeHi = 0U;
  c5_info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c5_info[93].name = "eml_ixamax";
  c5_info[93].dominantType = "int32";
  c5_info[93].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c5_info[93].fileTimeLo = 1299076770U;
  c5_info[93].fileTimeHi = 0U;
  c5_info[93].mFileTimeLo = 0U;
  c5_info[93].mFileTimeHi = 0U;
  c5_info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m!below_threshold";
  c5_info[94].name = "length";
  c5_info[94].dominantType = "double";
  c5_info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c5_info[94].fileTimeLo = 1303146206U;
  c5_info[94].fileTimeHi = 0U;
  c5_info[94].mFileTimeLo = 0U;
  c5_info[94].mFileTimeHi = 0U;
  c5_info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c5_info[95].name = "eml_refblas_ixamax";
  c5_info[95].dominantType = "int32";
  c5_info[95].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c5_info[95].fileTimeLo = 1299076770U;
  c5_info[95].fileTimeHi = 0U;
  c5_info[95].mFileTimeLo = 0U;
  c5_info[95].mFileTimeHi = 0U;
  c5_info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c5_info[96].name = "eml_xcabs1";
  c5_info[96].dominantType = "double";
  c5_info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c5_info[96].fileTimeLo = 1286818706U;
  c5_info[96].fileTimeHi = 0U;
  c5_info[96].mFileTimeLo = 0U;
  c5_info[96].mFileTimeHi = 0U;
  c5_info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c5_info[97].name = "eml_xswap";
  c5_info[97].dominantType = "int32";
  c5_info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c5_info[97].fileTimeLo = 1299076778U;
  c5_info[97].fileTimeHi = 0U;
  c5_info[97].mFileTimeLo = 0U;
  c5_info[97].mFileTimeHi = 0U;
  c5_info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c5_info[98].name = "eml_refblas_xswap";
  c5_info[98].dominantType = "int32";
  c5_info[98].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c5_info[98].fileTimeLo = 1299076786U;
  c5_info[98].fileTimeHi = 0U;
  c5_info[98].mFileTimeLo = 0U;
  c5_info[98].mFileTimeHi = 0U;
  c5_info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c5_info[99].name = "eml_xgeru";
  c5_info[99].dominantType = "int32";
  c5_info[99].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c5_info[99].fileTimeLo = 1299076774U;
  c5_info[99].fileTimeHi = 0U;
  c5_info[99].mFileTimeLo = 0U;
  c5_info[99].mFileTimeHi = 0U;
  c5_info[100].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgeru.m";
  c5_info[100].name = "eml_refblas_xgeru";
  c5_info[100].dominantType = "int32";
  c5_info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c5_info[100].fileTimeLo = 1299076776U;
  c5_info[100].fileTimeHi = 0U;
  c5_info[100].mFileTimeLo = 0U;
  c5_info[100].mFileTimeHi = 0U;
  c5_info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgeru.m";
  c5_info[101].name = "eml_refblas_xgerx";
  c5_info[101].dominantType = "int32";
  c5_info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c5_info[101].fileTimeLo = 1299076778U;
  c5_info[101].fileTimeHi = 0U;
  c5_info[101].mFileTimeLo = 0U;
  c5_info[101].mFileTimeHi = 0U;
  c5_info[102].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c5_info[102].name = "eml_xtrsm";
  c5_info[102].dominantType = "int32";
  c5_info[102].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c5_info[102].fileTimeLo = 1299076778U;
  c5_info[102].fileTimeHi = 0U;
  c5_info[102].mFileTimeLo = 0U;
  c5_info[102].mFileTimeHi = 0U;
  c5_info[103].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c5_info[103].name = "eml_refblas_xtrsm";
  c5_info[103].dominantType = "int32";
  c5_info[103].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c5_info[103].fileTimeLo = 1299076786U;
  c5_info[103].fileTimeHi = 0U;
  c5_info[103].mFileTimeLo = 0U;
  c5_info[103].mFileTimeHi = 0U;
  c5_info[104].context = "";
  c5_info[104].name = "log10";
  c5_info[104].dominantType = "double";
  c5_info[104].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c5_info[104].fileTimeLo = 1286818744U;
  c5_info[104].fileTimeHi = 0U;
  c5_info[104].mFileTimeLo = 0U;
  c5_info[104].mFileTimeHi = 0U;
  c5_info[105].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log10.m";
  c5_info[105].name = "eml_scalar_log10";
  c5_info[105].dominantType = "double";
  c5_info[105].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log10.m";
  c5_info[105].fileTimeLo = 1286818728U;
  c5_info[105].fileTimeHi = 0U;
  c5_info[105].mFileTimeLo = 0U;
  c5_info[105].mFileTimeHi = 0U;
}

static real_T c5_mrdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_A, real_T c5_B)
{
  return c5_A / c5_B;
}

static real_T c5_rdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_x, real_T c5_y)
{
  return c5_x / c5_y;
}

static void c5_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   real_T c5_I[63504])
{
  int32_T c5_i;
  int32_T c5_b_i;
  for (c5_i = 0; c5_i < 63504; c5_i++) {
    c5_I[c5_i] = 0.0;
  }

  c5_i = 0;
  for (c5_b_i = 0; c5_b_i < 252; c5_b_i++) {
    c5_I[c5_i] = 1.0;
    c5_i += 253;
  }
}

static void c5_b_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[100])
{
  int32_T c5_i;
  int32_T c5_b_i;
  for (c5_i = 0; c5_i < 100; c5_i++) {
    c5_I[c5_i] = 0.0;
  }

  c5_i = 0;
  for (c5_b_i = 0; c5_b_i < 10; c5_b_i++) {
    c5_I[c5_i] = 1.0;
    c5_i += 11;
  }
}

static void c5_c_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[64])
{
  int32_T c5_i;
  int32_T c5_b_i;
  for (c5_i = 0; c5_i < 64; c5_i++) {
    c5_I[c5_i] = 0.0;
  }

  c5_i = 0;
  for (c5_b_i = 0; c5_b_i < 8; c5_b_i++) {
    c5_I[c5_i] = 1.0;
    c5_i += 9;
  }
}

static void c5_logspace(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  real_T c5_d1, real_T c5_d2, real_T c5_y[300])
{
  real_T c5_b_y[300];
  real_T c5_delta1;
  real_T c5_delta2;
  int32_T c5_b_k;
  if (muDoubleScalarAbs(c5_d2 - 3.1415926535897931) < 8.8817841970012523E-16) {
    c5_eml_warning(chartInstance);
  }

  c5_b_y[299] = c5_d2;
  c5_b_y[0] = c5_d1;
  if (((c5_d1 < 0.0 != c5_d2 < 0.0) && (muDoubleScalarAbs(c5_d1) >
        8.9884656743115785E+307)) || (muDoubleScalarAbs(c5_d2) >
       8.9884656743115785E+307)) {
    c5_delta1 = c5_d1 / 299.0;
    c5_delta2 = c5_d2 / 299.0;
    for (c5_b_k = 0; c5_b_k < 298; c5_b_k++) {
      c5_b_y[c5_b_k + 1] = (c5_d1 + c5_delta2 * (1.0 + (real_T)c5_b_k)) -
        c5_delta1 * (1.0 + (real_T)c5_b_k);
    }
  } else {
    c5_delta1 = (c5_d2 - c5_d1) / 299.0;
    for (c5_b_k = 0; c5_b_k < 298; c5_b_k++) {
      c5_b_y[c5_b_k + 1] = c5_d1 + (1.0 + (real_T)c5_b_k) * c5_delta1;
    }
  }

  for (c5_b_k = 0; c5_b_k < 300; c5_b_k++) {
    c5_y[c5_b_k] = muDoubleScalarPower(10.0, c5_b_y[c5_b_k]);
  }
}

static void c5_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i5;
  static char_T c5_varargin_1[32] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'p', 'i', 'S', 'p', 'e', 'c', 'i', 'a', 'l',
    'S', 'u', 'p', 'p', 'o', 'r', 't', 'e', 'd' };

  char_T c5_u[32];
  const mxArray *c5_y = NULL;
  for (c5_i5 = 0; c5_i5 < 32; c5_i5++) {
    c5_u[c5_i5] = c5_varargin_1[c5_i5];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 32), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static void c5_b_rdivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_x, real_T c5_y[2], real_T c5_z[2])
{
  int32_T c5_i6;
  for (c5_i6 = 0; c5_i6 < 2; c5_i6++) {
    c5_z[c5_i6] = c5_x / c5_y[c5_i6];
  }
}

static void c5_diag(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    real_T c5_v[2], real_T c5_d[4])
{
  int32_T c5_j;
  int32_T c5_b_j;
  for (c5_j = 0; c5_j < 4; c5_j++) {
    c5_d[c5_j] = 0.0;
  }

  c5_j = 0;
  for (c5_b_j = 0; c5_b_j < 2; c5_b_j++) {
    c5_d[c5_j] = c5_v[c5_b_j];
    c5_j += 3;
  }
}

static void c5_eml_error(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i7;
  static char_T c5_varargin_1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 's', 'q', 'r', 't', '_', 'd', 'o', 'm', 'a',
    'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c5_u[30];
  const mxArray *c5_y = NULL;
  for (c5_i7 = 0; c5_i7 < 30; c5_i7++) {
    c5_u[c5_i7] = c5_varargin_1[c5_i7];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 30), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static void c5_damp(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    real_T c5_a[49], real_T c5_h, real_T c5_wn[7], real_T c5_z[7])
{
  creal_T c5_r[7];
  int32_T c5_i;
  real_T c5_x[7];
  int8_T c5_idx[7];
  int32_T c5_b_k;
  boolean_T c5_y;
  int8_T c5_idx0[7];
  int32_T c5_i2;
  int32_T c5_j;
  int32_T c5_pEnd;
  int32_T c5_p;
  int32_T c5_q;
  int32_T c5_qEnd;
  int32_T c5_kEnd;
  creal_T c5_b_r[7];
  boolean_T c5_exitg1;
  boolean_T c5_b0;
  real_T c5_b_a;
  real_T c5_b;
  c5_eig(chartInstance, c5_a, c5_r);
  for (c5_i = 0; c5_i < 7; c5_i++) {
    c5_x[c5_i] = c5_r[c5_i].re;
    c5_idx[c5_i] = (int8_T)(c5_i + 1);
  }

  for (c5_b_k = 0; c5_b_k < 5; c5_b_k += 2) {
    c5_i = c5_b_k + 1;
    if ((c5_x[c5_b_k] <= c5_x[c5_i]) || muDoubleScalarIsNaN(c5_x[c5_i])) {
      c5_y = TRUE;
    } else {
      c5_y = FALSE;
    }

    if (!c5_y) {
      c5_idx[c5_b_k] = (int8_T)(c5_b_k + 2);
      c5_idx[c5_b_k + 1] = (int8_T)(c5_b_k + 1);
    }
  }

  for (c5_i = 0; c5_i < 7; c5_i++) {
    c5_idx0[c5_i] = 1;
  }

  c5_i = 2;
  while (c5_i < 7) {
    c5_i2 = c5_i << 1;
    c5_j = 1;
    for (c5_pEnd = 1 + c5_i; c5_pEnd < 8; c5_pEnd = c5_qEnd + c5_i) {
      c5_p = c5_j;
      c5_q = c5_pEnd;
      c5_qEnd = c5_j + c5_i2;
      if (c5_qEnd > 8) {
        c5_qEnd = 8;
      }

      c5_b_k = 1;
      c5_kEnd = c5_qEnd - c5_j;
      while (c5_b_k <= c5_kEnd) {
        sf_mex_lw_bounds_check(c5_p, 1, 7);
        sf_mex_lw_bounds_check(c5_q, 1, 7);
        if ((c5_x[c5_idx[c5_p - 1] - 1] <= c5_x[c5_idx[c5_q - 1] - 1]) ||
            muDoubleScalarIsNaN(c5_x[c5_idx[c5_q - 1] - 1])) {
          c5_y = TRUE;
        } else {
          c5_y = FALSE;
        }

        if (c5_y) {
          c5_idx0[sf_mex_lw_bounds_check(c5_b_k, 1, 7) - 1] = c5_idx[c5_p - 1];
          c5_p++;
          if (c5_p == c5_pEnd) {
            while (c5_q < c5_qEnd) {
              c5_b_k++;
              c5_idx0[sf_mex_lw_bounds_check(c5_b_k, 1, 7) - 1] = c5_idx[c5_q -
                1];
              c5_q++;
            }
          }
        } else {
          c5_idx0[sf_mex_lw_bounds_check(c5_b_k, 1, 7) - 1] = c5_idx[c5_q - 1];
          c5_q++;
          if (c5_q == c5_qEnd) {
            while (c5_p < c5_pEnd) {
              c5_b_k++;
              c5_idx0[sf_mex_lw_bounds_check(c5_b_k, 1, 7) - 1] = c5_idx[c5_p -
                1];
              c5_p++;
            }
          }
        }

        c5_b_k++;
      }

      for (c5_b_k = 1; c5_b_k <= c5_kEnd; c5_b_k++) {
        c5_idx[sf_mex_lw_bounds_check((c5_j + c5_b_k) - 1, 1, 7) - 1] =
          c5_idx0[c5_b_k - 1];
      }

      c5_j = c5_qEnd;
    }

    c5_i = c5_i2;
  }

  for (c5_i = 0; c5_i < 7; c5_i++) {
    c5_b_r[c5_i] = c5_r[c5_idx[c5_i] - 1];
  }

  for (c5_i = 0; c5_i < 7; c5_i++) {
    c5_r[c5_i] = c5_b_r[c5_i];
  }

  c5_y = FALSE;
  c5_b_k = 0;
  c5_exitg1 = 0U;
  while ((c5_exitg1 == 0U) && (c5_b_k < 7)) {
    if (((c5_r[c5_b_k].re == 0.0) && (c5_r[c5_b_k].im == 0.0)) ||
        (muDoubleScalarIsNaN(c5_r[c5_b_k].re) || muDoubleScalarIsNaN(c5_r[c5_b_k]
          .im))) {
      c5_b0 = TRUE;
    } else {
      c5_b0 = FALSE;
    }

    if (!c5_b0) {
      c5_y = TRUE;
      c5_exitg1 = 1U;
    } else {
      c5_b_k++;
    }
  }

  if (c5_y == 0) {
    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_z[c5_i] = 0.0;
      c5_wn[c5_i] = 0.0;
    }
  } else {
    for (c5_b_k = 0; c5_b_k < 7; c5_b_k++) {
      c5_b_a = c5_r[c5_b_k].re;
      c5_b = c5_r[c5_b_k].im;
      if ((c5_r[c5_b_k].im == 0.0) && muDoubleScalarIsNaN(c5_r[c5_b_k].re)) {
      } else if ((muDoubleScalarAbs(c5_r[c5_b_k].re) > 8.9884656743115785E+307) ||
                 (muDoubleScalarAbs(c5_r[c5_b_k].im) > 8.9884656743115785E+307))
      {
        c5_b_a = muDoubleScalarAbs(c5_r[c5_b_k].re / 2.0);
        c5_b = muDoubleScalarAbs(c5_r[c5_b_k].im / 2.0);
        if (c5_b_a < c5_b) {
          c5_b_a /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_b_a * c5_b_a + 1.0);
        } else if (c5_b_a > c5_b) {
          c5_b /= c5_b_a;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_b_a * 1.4142135623730951;
          }
        }

        c5_b_a = muDoubleScalarLog(c5_b) + 0.69314718055994529;
        c5_b = muDoubleScalarAtan2(c5_r[c5_b_k].im, c5_r[c5_b_k].re);
      } else {
        c5_b_a = muDoubleScalarAbs(c5_r[c5_b_k].re);
        c5_b = muDoubleScalarAbs(c5_r[c5_b_k].im);
        if (c5_b_a < c5_b) {
          c5_b_a /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_b_a * c5_b_a + 1.0);
        } else if (c5_b_a > c5_b) {
          c5_b /= c5_b_a;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_b_a;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_b_a * 1.4142135623730951;
          }
        }

        c5_b_a = muDoubleScalarLog(c5_b);
        c5_b = muDoubleScalarAtan2(c5_r[c5_b_k].im, c5_r[c5_b_k].re);
      }

      c5_r[c5_b_k].re = c5_b_a;
      c5_r[c5_b_k].im = c5_b;
    }

    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_b_a = c5_r[c5_i].im;
      if (c5_r[c5_i].im == 0.0) {
        c5_r[c5_i].re /= c5_h;
        c5_r[c5_i].im = 0.0;
      } else if (c5_r[c5_i].re == 0.0) {
        c5_r[c5_i].re = 0.0;
        c5_r[c5_i].im = c5_b_a / c5_h;
      } else {
        c5_r[c5_i].re /= c5_h;
        c5_r[c5_i].im = c5_b_a / c5_h;
      }
    }

    for (c5_b_k = 0; c5_b_k < 7; c5_b_k++) {
      c5_b_a = muDoubleScalarAbs(c5_r[c5_b_k].re);
      c5_b = muDoubleScalarAbs(c5_r[c5_b_k].im);
      if (c5_b_a < c5_b) {
        c5_b_a /= c5_b;
        c5_b *= muDoubleScalarSqrt(c5_b_a * c5_b_a + 1.0);
      } else if (c5_b_a > c5_b) {
        c5_b /= c5_b_a;
        c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_b_a;
      } else {
        if (!muDoubleScalarIsNaN(c5_b)) {
          c5_b = c5_b_a * 1.4142135623730951;
        }
      }

      c5_wn[c5_b_k] = c5_b;
    }

    for (c5_i = 0; c5_i < 7; c5_i++) {
      c5_z[c5_i] = -c5_r[c5_i].re / c5_wn[c5_i];
    }
  }
}

static void c5_eig(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   real_T c5_A[49], creal_T c5_V[7])
{
  int32_T c5_ii;
  creal_T c5_b_A[49];
  real_T c5_info;
  real_T c5_anrm;
  boolean_T c5_exitg7;
  real_T c5_a;
  real_T c5_b;
  creal_T c5_beta1[7];
  boolean_T c5_ilascl;
  real_T c5_anrmto;
  boolean_T c5_notdone;
  real_T c5_cfrom1;
  real_T c5_cto1;
  real_T c5_mul;
  int32_T c5_ilo;
  int32_T c5_ihi;
  int32_T c5_exitg2;
  int32_T c5_i;
  int32_T c5_j;
  boolean_T c5_exitg5;
  int32_T c5_nzcount;
  int32_T c5_jrow;
  boolean_T c5_exitg6;
  boolean_T c5_c_A;
  boolean_T c5_guard2 = FALSE;
  creal_T c5_atmp;
  int32_T c5_exitg1;
  boolean_T c5_exitg3;
  boolean_T c5_exitg4;
  boolean_T c5_guard1 = FALSE;
  int32_T c5_jrowm1;
  creal_T c5_s;
  static creal_T c5_dc0 = { 0.0, 0.0 };

  creal_T c5_d_A;
  creal_T c5_e_A;
  for (c5_ii = 0; c5_ii < 49; c5_ii++) {
    c5_b_A[c5_ii].re = c5_A[c5_ii];
    c5_b_A[c5_ii].im = 0.0;
  }

  c5_info = 0.0;
  c5_anrm = 0.0;
  c5_ii = 0;
  c5_exitg7 = 0U;
  while ((c5_exitg7 == 0U) && (c5_ii < 49)) {
    c5_a = muDoubleScalarAbs(c5_b_A[c5_ii].re);
    if (c5_a > 0.0) {
      c5_b = 0.0 / c5_a;
      c5_a *= muDoubleScalarSqrt(c5_b * c5_b + 1.0);
    } else {
      c5_a *= 1.4142135623730951;
    }

    if (muDoubleScalarIsNaN(c5_a)) {
      c5_anrm = rtNaN;
      c5_exitg7 = 1U;
    } else {
      if (c5_a > c5_anrm) {
        c5_anrm = c5_a;
      }

      c5_ii++;
    }
  }

  if (!((!muDoubleScalarIsInf(c5_anrm)) && (!muDoubleScalarIsNaN(c5_anrm)))) {
    for (c5_ii = 0; c5_ii < 7; c5_ii++) {
      c5_V[c5_ii].re = rtNaN;
      c5_V[c5_ii].im = 0.0;
      c5_beta1[c5_ii].re = rtNaN;
      c5_beta1[c5_ii].im = 0.0;
    }
  } else {
    c5_ilascl = FALSE;
    c5_anrmto = c5_anrm;
    if ((c5_anrm > 0.0) && (c5_anrm < 6.7178761075670888E-139)) {
      c5_anrmto = 6.7178761075670888E-139;
      c5_ilascl = TRUE;
    } else {
      if (c5_anrm > 1.4885657073574029E+138) {
        c5_anrmto = 1.4885657073574029E+138;
        c5_ilascl = TRUE;
      }
    }

    if (c5_ilascl) {
      c5_a = c5_anrm;
      c5_b = c5_anrmto;
      c5_notdone = TRUE;
      while (c5_notdone) {
        c5_cfrom1 = c5_a * 2.0041683600089728E-292;
        c5_cto1 = c5_b / 4.9896007738368E+291;
        if ((c5_cfrom1 > c5_b) && (c5_b != 0.0)) {
          c5_mul = 2.0041683600089728E-292;
          c5_a = c5_cfrom1;
        } else if (c5_cto1 > c5_a) {
          c5_mul = 4.9896007738368E+291;
          c5_b = c5_cto1;
        } else {
          c5_mul = c5_b / c5_a;
          c5_notdone = FALSE;
        }

        for (c5_ii = 0; c5_ii < 49; c5_ii++) {
          c5_b_A[c5_ii].re *= c5_mul;
          c5_b_A[c5_ii].im *= c5_mul;
        }
      }
    }

    c5_ilo = 1;
    c5_ihi = 7;
    do {
      c5_exitg2 = 0U;
      c5_i = 0;
      c5_j = 0;
      c5_notdone = FALSE;
      c5_ii = c5_ihi;
      c5_exitg5 = 0U;
      while ((c5_exitg5 == 0U) && (c5_ii > 0)) {
        c5_nzcount = 0;
        c5_i = c5_ii;
        c5_j = c5_ihi;
        c5_jrow = 1;
        c5_exitg6 = 0U;
        while ((c5_exitg6 == 0U) && (c5_jrow <= c5_ihi)) {
          c5_c_A = ((c5_b_A[(sf_mex_lw_bounds_check(c5_ii, 1, 7) + 7 * (c5_jrow
            - 1)) - 1].re != 0.0) || (c5_b_A[(sf_mex_lw_bounds_check(c5_ii, 1, 7)
                      + 7 * (c5_jrow - 1)) - 1].im != 0.0));
          c5_guard2 = FALSE;
          if (c5_c_A || (c5_ii == c5_jrow)) {
            if (c5_nzcount == 0) {
              c5_j = c5_jrow;
              c5_nzcount = 1;
              c5_guard2 = TRUE;
            } else {
              c5_nzcount = 2;
              c5_exitg6 = 1U;
            }
          } else {
            c5_guard2 = TRUE;
          }

          if (c5_guard2 == TRUE) {
            c5_jrow++;
          }
        }

        if (c5_nzcount < 2) {
          c5_notdone = TRUE;
          c5_exitg5 = 1U;
        } else {
          c5_ii--;
        }
      }

      if (!c5_notdone) {
        c5_exitg2 = 2U;
      } else {
        if (c5_i != c5_ihi) {
          for (c5_ii = 0; c5_ii < 7; c5_ii++) {
            c5_atmp = c5_b_A[(sf_mex_lw_bounds_check(c5_i, 1, 7) + 7 * c5_ii) -
              1];
            c5_b_A[(c5_i + 7 * c5_ii) - 1] = c5_b_A[(sf_mex_lw_bounds_check
              (c5_ihi, 1, 7) + 7 * c5_ii) - 1];
            c5_b_A[(c5_ihi + 7 * c5_ii) - 1] = c5_atmp;
          }
        }

        if (c5_j != c5_ihi) {
          for (c5_ii = 0; c5_ii + 1 <= c5_ihi; c5_ii++) {
            c5_atmp = c5_b_A[c5_ii + 7 * (sf_mex_lw_bounds_check(c5_j, 1, 7) - 1)];
            c5_b_A[c5_ii + 7 * (c5_j - 1)] = c5_b_A[c5_ii + 7 * (c5_ihi - 1)];
            c5_b_A[c5_ii + 7 * (c5_ihi - 1)] = c5_atmp;
          }
        }

        sf_mex_lw_bounds_check(c5_ihi, 1, 7);
        c5_ihi--;
        if (c5_ihi == 1) {
          c5_exitg2 = 1U;
        }
      }
    } while (c5_exitg2 == 0U);

    if (c5_exitg2 == 1U) {
    } else {
      do {
        c5_exitg1 = 0U;
        c5_i = 0;
        c5_j = 0;
        c5_notdone = FALSE;
        c5_jrow = c5_ilo;
        c5_exitg3 = 0U;
        while ((c5_exitg3 == 0U) && (c5_jrow <= c5_ihi)) {
          c5_nzcount = 0;
          c5_i = c5_ihi;
          c5_j = c5_jrow;
          c5_ii = c5_ilo;
          c5_exitg4 = 0U;
          while ((c5_exitg4 == 0U) && (c5_ii <= c5_ihi)) {
            c5_c_A = ((c5_b_A[(c5_ii + 7 * (sf_mex_lw_bounds_check(c5_jrow, 1, 7)
              - 1)) - 1].re != 0.0) || (c5_b_A[(c5_ii + 7 *
                        (sf_mex_lw_bounds_check(c5_jrow, 1, 7) - 1)) - 1].im !=
                       0.0));
            c5_guard1 = FALSE;
            if (c5_c_A || (c5_ii == c5_jrow)) {
              if (c5_nzcount == 0) {
                c5_i = c5_ii;
                c5_nzcount = 1;
                c5_guard1 = TRUE;
              } else {
                c5_nzcount = 2;
                c5_exitg4 = 1U;
              }
            } else {
              c5_guard1 = TRUE;
            }

            if (c5_guard1 == TRUE) {
              c5_ii++;
            }
          }

          if (c5_nzcount < 2) {
            c5_notdone = TRUE;
            c5_exitg3 = 1U;
          } else {
            c5_jrow++;
          }
        }

        if (!c5_notdone) {
          c5_exitg1 = 1U;
        } else {
          if (c5_i != c5_ilo) {
            for (c5_ii = c5_ilo - 1; c5_ii + 1 < 8; c5_ii++) {
              c5_atmp = c5_b_A[(sf_mex_lw_bounds_check(c5_i, 1, 7) + 7 * c5_ii)
                - 1];
              c5_b_A[(c5_i + 7 * c5_ii) - 1] = c5_b_A[(sf_mex_lw_bounds_check
                (c5_ilo, 1, 7) + 7 * c5_ii) - 1];
              c5_b_A[(c5_ilo + 7 * c5_ii) - 1] = c5_atmp;
            }
          }

          if (c5_j != c5_ilo) {
            for (c5_ii = 0; c5_ii + 1 <= c5_ihi; c5_ii++) {
              c5_atmp = c5_b_A[c5_ii + 7 * (sf_mex_lw_bounds_check(c5_j, 1, 7) -
                1)];
              c5_b_A[c5_ii + 7 * (c5_j - 1)] = c5_b_A[c5_ii + 7 *
                (sf_mex_lw_bounds_check(c5_ilo, 1, 7) - 1)];
              c5_b_A[c5_ii + 7 * (c5_ilo - 1)] = c5_atmp;
            }
          }

          sf_mex_lw_bounds_check(c5_ilo, 1, 7);
          c5_ilo++;
          if (c5_ilo == c5_ihi) {
            c5_exitg1 = 1U;
          }
        }
      } while (c5_exitg1 == 0U);
    }

    if (!(c5_ihi < c5_ilo + 2)) {
      c5_ii = c5_ilo - 1;
      while (c5_ii + 1 < c5_ihi - 1) {
        c5_nzcount = c5_ii + 1;
        c5_jrow = c5_ihi - 1;
        while (c5_jrow + 1 > c5_nzcount + 1) {
          c5_jrowm1 = c5_jrow - 1;
          c5_eml_matlab_zlartg(chartInstance, c5_b_A[c5_jrowm1 + 7 * c5_ii],
                               c5_b_A[c5_jrow + 7 * c5_ii], &c5_a, &c5_s,
                               &c5_atmp);
          c5_b_A[c5_jrowm1 + 7 * c5_ii] = c5_atmp;
          c5_b_A[c5_jrow + 7 * c5_ii] = c5_dc0;
          for (c5_j = c5_nzcount; c5_j + 1 <= c5_ihi; c5_j++) {
            c5_atmp.re = c5_a * c5_b_A[c5_jrowm1 + 7 * c5_j].re;
            c5_atmp.im = c5_a * c5_b_A[c5_jrowm1 + 7 * c5_j].im;
            c5_b = c5_s.re * c5_b_A[c5_jrow + 7 * c5_j].re - c5_s.im *
              c5_b_A[c5_jrow + 7 * c5_j].im;
            c5_cfrom1 = c5_s.re * c5_b_A[c5_jrow + 7 * c5_j].im + c5_s.im *
              c5_b_A[c5_jrow + 7 * c5_j].re;
            c5_d_A = c5_b_A[c5_jrowm1 + 7 * c5_j];
            c5_e_A = c5_b_A[c5_jrowm1 + 7 * c5_j];
            c5_b_A[c5_jrow + 7 * c5_j].re = c5_a * c5_b_A[c5_jrow + 7 * c5_j].re
              - (c5_s.re * c5_b_A[c5_jrowm1 + 7 * c5_j].re + c5_s.im *
                 c5_b_A[c5_jrowm1 + 7 * c5_j].im);
            c5_b_A[c5_jrow + 7 * c5_j].im = c5_a * c5_b_A[c5_jrow + 7 * c5_j].im
              - (c5_s.re * c5_d_A.im - c5_s.im * c5_e_A.re);
            c5_b_A[c5_jrowm1 + 7 * c5_j].re = c5_atmp.re + c5_b;
            c5_b_A[c5_jrowm1 + 7 * c5_j].im = c5_atmp.im + c5_cfrom1;
          }

          c5_s.re = -c5_s.re;
          c5_s.im = -c5_s.im;
          for (c5_i = c5_ilo - 1; c5_i + 1 <= c5_ihi; c5_i++) {
            c5_atmp.re = c5_a * c5_b_A[c5_i + 7 * c5_jrow].re;
            c5_atmp.im = c5_a * c5_b_A[c5_i + 7 * c5_jrow].im;
            c5_b = c5_s.re * c5_b_A[c5_i + 7 * c5_jrowm1].re - c5_s.im *
              c5_b_A[c5_i + 7 * c5_jrowm1].im;
            c5_cfrom1 = c5_s.re * c5_b_A[c5_i + 7 * c5_jrowm1].im + c5_s.im *
              c5_b_A[c5_i + 7 * c5_jrowm1].re;
            c5_d_A = c5_b_A[c5_i + 7 * c5_jrow];
            c5_e_A = c5_b_A[c5_i + 7 * c5_jrow];
            c5_b_A[c5_i + 7 * c5_jrowm1].re = c5_a * c5_b_A[c5_i + 7 * c5_jrowm1]
              .re - (c5_s.re * c5_b_A[c5_i + 7 * c5_jrow].re + c5_s.im *
                     c5_b_A[c5_i + 7 * c5_jrow].im);
            c5_b_A[c5_i + 7 * c5_jrowm1].im = c5_a * c5_b_A[c5_i + 7 * c5_jrowm1]
              .im - (c5_s.re * c5_d_A.im - c5_s.im * c5_e_A.re);
            c5_b_A[c5_i + 7 * c5_jrow].re = c5_atmp.re + c5_b;
            c5_b_A[c5_i + 7 * c5_jrow].im = c5_atmp.im + c5_cfrom1;
          }

          c5_jrow = c5_jrowm1;
        }

        c5_ii = c5_nzcount;
      }
    }

    c5_eml_matlab_zhgeqz(chartInstance, c5_b_A, c5_ilo, c5_ihi, &c5_info, c5_V,
                         c5_beta1);
    if ((!(c5_info != 0.0)) && c5_ilascl) {
      c5_notdone = TRUE;
      while (c5_notdone) {
        c5_cfrom1 = c5_anrmto * 2.0041683600089728E-292;
        c5_cto1 = c5_anrm / 4.9896007738368E+291;
        if ((c5_cfrom1 > c5_anrm) && (c5_anrm != 0.0)) {
          c5_mul = 2.0041683600089728E-292;
          c5_anrmto = c5_cfrom1;
        } else if (c5_cto1 > c5_anrmto) {
          c5_mul = 4.9896007738368E+291;
          c5_anrm = c5_cto1;
        } else {
          c5_mul = c5_anrm / c5_anrmto;
          c5_notdone = FALSE;
        }

        for (c5_ii = 0; c5_ii < 7; c5_ii++) {
          c5_V[c5_ii].re *= c5_mul;
          c5_V[c5_ii].im *= c5_mul;
        }
      }
    }
  }

  for (c5_ii = 0; c5_ii < 7; c5_ii++) {
    c5_cto1 = c5_V[c5_ii].re;
    c5_mul = c5_V[c5_ii].im;
    if (c5_beta1[c5_ii].im == 0.0) {
      if (c5_V[c5_ii].im == 0.0) {
        c5_V[c5_ii].re /= c5_beta1[c5_ii].re;
        c5_V[c5_ii].im = 0.0;
      } else if (c5_V[c5_ii].re == 0.0) {
        c5_V[c5_ii].re = 0.0;
        c5_V[c5_ii].im = c5_mul / c5_beta1[c5_ii].re;
      } else {
        c5_V[c5_ii].re /= c5_beta1[c5_ii].re;
        c5_V[c5_ii].im = c5_mul / c5_beta1[c5_ii].re;
      }
    } else if (c5_beta1[c5_ii].re == 0.0) {
      if (c5_V[c5_ii].re == 0.0) {
        c5_V[c5_ii].re = c5_V[c5_ii].im / c5_beta1[c5_ii].im;
        c5_V[c5_ii].im = 0.0;
      } else if (c5_V[c5_ii].im == 0.0) {
        c5_V[c5_ii].re = 0.0;
        c5_V[c5_ii].im = -(c5_cto1 / c5_beta1[c5_ii].im);
      } else {
        c5_V[c5_ii].re = c5_V[c5_ii].im / c5_beta1[c5_ii].im;
        c5_V[c5_ii].im = -(c5_cto1 / c5_beta1[c5_ii].im);
      }
    } else {
      c5_cfrom1 = muDoubleScalarAbs(c5_beta1[c5_ii].re);
      c5_a = muDoubleScalarAbs(c5_beta1[c5_ii].im);
      if (c5_cfrom1 > c5_a) {
        c5_a = c5_beta1[c5_ii].im / c5_beta1[c5_ii].re;
        c5_b = c5_beta1[c5_ii].re + c5_a * c5_beta1[c5_ii].im;
        c5_V[c5_ii].re = (c5_V[c5_ii].re + c5_a * c5_V[c5_ii].im) / c5_b;
        c5_V[c5_ii].im = (c5_mul - c5_a * c5_cto1) / c5_b;
      } else if (c5_a == c5_cfrom1) {
        c5_a = c5_beta1[c5_ii].re > 0.0 ? 0.5 : -0.5;
        c5_b = c5_beta1[c5_ii].im > 0.0 ? 0.5 : -0.5;
        c5_V[c5_ii].re = (c5_V[c5_ii].re * c5_a + c5_V[c5_ii].im * c5_b) /
          c5_cfrom1;
        c5_V[c5_ii].im = (c5_mul * c5_a - c5_cto1 * c5_b) / c5_cfrom1;
      } else {
        c5_a = c5_beta1[c5_ii].re / c5_beta1[c5_ii].im;
        c5_b = c5_beta1[c5_ii].im + c5_a * c5_beta1[c5_ii].re;
        c5_V[c5_ii].re = (c5_a * c5_V[c5_ii].re + c5_V[c5_ii].im) / c5_b;
        c5_V[c5_ii].im = (c5_a * c5_mul - c5_cto1) / c5_b;
      }
    }
  }

  if (c5_info < 0.0) {
    c5_b_eml_warning(chartInstance);
  } else {
    if (c5_info > 0.0) {
      c5_c_eml_warning(chartInstance);
    }
  }
}

static void c5_eml_matlab_zlartg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn,
  creal_T *c5_r)
{
  real_T c5_scale;
  real_T c5_dr;
  real_T c5_f2;
  real_T c5_fs_re;
  real_T c5_fs_im;
  real_T c5_gs_re;
  real_T c5_gs_im;
  int32_T c5_count;
  int32_T c5_rescaledir;
  boolean_T c5_guard1 = FALSE;
  real_T c5_b;
  real_T c5_g2s;
  c5_scale = muDoubleScalarAbs(c5_f.re);
  c5_dr = muDoubleScalarAbs(c5_f.im);
  if (c5_dr > c5_scale) {
    c5_scale = c5_dr;
  }

  c5_dr = muDoubleScalarAbs(c5_g.re);
  c5_f2 = muDoubleScalarAbs(c5_g.im);
  if (c5_f2 > c5_dr) {
    c5_dr = c5_f2;
  }

  if (c5_dr > c5_scale) {
    c5_scale = c5_dr;
  }

  c5_fs_re = c5_f.re;
  c5_fs_im = c5_f.im;
  c5_gs_re = c5_g.re;
  c5_gs_im = c5_g.im;
  c5_count = 0;
  c5_rescaledir = 0;
  c5_guard1 = FALSE;
  if (c5_scale >= 7.4428285367870146E+137) {
    do {
      c5_count++;
      c5_fs_re *= 1.3435752215134178E-138;
      c5_fs_im *= 1.3435752215134178E-138;
      c5_gs_re *= 1.3435752215134178E-138;
      c5_gs_im *= 1.3435752215134178E-138;
      c5_scale *= 1.3435752215134178E-138;
    } while (!(c5_scale < 7.4428285367870146E+137));

    c5_rescaledir = 1;
    c5_guard1 = TRUE;
  } else if (c5_scale <= 1.3435752215134178E-138) {
    if ((c5_g.re == 0.0) && (c5_g.im == 0.0)) {
      *c5_cs = 1.0;
      c5_sn->re = 0.0;
      c5_sn->im = 0.0;
      *c5_r = c5_f;
    } else {
      do {
        c5_count++;
        c5_fs_re *= 7.4428285367870146E+137;
        c5_fs_im *= 7.4428285367870146E+137;
        c5_gs_re *= 7.4428285367870146E+137;
        c5_gs_im *= 7.4428285367870146E+137;
        c5_scale *= 7.4428285367870146E+137;
      } while (!(c5_scale > 1.3435752215134178E-138));

      c5_rescaledir = -1;
      c5_guard1 = TRUE;
    }
  } else {
    c5_guard1 = TRUE;
  }

  if (c5_guard1 == TRUE) {
    c5_f2 = c5_fs_re * c5_fs_re + c5_fs_im * c5_fs_im;
    c5_scale = c5_gs_re * c5_gs_re + c5_gs_im * c5_gs_im;
    c5_dr = c5_scale;
    if (1.0 > c5_scale) {
      c5_dr = 1.0;
    }

    if (c5_f2 <= c5_dr * 2.0041683600089728E-292) {
      if ((c5_f.re == 0.0) && (c5_f.im == 0.0)) {
        *c5_cs = 0.0;
        c5_f2 = muDoubleScalarAbs(c5_g.re);
        c5_b = muDoubleScalarAbs(c5_g.im);
        if (c5_f2 < c5_b) {
          c5_f2 /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
        } else if (c5_f2 > c5_b) {
          c5_b /= c5_f2;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_f2 * 1.4142135623730951;
          }
        }

        c5_r->re = c5_b;
        c5_r->im = 0.0;
        c5_f2 = muDoubleScalarAbs(c5_gs_re);
        c5_b = muDoubleScalarAbs(c5_gs_im);
        if (c5_f2 < c5_b) {
          c5_f2 /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
        } else if (c5_f2 > c5_b) {
          c5_b /= c5_f2;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_f2 * 1.4142135623730951;
          }
        }

        c5_sn->re = c5_gs_re / c5_b;
        c5_sn->im = -c5_gs_im / c5_b;
      } else {
        c5_f2 = muDoubleScalarAbs(c5_fs_re);
        c5_b = muDoubleScalarAbs(c5_fs_im);
        if (c5_f2 < c5_b) {
          c5_f2 /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
        } else if (c5_f2 > c5_b) {
          c5_b /= c5_f2;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_f2 * 1.4142135623730951;
          }
        }

        if (c5_scale < 0.0) {
          c5_eml_error(chartInstance);
        }

        c5_g2s = muDoubleScalarSqrt(c5_scale);
        *c5_cs = c5_b / c5_g2s;
        c5_dr = muDoubleScalarAbs(c5_f.re);
        c5_f2 = muDoubleScalarAbs(c5_f.im);
        if (c5_f2 > c5_dr) {
          c5_dr = c5_f2;
        }

        if (c5_dr > 1.0) {
          c5_f2 = muDoubleScalarAbs(c5_f.re);
          c5_b = muDoubleScalarAbs(c5_f.im);
          if (c5_f2 < c5_b) {
            c5_f2 /= c5_b;
            c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
          } else if (c5_f2 > c5_b) {
            c5_b /= c5_f2;
            c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
          } else {
            if (!muDoubleScalarIsNaN(c5_b)) {
              c5_b = c5_f2 * 1.4142135623730951;
            }
          }

          c5_fs_re = c5_f.re / c5_b;
          c5_fs_im = c5_f.im / c5_b;
        } else {
          c5_dr = 7.4428285367870146E+137 * c5_f.re;
          c5_scale = 7.4428285367870146E+137 * c5_f.im;
          c5_f2 = muDoubleScalarAbs(c5_dr);
          c5_b = muDoubleScalarAbs(c5_scale);
          if (c5_f2 < c5_b) {
            c5_f2 /= c5_b;
            c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
          } else if (c5_f2 > c5_b) {
            c5_b /= c5_f2;
            c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
          } else {
            if (!muDoubleScalarIsNaN(c5_b)) {
              c5_b = c5_f2 * 1.4142135623730951;
            }
          }

          c5_fs_re = c5_dr / c5_b;
          c5_fs_im = c5_scale / c5_b;
        }

        c5_gs_re /= c5_g2s;
        c5_gs_im = -c5_gs_im / c5_g2s;
        c5_sn->re = c5_fs_re * c5_gs_re - c5_fs_im * c5_gs_im;
        c5_sn->im = c5_fs_re * c5_gs_im + c5_fs_im * c5_gs_re;
        c5_r->re = *c5_cs * c5_f.re + (c5_sn->re * c5_g.re - c5_sn->im * c5_g.im);
        c5_r->im = *c5_cs * c5_f.im + (c5_sn->re * c5_g.im + c5_sn->im * c5_g.re);
      }
    } else {
      c5_dr = 1.0 + c5_scale / c5_f2;
      if (c5_dr < 0.0) {
        c5_eml_error(chartInstance);
      }

      c5_b = muDoubleScalarSqrt(c5_dr);
      c5_r->re = c5_b * c5_fs_re;
      c5_r->im = c5_b * c5_fs_im;
      *c5_cs = 1.0 / c5_b;
      c5_b = c5_f2 + c5_scale;
      c5_f2 = c5_r->re / c5_b;
      c5_dr = c5_r->im / c5_b;
      c5_sn->re = c5_f2 * c5_gs_re - c5_dr * -c5_gs_im;
      c5_sn->im = c5_f2 * -c5_gs_im + c5_dr * c5_gs_re;
      if (c5_rescaledir > 0) {
        for (c5_rescaledir = 1; c5_rescaledir <= c5_count; c5_rescaledir++) {
          c5_r->re *= 7.4428285367870146E+137;
          c5_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c5_rescaledir < 0) {
          for (c5_rescaledir = 1; c5_rescaledir <= c5_count; c5_rescaledir++) {
            c5_r->re *= 1.3435752215134178E-138;
            c5_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c5_eml_matlab_zhgeqz(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_A[49], int32_T c5_ilo, int32_T c5_ihi, real_T
  *c5_info, creal_T c5_alpha1[7], creal_T c5_beta1[7])
{
  static creal_T c5_dc1 = { 0.0, 0.0 };

  int32_T c5_jm1;
  creal_T c5_b_A[49];
  real_T c5_eshift_re;
  real_T c5_eshift_im;
  static creal_T c5_dc2 = { 0.0, 0.0 };

  creal_T c5_ctemp;
  real_T c5_rho_re;
  real_T c5_rho_im;
  real_T c5_anorm;
  real_T c5_a;
  real_T c5_atol;
  real_T c5_ascale;
  boolean_T c5_failed;
  int32_T c5_j;
  boolean_T c5_guard1 = FALSE;
  boolean_T c5_guard2 = FALSE;
  int32_T c5_ifirst;
  int32_T c5_istart;
  int32_T c5_ilast;
  int32_T c5_ilastm1;
  int32_T c5_ifrstm;
  int32_T c5_iiter;
  int32_T c5_maxit;
  boolean_T c5_goto60;
  boolean_T c5_goto70;
  boolean_T c5_goto90;
  int32_T c5_jiter;
  int32_T c5_exitg1;
  boolean_T c5_exitg3;
  boolean_T c5_ilazro;
  boolean_T c5_b_guard1 = FALSE;
  creal_T c5_c_A;
  creal_T c5_t1;
  creal_T c5_d;
  creal_T c5_sigma1;
  real_T c5_sigma2_re;
  real_T c5_sigma2_im;
  real_T c5_b;
  int32_T c5_jp1;
  boolean_T c5_exitg2;
  creal_T c5_d_A;
  int32_T c5_i;
  c5_dc1.re = rtNaN;
  for (c5_jm1 = 0; c5_jm1 < 49; c5_jm1++) {
    c5_b_A[c5_jm1] = c5_A[c5_jm1];
  }

  for (c5_jm1 = 0; c5_jm1 < 7; c5_jm1++) {
    c5_alpha1[c5_jm1].re = 0.0;
    c5_alpha1[c5_jm1].im = 0.0;
    c5_beta1[c5_jm1].re = 1.0;
    c5_beta1[c5_jm1].im = 0.0;
  }

  c5_eshift_re = 0.0;
  c5_eshift_im = 0.0;
  c5_ctemp = c5_dc2;
  c5_rho_re = 0.0;
  c5_rho_im = 0.0;
  c5_anorm = c5_eml_matlab_zlanhs(chartInstance, c5_A, c5_ilo, c5_ihi);
  c5_a = 2.2204460492503131E-16 * c5_anorm;
  c5_atol = 2.2250738585072014E-308;
  if (c5_a > 2.2250738585072014E-308) {
    c5_atol = c5_a;
  }

  c5_a = 2.2250738585072014E-308;
  if (c5_anorm > 2.2250738585072014E-308) {
    c5_a = c5_anorm;
  }

  c5_ascale = 1.0 / c5_a;
  c5_failed = TRUE;
  for (c5_j = c5_ihi + 1; c5_j < 8; c5_j++) {
    c5_alpha1[sf_mex_lw_bounds_check(c5_j, 1, 7) - 1] = c5_A
      [(sf_mex_lw_bounds_check(c5_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c5_j, 1,
          7) - 1)) - 1];
  }

  c5_guard1 = FALSE;
  c5_guard2 = FALSE;
  if (c5_ihi >= c5_ilo) {
    c5_ifirst = c5_ilo;
    c5_istart = c5_ilo;
    c5_ilast = c5_ihi;
    c5_ilastm1 = c5_ihi - 2;
    c5_ifrstm = c5_ilo;
    c5_iiter = 0;
    c5_maxit = 30 * ((c5_ihi - c5_ilo) + 1);
    c5_goto60 = FALSE;
    c5_goto70 = FALSE;
    c5_goto90 = FALSE;
    c5_jiter = 1;
    do {
      c5_exitg1 = 0U;
      if (c5_jiter <= c5_maxit) {
        if (c5_ilast == c5_ilo) {
          c5_goto60 = TRUE;
        } else {
          sf_mex_lw_bounds_check(c5_ilastm1 + 1, 1, 7);
          sf_mex_lw_bounds_check(c5_ilast, 1, 7);
          if (muDoubleScalarAbs(c5_b_A[(c5_ilast + 7 * c5_ilastm1) - 1].re) +
              muDoubleScalarAbs(c5_b_A[(c5_ilast + 7 * c5_ilastm1) - 1].im) <=
              c5_atol) {
            c5_b_A[(c5_ilast + 7 * c5_ilastm1) - 1] = c5_dc2;
            c5_goto60 = TRUE;
          } else {
            c5_j = c5_ilastm1;
            c5_exitg3 = 0U;
            while ((c5_exitg3 == 0U) && (c5_j + 1 >= c5_ilo)) {
              c5_jm1 = c5_j - 1;
              if (c5_j + 1 == c5_ilo) {
                c5_ilazro = TRUE;
              } else {
                sf_mex_lw_bounds_check(c5_jm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c5_j + 1, 1, 7);
                if (muDoubleScalarAbs(c5_b_A[c5_j + 7 * c5_jm1].re) +
                    muDoubleScalarAbs(c5_b_A[c5_j + 7 * c5_jm1].im) <= c5_atol)
                {
                  c5_b_A[c5_j + 7 * c5_jm1] = c5_dc2;
                  c5_ilazro = TRUE;
                } else {
                  c5_ilazro = FALSE;
                }
              }

              if (c5_ilazro) {
                c5_ifirst = c5_j + 1;
                c5_goto70 = TRUE;
                c5_exitg3 = 1U;
              } else {
                c5_j = c5_jm1;
              }
            }
          }
        }

        if (c5_goto60 || c5_goto70) {
          c5_ilazro = TRUE;
        } else {
          c5_ilazro = FALSE;
        }

        if (!c5_ilazro) {
          for (c5_jm1 = 0; c5_jm1 < 7; c5_jm1++) {
            c5_alpha1[c5_jm1] = c5_dc1;
            c5_beta1[c5_jm1] = c5_dc1;
          }

          *c5_info = -1.0;
          c5_exitg1 = 1U;
        } else {
          c5_b_guard1 = FALSE;
          if (c5_goto60) {
            c5_goto60 = FALSE;
            c5_alpha1[sf_mex_lw_bounds_check(c5_ilast, 1, 7) - 1] = c5_b_A
              [(sf_mex_lw_bounds_check(c5_ilast, 1, 7) + 7 *
                (sf_mex_lw_bounds_check(c5_ilast, 1, 7) - 1)) - 1];
            c5_ilast = c5_ilastm1 + 1;
            c5_ilastm1--;
            if (c5_ilast < c5_ilo) {
              c5_failed = FALSE;
              c5_guard2 = TRUE;
              c5_exitg1 = 1U;
            } else {
              c5_iiter = 0;
              c5_eshift_re = 0.0;
              c5_eshift_im = 0.0;
              if (c5_ifrstm > c5_ilast) {
                c5_ifrstm = c5_ilo;
              }

              c5_b_guard1 = TRUE;
            }
          } else {
            if (c5_goto70) {
              c5_goto70 = FALSE;
              c5_iiter++;
              c5_ifrstm = c5_ifirst;
              if (c5_iiter - c5_div_s32_floor(chartInstance, c5_iiter, 10) * 10
                  != 0) {
                sf_mex_lw_bounds_check(c5_ilastm1 + 1, 1, 7);
                sf_mex_lw_bounds_check(c5_ilast, 1, 7);
                c5_c_A.re = -(c5_b_A[(c5_ilast + 7 * (c5_ilast - 1)) - 1].re -
                              c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].re);
                c5_c_A.im = -(c5_b_A[(c5_ilast + 7 * (c5_ilast - 1)) - 1].im -
                              c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].im);
                c5_t1 = c5_eml_div(chartInstance, c5_c_A, 2.0);
                c5_d.re = (c5_t1.re * c5_t1.re - c5_t1.im * c5_t1.im) +
                  (c5_b_A[c5_ilastm1 + 7 * (c5_ilast - 1)].re * c5_b_A[(c5_ilast
                    + 7 * c5_ilastm1) - 1].re - c5_b_A[c5_ilastm1 + 7 *
                   (c5_ilast - 1)].im * c5_b_A[(c5_ilast + 7 * c5_ilastm1) - 1].
                   im);
                c5_d.im = (c5_t1.re * c5_t1.im + c5_t1.im * c5_t1.re) +
                  (c5_b_A[c5_ilastm1 + 7 * (c5_ilast - 1)].re * c5_b_A[(c5_ilast
                    + 7 * c5_ilastm1) - 1].im + c5_b_A[c5_ilastm1 + 7 *
                   (c5_ilast - 1)].im * c5_b_A[(c5_ilast + 7 * c5_ilastm1) - 1].
                   re);
                c5_sqrt(chartInstance, &c5_d);
                c5_sigma1.re = c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].re -
                  (c5_t1.re - c5_d.re);
                c5_sigma1.im = c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].im -
                  (c5_t1.im - c5_d.im);
                c5_sigma2_re = c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].re -
                  (c5_t1.re + c5_d.re);
                c5_sigma2_im = c5_b_A[c5_ilastm1 + 7 * c5_ilastm1].im -
                  (c5_t1.im + c5_d.im);
                c5_a = muDoubleScalarAbs(c5_sigma1.re - c5_b_A[(c5_ilast + 7 *
                  (c5_ilast - 1)) - 1].re);
                c5_anorm = muDoubleScalarAbs(c5_sigma1.im - c5_b_A[(c5_ilast + 7
                  * (c5_ilast - 1)) - 1].im);
                if (c5_a < c5_anorm) {
                  c5_a /= c5_anorm;
                  c5_anorm *= muDoubleScalarSqrt(c5_a * c5_a + 1.0);
                } else if (c5_a > c5_anorm) {
                  c5_anorm /= c5_a;
                  c5_anorm = muDoubleScalarSqrt(c5_anorm * c5_anorm + 1.0) *
                    c5_a;
                } else {
                  if (!muDoubleScalarIsNaN(c5_anorm)) {
                    c5_anorm = c5_a * 1.4142135623730951;
                  }
                }

                c5_a = muDoubleScalarAbs(c5_sigma2_re - c5_b_A[(c5_ilast + 7 *
                  (c5_ilast - 1)) - 1].re);
                c5_b = muDoubleScalarAbs(c5_sigma2_im - c5_b_A[(c5_ilast + 7 *
                  (c5_ilast - 1)) - 1].im);
                if (c5_a < c5_b) {
                  c5_a /= c5_b;
                  c5_b *= muDoubleScalarSqrt(c5_a * c5_a + 1.0);
                } else if (c5_a > c5_b) {
                  c5_b /= c5_a;
                  c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_a;
                } else {
                  if (!muDoubleScalarIsNaN(c5_b)) {
                    c5_b = c5_a * 1.4142135623730951;
                  }
                }

                if (c5_anorm <= c5_b) {
                  c5_sigma2_re = c5_sigma1.re;
                  c5_sigma2_im = c5_sigma1.im;
                  c5_rho_re = c5_t1.re - c5_d.re;
                  c5_rho_im = c5_t1.im - c5_d.im;
                } else {
                  c5_rho_re = c5_t1.re + c5_d.re;
                  c5_rho_im = c5_t1.im + c5_d.im;
                }
              } else {
                c5_eshift_re += c5_b_A[(sf_mex_lw_bounds_check(c5_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c5_ilastm1 + 1, 1, 7) - 1)) - 1].
                  re;
                c5_eshift_im += c5_b_A[(sf_mex_lw_bounds_check(c5_ilast, 1, 7) +
                  7 * (sf_mex_lw_bounds_check(c5_ilastm1 + 1, 1, 7) - 1)) - 1].
                  im;
                c5_sigma2_re = c5_eshift_re;
                c5_sigma2_im = c5_eshift_im;
              }

              c5_j = c5_ilastm1;
              c5_jp1 = c5_ilastm1 + 1;
              c5_exitg2 = 0U;
              while ((c5_exitg2 == 0U) && (c5_j + 1 > c5_ifirst)) {
                c5_jm1 = c5_j - 1;
                c5_istart = c5_j + 1;
                c5_ctemp.re = c5_b_A[c5_j + 7 * c5_j].re - c5_sigma2_re;
                c5_ctemp.im = c5_b_A[c5_j + 7 * c5_j].im - c5_sigma2_im;
                c5_anorm = c5_ascale * (muDoubleScalarAbs(c5_ctemp.re) +
                  muDoubleScalarAbs(c5_ctemp.im));
                sf_mex_lw_bounds_check(c5_jp1 + 1, 1, 7);
                c5_a = c5_ascale * (muDoubleScalarAbs(c5_b_A[c5_jp1 + 7 * c5_j].
                  re) + muDoubleScalarAbs(c5_b_A[c5_jp1 + 7 * c5_j].im));
                c5_b = c5_anorm;
                if (c5_a > c5_anorm) {
                  c5_b = c5_a;
                }

                if ((c5_b < 1.0) && (c5_b != 0.0)) {
                  c5_anorm /= c5_b;
                  c5_a /= c5_b;
                }

                sf_mex_lw_bounds_check(c5_jm1 + 1, 1, 7);
                if ((muDoubleScalarAbs(c5_b_A[c5_j + 7 * c5_jm1].re) +
                     muDoubleScalarAbs(c5_b_A[c5_j + 7 * c5_jm1].im)) * c5_a <=
                    c5_anorm * c5_atol) {
                  c5_goto90 = TRUE;
                  c5_exitg2 = 1U;
                } else {
                  c5_jp1 = c5_j;
                  c5_j = c5_jm1;
                }
              }

              if (!c5_goto90) {
                c5_istart = c5_ifirst;
                if (c5_ifirst == c5_ilastm1 + 1) {
                  c5_ctemp.re = c5_rho_re;
                  c5_ctemp.im = c5_rho_im;
                } else {
                  c5_ctemp.re = c5_b_A[(sf_mex_lw_bounds_check(c5_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c5_ifirst,
                    1, 7) - 1)) - 1].re - c5_sigma2_re;
                  c5_ctemp.im = c5_b_A[(sf_mex_lw_bounds_check(c5_ifirst, 1, 7)
                                        + 7 * (sf_mex_lw_bounds_check(c5_ifirst,
                    1, 7) - 1)) - 1].im - c5_sigma2_im;
                }

                c5_goto90 = TRUE;
              }
            }

            if (c5_goto90) {
              c5_goto90 = FALSE;
              sf_mex_lw_bounds_check(c5_istart, 1, 7);
              sf_mex_lw_bounds_check(c5_istart + 1, 1, 7);
              c5_b_eml_matlab_zlartg(chartInstance, c5_ctemp, c5_b_A[c5_istart +
                7 * (c5_istart - 1)], &c5_a, &c5_sigma1);
              c5_j = c5_istart - 1;
              c5_jm1 = c5_istart - 2;
              while (c5_j + 1 < c5_ilast) {
                c5_jp1 = c5_j + 1;
                if (c5_j + 1 > c5_istart) {
                  c5_c_A = c5_b_A[(sf_mex_lw_bounds_check(c5_j + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c5_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c5_d_A = c5_b_A[(sf_mex_lw_bounds_check(c5_jp1 + 1, 1, 7) + 7 *
                                   (sf_mex_lw_bounds_check(c5_jm1 + 1, 1, 7) - 1))
                    - 1];
                  c5_eml_matlab_zlartg(chartInstance, c5_c_A, c5_d_A, &c5_a,
                                       &c5_sigma1, &c5_t1);
                  c5_b_A[c5_j + 7 * c5_jm1] = c5_t1;
                  c5_b_A[c5_jp1 + 7 * c5_jm1] = c5_dc2;
                }

                for (c5_jm1 = c5_j; c5_jm1 + 1 <= c5_ilast; c5_jm1++) {
                  sf_mex_lw_bounds_check(c5_jm1 + 1, 1, 7);
                  c5_t1.re = c5_a * c5_b_A[c5_j + 7 * c5_jm1].re;
                  c5_t1.im = c5_a * c5_b_A[c5_j + 7 * c5_jm1].im;
                  sf_mex_lw_bounds_check(c5_jp1 + 1, 1, 7);
                  c5_d.re = c5_sigma1.re * c5_b_A[c5_jp1 + 7 * c5_jm1].re -
                    c5_sigma1.im * c5_b_A[c5_jp1 + 7 * c5_jm1].im;
                  c5_d.im = c5_sigma1.re * c5_b_A[c5_jp1 + 7 * c5_jm1].im +
                    c5_sigma1.im * c5_b_A[c5_jp1 + 7 * c5_jm1].re;
                  c5_c_A = c5_b_A[c5_j + 7 * c5_jm1];
                  c5_d_A = c5_b_A[c5_j + 7 * c5_jm1];
                  c5_b_A[c5_jp1 + 7 * c5_jm1].re = c5_a * c5_b_A[c5_jp1 + 7 *
                    c5_jm1].re - (c5_sigma1.re * c5_b_A[c5_j + 7 * c5_jm1].re +
                                  c5_sigma1.im * c5_b_A[c5_j + 7 * c5_jm1].im);
                  c5_b_A[c5_jp1 + 7 * c5_jm1].im = c5_a * c5_b_A[c5_jp1 + 7 *
                    c5_jm1].im - (c5_sigma1.re * c5_c_A.im - c5_sigma1.im *
                                  c5_d_A.re);
                  c5_b_A[c5_j + 7 * c5_jm1].re = c5_t1.re + c5_d.re;
                  c5_b_A[c5_j + 7 * c5_jm1].im = c5_t1.im + c5_d.im;
                }

                c5_sigma1.re = -c5_sigma1.re;
                c5_sigma1.im = -c5_sigma1.im;
                c5_jm1 = c5_jp1 + 2;
                if (c5_ilast < c5_jm1) {
                  c5_jm1 = c5_ilast;
                }

                for (c5_i = c5_ifrstm - 1; c5_i + 1 <= c5_jm1; c5_i++) {
                  sf_mex_lw_bounds_check(c5_jp1 + 1, 1, 7);
                  sf_mex_lw_bounds_check(c5_i + 1, 1, 7);
                  c5_t1.re = c5_a * c5_b_A[c5_i + 7 * c5_jp1].re;
                  c5_t1.im = c5_a * c5_b_A[c5_i + 7 * c5_jp1].im;
                  c5_d.re = c5_sigma1.re * c5_b_A[c5_i + 7 * c5_j].re -
                    c5_sigma1.im * c5_b_A[c5_i + 7 * c5_j].im;
                  c5_d.im = c5_sigma1.re * c5_b_A[c5_i + 7 * c5_j].im +
                    c5_sigma1.im * c5_b_A[c5_i + 7 * c5_j].re;
                  c5_c_A = c5_b_A[c5_i + 7 * c5_jp1];
                  c5_d_A = c5_b_A[c5_i + 7 * c5_jp1];
                  c5_b_A[c5_i + 7 * c5_j].re = c5_a * c5_b_A[c5_i + 7 * c5_j].re
                    - (c5_sigma1.re * c5_b_A[c5_i + 7 * c5_jp1].re +
                       c5_sigma1.im * c5_b_A[c5_i + 7 * c5_jp1].im);
                  c5_b_A[c5_i + 7 * c5_j].im = c5_a * c5_b_A[c5_i + 7 * c5_j].im
                    - (c5_sigma1.re * c5_c_A.im - c5_sigma1.im * c5_d_A.re);
                  c5_b_A[c5_i + 7 * c5_jp1].re = c5_t1.re + c5_d.re;
                  c5_b_A[c5_i + 7 * c5_jp1].im = c5_t1.im + c5_d.im;
                }

                c5_jm1 = c5_j;
                c5_j = c5_jp1;
              }
            }

            c5_b_guard1 = TRUE;
          }

          if (c5_b_guard1 == TRUE) {
            c5_jiter++;
          }
        }
      } else {
        c5_guard2 = TRUE;
        c5_exitg1 = 1U;
      }
    } while (c5_exitg1 == 0U);
  } else {
    c5_guard1 = TRUE;
  }

  if (c5_guard2 == TRUE) {
    if (c5_failed) {
      *c5_info = (real_T)c5_ilast;
      for (c5_jm1 = 1; c5_jm1 <= c5_ilast; c5_jm1++) {
        c5_alpha1[sf_mex_lw_bounds_check(c5_jm1, 1, 7) - 1] = c5_dc1;
        c5_beta1[c5_jm1 - 1] = c5_dc1;
      }
    } else {
      c5_guard1 = TRUE;
    }
  }

  if (c5_guard1 == TRUE) {
    for (c5_j = 1; c5_j <= c5_ilo - 1; c5_j++) {
      c5_alpha1[sf_mex_lw_bounds_check(c5_j, 1, 7) - 1] = c5_b_A
        [(sf_mex_lw_bounds_check(c5_j, 1, 7) + 7 * (sf_mex_lw_bounds_check(c5_j,
            1, 7) - 1)) - 1];
    }

    *c5_info = 0.0;
  }
}

static real_T c5_eml_matlab_zlanhs(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_A[49], int32_T c5_ilo, int32_T c5_ihi)
{
  real_T c5_f;
  real_T c5_scale;
  real_T c5_sumsq;
  boolean_T c5_firstNonZero;
  int32_T c5_j;
  int32_T c5_c;
  int32_T c5_i;
  real_T c5_temp1;
  real_T c5_temp2;
  c5_f = 0.0;
  if (!(c5_ilo > c5_ihi)) {
    c5_scale = 0.0;
    c5_sumsq = 0.0;
    c5_firstNonZero = TRUE;
    for (c5_j = c5_ilo; c5_j <= c5_ihi; c5_j++) {
      c5_c = c5_j + 1;
      if (c5_ihi < c5_c) {
        c5_c = c5_ihi;
      }

      for (c5_i = c5_ilo; c5_i <= c5_c; c5_i++) {
        sf_mex_lw_bounds_check(c5_i, 1, 7);
        sf_mex_lw_bounds_check(c5_j, 1, 7);
        if (c5_A[(c5_i + 7 * (c5_j - 1)) - 1].re != 0.0) {
          c5_temp1 = muDoubleScalarAbs(c5_A[(c5_i + 7 * (c5_j - 1)) - 1].re);
          if (c5_firstNonZero) {
            c5_sumsq = 1.0;
            c5_scale = c5_temp1;
            c5_firstNonZero = FALSE;
          } else if (c5_scale < c5_temp1) {
            c5_temp2 = c5_scale / c5_temp1;
            c5_sumsq = 1.0 + c5_sumsq * c5_temp2 * c5_temp2;
            c5_scale = c5_temp1;
          } else {
            c5_temp2 = c5_temp1 / c5_scale;
            c5_sumsq += c5_temp2 * c5_temp2;
          }
        }

        if (c5_A[(c5_i + 7 * (c5_j - 1)) - 1].im != 0.0) {
          c5_temp1 = muDoubleScalarAbs(c5_A[(c5_i + 7 * (c5_j - 1)) - 1].im);
          if (c5_firstNonZero) {
            c5_sumsq = 1.0;
            c5_scale = c5_temp1;
            c5_firstNonZero = FALSE;
          } else if (c5_scale < c5_temp1) {
            c5_temp2 = c5_scale / c5_temp1;
            c5_sumsq = 1.0 + c5_sumsq * c5_temp2 * c5_temp2;
            c5_scale = c5_temp1;
          } else {
            c5_temp2 = c5_temp1 / c5_scale;
            c5_sumsq += c5_temp2 * c5_temp2;
          }
        }
      }
    }

    if (c5_sumsq < 0.0) {
      c5_eml_error(chartInstance);
    }

    c5_f = c5_scale * muDoubleScalarSqrt(c5_sumsq);
  }

  return c5_f;
}

static creal_T c5_eml_div(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_x, real_T c5_y)
{
  creal_T c5_z;
  if (c5_x.im == 0.0) {
    c5_z.re = c5_x.re / c5_y;
    c5_z.im = 0.0;
  } else if (c5_x.re == 0.0) {
    c5_z.re = 0.0;
    c5_z.im = c5_x.im / c5_y;
  } else {
    c5_z.re = c5_x.re / c5_y;
    c5_z.im = c5_x.im / c5_y;
  }

  return c5_z;
}

static void c5_b_eml_matlab_zlartg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn)
{
  real_T c5_scale;
  real_T c5_b;
  real_T c5_f2;
  real_T c5_fs_re;
  real_T c5_fs_im;
  real_T c5_gs_re;
  real_T c5_gs_im;
  boolean_T c5_guard1 = FALSE;
  real_T c5_b_b;
  real_T c5_g2s;
  c5_scale = muDoubleScalarAbs(c5_f.re);
  c5_b = muDoubleScalarAbs(c5_f.im);
  if (c5_b > c5_scale) {
    c5_scale = c5_b;
  }

  c5_b = muDoubleScalarAbs(c5_g.re);
  c5_f2 = muDoubleScalarAbs(c5_g.im);
  if (c5_f2 > c5_b) {
    c5_b = c5_f2;
  }

  if (c5_b > c5_scale) {
    c5_scale = c5_b;
  }

  c5_fs_re = c5_f.re;
  c5_fs_im = c5_f.im;
  c5_gs_re = c5_g.re;
  c5_gs_im = c5_g.im;
  c5_guard1 = FALSE;
  if (c5_scale >= 7.4428285367870146E+137) {
    do {
      c5_fs_re *= 1.3435752215134178E-138;
      c5_fs_im *= 1.3435752215134178E-138;
      c5_gs_re *= 1.3435752215134178E-138;
      c5_gs_im *= 1.3435752215134178E-138;
      c5_scale *= 1.3435752215134178E-138;
    } while (!(c5_scale < 7.4428285367870146E+137));

    c5_guard1 = TRUE;
  } else if (c5_scale <= 1.3435752215134178E-138) {
    if ((c5_g.re == 0.0) && (c5_g.im == 0.0)) {
      *c5_cs = 1.0;
      c5_sn->re = 0.0;
      c5_sn->im = 0.0;
    } else {
      do {
        c5_fs_re *= 7.4428285367870146E+137;
        c5_fs_im *= 7.4428285367870146E+137;
        c5_gs_re *= 7.4428285367870146E+137;
        c5_gs_im *= 7.4428285367870146E+137;
        c5_scale *= 7.4428285367870146E+137;
      } while (!(c5_scale > 1.3435752215134178E-138));

      c5_guard1 = TRUE;
    }
  } else {
    c5_guard1 = TRUE;
  }

  if (c5_guard1 == TRUE) {
    c5_f2 = c5_fs_re * c5_fs_re + c5_fs_im * c5_fs_im;
    c5_scale = c5_gs_re * c5_gs_re + c5_gs_im * c5_gs_im;
    c5_b = c5_scale;
    if (1.0 > c5_scale) {
      c5_b = 1.0;
    }

    if (c5_f2 <= c5_b * 2.0041683600089728E-292) {
      if ((c5_f.re == 0.0) && (c5_f.im == 0.0)) {
        *c5_cs = 0.0;
        c5_f2 = muDoubleScalarAbs(c5_gs_re);
        c5_b_b = muDoubleScalarAbs(c5_gs_im);
        if (c5_f2 < c5_b_b) {
          c5_f2 /= c5_b_b;
          c5_b_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
        } else if (c5_f2 > c5_b_b) {
          c5_b_b /= c5_f2;
          c5_b_b = muDoubleScalarSqrt(c5_b_b * c5_b_b + 1.0) * c5_f2;
        } else {
          if (!muDoubleScalarIsNaN(c5_b_b)) {
            c5_b_b = c5_f2 * 1.4142135623730951;
          }
        }

        c5_sn->re = c5_gs_re / c5_b_b;
        c5_sn->im = -c5_gs_im / c5_b_b;
      } else {
        c5_f2 = muDoubleScalarAbs(c5_fs_re);
        c5_b = muDoubleScalarAbs(c5_fs_im);
        if (c5_f2 < c5_b) {
          c5_f2 /= c5_b;
          c5_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
        } else if (c5_f2 > c5_b) {
          c5_b /= c5_f2;
          c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_f2;
        } else {
          if (!muDoubleScalarIsNaN(c5_b)) {
            c5_b = c5_f2 * 1.4142135623730951;
          }
        }

        if (c5_scale < 0.0) {
          c5_eml_error(chartInstance);
        }

        c5_g2s = muDoubleScalarSqrt(c5_scale);
        *c5_cs = c5_b / c5_g2s;
        c5_b = muDoubleScalarAbs(c5_f.re);
        c5_f2 = muDoubleScalarAbs(c5_f.im);
        if (c5_f2 > c5_b) {
          c5_b = c5_f2;
        }

        if (c5_b > 1.0) {
          c5_f2 = muDoubleScalarAbs(c5_f.re);
          c5_b_b = muDoubleScalarAbs(c5_f.im);
          if (c5_f2 < c5_b_b) {
            c5_f2 /= c5_b_b;
            c5_b_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
          } else if (c5_f2 > c5_b_b) {
            c5_b_b /= c5_f2;
            c5_b_b = muDoubleScalarSqrt(c5_b_b * c5_b_b + 1.0) * c5_f2;
          } else {
            if (!muDoubleScalarIsNaN(c5_b_b)) {
              c5_b_b = c5_f2 * 1.4142135623730951;
            }
          }

          c5_fs_re = c5_f.re / c5_b_b;
          c5_fs_im = c5_f.im / c5_b_b;
        } else {
          c5_b = 7.4428285367870146E+137 * c5_f.re;
          c5_scale = 7.4428285367870146E+137 * c5_f.im;
          c5_f2 = muDoubleScalarAbs(c5_b);
          c5_b_b = muDoubleScalarAbs(c5_scale);
          if (c5_f2 < c5_b_b) {
            c5_f2 /= c5_b_b;
            c5_b_b *= muDoubleScalarSqrt(c5_f2 * c5_f2 + 1.0);
          } else if (c5_f2 > c5_b_b) {
            c5_b_b /= c5_f2;
            c5_b_b = muDoubleScalarSqrt(c5_b_b * c5_b_b + 1.0) * c5_f2;
          } else {
            if (!muDoubleScalarIsNaN(c5_b_b)) {
              c5_b_b = c5_f2 * 1.4142135623730951;
            }
          }

          c5_fs_re = c5_b / c5_b_b;
          c5_fs_im = c5_scale / c5_b_b;
        }

        c5_gs_re /= c5_g2s;
        c5_gs_im = -c5_gs_im / c5_g2s;
        c5_sn->re = c5_fs_re * c5_gs_re - c5_fs_im * c5_gs_im;
        c5_sn->im = c5_fs_re * c5_gs_im + c5_fs_im * c5_gs_re;
      }
    } else {
      c5_b = 1.0 + c5_scale / c5_f2;
      if (c5_b < 0.0) {
        c5_eml_error(chartInstance);
      }

      c5_b = muDoubleScalarSqrt(c5_b);
      *c5_cs = 1.0 / c5_b;
      c5_b_b = c5_f2 + c5_scale;
      c5_fs_re = c5_b * c5_fs_re / c5_b_b;
      c5_fs_im = c5_b * c5_fs_im / c5_b_b;
      c5_sn->re = c5_fs_re * c5_gs_re - c5_fs_im * -c5_gs_im;
      c5_sn->im = c5_fs_re * -c5_gs_im + c5_fs_im * c5_gs_re;
    }
  }
}

static void c5_b_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i8;
  static char_T c5_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c5_u[26];
  const mxArray *c5_y = NULL;
  for (c5_i8 = 0; c5_i8 < 26; c5_i8++) {
    c5_u[c5_i8] = c5_varargin_1[c5_i8];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 26), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static void c5_c_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i9;
  static char_T c5_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c5_u[34];
  const mxArray *c5_y = NULL;
  for (c5_i9 = 0; c5_i9 < 34; c5_i9++) {
    c5_u[c5_i9] = c5_varargin_1[c5_i9];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 34), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static void c5_d_eye(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T c5_I[49])
{
  int32_T c5_i;
  int32_T c5_b_i;
  for (c5_i = 0; c5_i < 49; c5_i++) {
    c5_I[c5_i] = 0.0;
  }

  c5_i = 0;
  for (c5_b_i = 0; c5_b_i < 7; c5_b_i++) {
    c5_I[c5_i] = 1.0;
    c5_i += 8;
  }
}

static void c5_mldivide(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  creal_T c5_A[49], real_T c5_B[7], creal_T c5_Y[7])
{
  int32_T c5_i10;
  creal_T c5_b_A[49];
  int8_T c5_ipiv[7];
  int32_T c5_info;
  int32_T c5_j;
  int32_T c5_jj;
  int32_T c5_jp1j;
  int32_T c5_iy;
  int32_T c5_ix;
  real_T c5_smax;
  int32_T c5_b_k;
  real_T c5_s;
  int32_T c5_jrow;
  creal_T c5_temp;
  int32_T c5_jA;
  for (c5_i10 = 0; c5_i10 < 49; c5_i10++) {
    c5_b_A[c5_i10] = c5_A[c5_i10];
  }

  for (c5_i10 = 0; c5_i10 < 7; c5_i10++) {
    c5_ipiv[c5_i10] = (int8_T)(1 + c5_i10);
  }

  c5_info = 0;
  for (c5_j = 0; c5_j < 6; c5_j++) {
    c5_jj = c5_j << 3;
    c5_jp1j = c5_jj + 2;
    c5_iy = 1;
    c5_ix = c5_jj;
    c5_smax = muDoubleScalarAbs(c5_b_A[c5_jj].re) + muDoubleScalarAbs
      (c5_b_A[c5_jj].im);
    for (c5_b_k = 2; c5_b_k <= 7 - c5_j; c5_b_k++) {
      c5_ix++;
      sf_mex_lw_bounds_check(c5_ix + 1, 1, 49);
      c5_s = muDoubleScalarAbs(c5_b_A[c5_ix].re) + muDoubleScalarAbs
        (c5_b_A[c5_ix].im);
      if (c5_s > c5_smax) {
        c5_iy = c5_b_k;
        c5_smax = c5_s;
      }
    }

    if ((c5_b_A[(c5_jj + c5_iy) - 1].re != 0.0) || (c5_b_A[(c5_jj + c5_iy) - 1].
         im != 0.0)) {
      if (c5_iy - 1 != 0) {
        c5_ipiv[c5_j] = (int8_T)(c5_j + c5_iy);
        c5_jrow = 1 + c5_j;
        c5_iy = (c5_jrow + c5_iy) - 1;
        for (c5_b_k = 0; c5_b_k < 7; c5_b_k++) {
          c5_temp = c5_b_A[sf_mex_lw_bounds_check(c5_jrow, 1, 49) - 1];
          c5_b_A[c5_jrow - 1] = c5_b_A[sf_mex_lw_bounds_check(c5_iy, 1, 49) - 1];
          c5_b_A[c5_iy - 1] = c5_temp;
          c5_jrow += 7;
          c5_iy += 7;
        }
      }

      c5_i10 = (c5_jp1j - c5_j) + 5;
      for (c5_jrow = c5_jp1j; c5_jrow <= c5_i10; c5_jrow++) {
        c5_b_A[c5_jrow - 1] = c5_b_eml_div(chartInstance, c5_b_A[c5_jrow - 1],
          c5_b_A[c5_jj]);
      }
    } else {
      c5_info = c5_j + 1;
    }

    c5_jA = c5_jj + 8;
    c5_iy = c5_jj + 7;
    for (c5_jrow = 1; c5_jrow <= 6 - c5_j; c5_jrow++) {
      sf_mex_lw_bounds_check(c5_iy + 1, 1, 49);
      if ((c5_b_A[c5_iy].re != 0.0) || (c5_b_A[c5_iy].im != 0.0)) {
        c5_temp.re = -c5_b_A[c5_iy].re - c5_b_A[c5_iy].im * 0.0;
        c5_temp.im = c5_b_A[c5_iy].re * 0.0 + -c5_b_A[c5_iy].im;
        c5_ix = c5_jp1j;
        c5_i10 = (c5_jA - c5_j) + 6;
        for (c5_b_k = 1 + c5_jA; c5_b_k <= c5_i10; c5_b_k++) {
          c5_smax = c5_b_A[sf_mex_lw_bounds_check(c5_ix, 1, 49) - 1].re *
            c5_temp.re - c5_b_A[sf_mex_lw_bounds_check(c5_ix, 1, 49) - 1].im *
            c5_temp.im;
          c5_s = c5_b_A[sf_mex_lw_bounds_check(c5_ix, 1, 49) - 1].re *
            c5_temp.im + c5_b_A[sf_mex_lw_bounds_check(c5_ix, 1, 49) - 1].im *
            c5_temp.re;
          c5_b_A[sf_mex_lw_bounds_check(c5_b_k, 1, 49) - 1].re =
            c5_b_A[sf_mex_lw_bounds_check(c5_b_k, 1, 49) - 1].re + c5_smax;
          c5_b_A[sf_mex_lw_bounds_check(c5_b_k, 1, 49) - 1].im =
            c5_b_A[sf_mex_lw_bounds_check(c5_b_k, 1, 49) - 1].im + c5_s;
          c5_ix++;
        }
      }

      c5_iy += 7;
      c5_jA += 7;
    }
  }

  if ((c5_info == 0) && (!((c5_b_A[48].re != 0.0) || (c5_b_A[48].im != 0.0)))) {
    c5_info = 7;
  }

  if (c5_info > 0) {
    c5_d_eml_warning(chartInstance);
  }

  for (c5_i10 = 0; c5_i10 < 7; c5_i10++) {
    c5_Y[c5_i10].re = c5_B[c5_i10];
    c5_Y[c5_i10].im = 0.0;
  }

  for (c5_jrow = 0; c5_jrow < 7; c5_jrow++) {
    if (c5_ipiv[c5_jrow] != c5_jrow + 1) {
      c5_temp = c5_Y[c5_jrow];
      c5_Y[c5_jrow] = c5_Y[sf_mex_lw_bounds_check((int32_T)c5_ipiv[c5_jrow], 1,
        7) - 1];
      c5_Y[sf_mex_lw_bounds_check((int32_T)c5_ipiv[c5_jrow], 1, 7) - 1] =
        c5_temp;
    }
  }

  for (c5_b_k = 0; c5_b_k < 7; c5_b_k++) {
    c5_iy = 7 * c5_b_k;
    if ((c5_Y[c5_b_k].re != 0.0) || (c5_Y[c5_b_k].im != 0.0)) {
      for (c5_jrow = c5_b_k + 2; c5_jrow < 8; c5_jrow++) {
        c5_smax = c5_Y[c5_b_k].re * c5_b_A[(c5_jrow + c5_iy) - 1].im +
          c5_Y[c5_b_k].im * c5_b_A[(c5_jrow + c5_iy) - 1].re;
        c5_Y[c5_jrow - 1].re -= c5_Y[c5_b_k].re * c5_b_A[(c5_jrow + c5_iy) - 1].
          re - c5_Y[c5_b_k].im * c5_b_A[(c5_jrow + c5_iy) - 1].im;
        c5_Y[c5_jrow - 1].im -= c5_smax;
      }
    }
  }

  for (c5_b_k = 6; c5_b_k > -1; c5_b_k += -1) {
    c5_iy = 7 * c5_b_k;
    if ((c5_Y[c5_b_k].re != 0.0) || (c5_Y[c5_b_k].im != 0.0)) {
      c5_Y[c5_b_k] = c5_b_eml_div(chartInstance, c5_Y[c5_b_k], c5_b_A[c5_b_k +
        c5_iy]);
      for (c5_jrow = 0; c5_jrow + 1 <= c5_b_k; c5_jrow++) {
        c5_smax = c5_Y[c5_b_k].re * c5_b_A[c5_jrow + c5_iy].im + c5_Y[c5_b_k].im
          * c5_b_A[c5_jrow + c5_iy].re;
        c5_Y[c5_jrow].re -= c5_Y[c5_b_k].re * c5_b_A[c5_jrow + c5_iy].re -
          c5_Y[c5_b_k].im * c5_b_A[c5_jrow + c5_iy].im;
        c5_Y[c5_jrow].im -= c5_smax;
      }
    }
  }
}

static creal_T c5_b_eml_div(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, creal_T c5_x, creal_T c5_y)
{
  creal_T c5_z;
  real_T c5_brm;
  real_T c5_bim;
  real_T c5_d;
  if (c5_y.im == 0.0) {
    if (c5_x.im == 0.0) {
      c5_z.re = c5_x.re / c5_y.re;
      c5_z.im = 0.0;
    } else if (c5_x.re == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = c5_x.im / c5_y.re;
    } else {
      c5_z.re = c5_x.re / c5_y.re;
      c5_z.im = c5_x.im / c5_y.re;
    }
  } else if (c5_y.re == 0.0) {
    if (c5_x.re == 0.0) {
      c5_z.re = c5_x.im / c5_y.im;
      c5_z.im = 0.0;
    } else if (c5_x.im == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = -(c5_x.re / c5_y.im);
    } else {
      c5_z.re = c5_x.im / c5_y.im;
      c5_z.im = -(c5_x.re / c5_y.im);
    }
  } else {
    c5_brm = muDoubleScalarAbs(c5_y.re);
    c5_bim = muDoubleScalarAbs(c5_y.im);
    if (c5_brm > c5_bim) {
      c5_bim = c5_y.im / c5_y.re;
      c5_d = c5_y.re + c5_bim * c5_y.im;
      c5_z.re = (c5_x.re + c5_bim * c5_x.im) / c5_d;
      c5_z.im = (c5_x.im - c5_bim * c5_x.re) / c5_d;
    } else if (c5_bim == c5_brm) {
      c5_bim = c5_y.re > 0.0 ? 0.5 : -0.5;
      c5_d = c5_y.im > 0.0 ? 0.5 : -0.5;
      c5_z.re = (c5_x.re * c5_bim + c5_x.im * c5_d) / c5_brm;
      c5_z.im = (c5_x.im * c5_bim - c5_x.re * c5_d) / c5_brm;
    } else {
      c5_bim = c5_y.re / c5_y.im;
      c5_d = c5_y.im + c5_bim * c5_y.re;
      c5_z.re = (c5_bim * c5_x.re + c5_x.im) / c5_d;
      c5_z.im = (c5_bim * c5_x.im - c5_x.re) / c5_d;
    }
  }

  return c5_z;
}

static void c5_d_eml_warning(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i11;
  static char_T c5_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c5_u[27];
  const mxArray *c5_y = NULL;
  for (c5_i11 = 0; c5_i11 < 27; c5_i11++) {
    c5_u[c5_i11] = c5_varargin_1[c5_i11];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call("warning", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static void c5_abs(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   creal_T c5_x[2], real_T c5_y[2])
{
  int32_T c5_b_k;
  real_T c5_a;
  real_T c5_b;
  for (c5_b_k = 0; c5_b_k < 2; c5_b_k++) {
    c5_a = muDoubleScalarAbs(c5_x[c5_b_k].re);
    c5_b = muDoubleScalarAbs(c5_x[c5_b_k].im);
    if (c5_a < c5_b) {
      c5_a /= c5_b;
      c5_b *= muDoubleScalarSqrt(c5_a * c5_a + 1.0);
    } else if (c5_a > c5_b) {
      c5_b /= c5_a;
      c5_b = muDoubleScalarSqrt(c5_b * c5_b + 1.0) * c5_a;
    } else {
      if (!muDoubleScalarIsNaN(c5_b)) {
        c5_b = c5_a * 1.4142135623730951;
      }
    }

    c5_y[c5_b_k] = c5_b;
  }
}

static void c5_b_eml_error(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
  int32_T c5_i12;
  static char_T c5_varargin_1[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'l', 'o', 'g', '1', '0', '_', 'd', 'o', 'm',
    'a', 'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c5_u[31];
  const mxArray *c5_y = NULL;
  for (c5_i12 = 0; c5_i12 < 31; c5_i12++) {
    c5_u[c5_i12] = c5_varargin_1[c5_i12];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 31), FALSE);
  sf_mex_call("error", 0U, 1U, 14, sf_mex_call("message", 1U, 1U, 14, c5_y));
}

static real_T c5_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_Fs1, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_Fs1), &c5_thisId);
  sf_mex_destroy(&c5_Fs1);
  return c5_y;
}

static real_T c5_b_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d1;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d1, 1, 0, 0U, 0, 0U, 0);
  c5_y = c5_d1;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_c_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_Wn, const char_T *c5_identifier, real_T
  c5_y[7])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_Wn), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_Wn);
}

static void c5_d_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7])
{
  real_T c5_dv19[7];
  int32_T c5_i13;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv19, 1, 0, 0U, 1, 0U, 1, 7);
  for (c5_i13 = 0; c5_i13 < 7; c5_i13++) {
    c5_y[c5_i13] = c5_dv19[c5_i13];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_e_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_ABK, const char_T *c5_identifier, real_T
  c5_y[70])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_ABK), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_ABK);
}

static void c5_f_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[70])
{
  real_T c5_dv20[70];
  int32_T c5_i14;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_ABK_not_empty = FALSE;
  } else {
    chartInstance->c5_ABK_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv20, 1, 0, 0U, 1, 0U, 2, 7,
                  10);
    for (c5_i14 = 0; c5_i14 < 70; c5_i14++) {
      c5_y[c5_i14] = c5_dv20[c5_i14];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_g_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_CD, const char_T *c5_identifier, real_T
  c5_y[16])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_CD), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_CD);
}

static void c5_h_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[16])
{
  real_T c5_dv21[16];
  int32_T c5_i15;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_CD_not_empty = FALSE;
  } else {
    chartInstance->c5_CD_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv21, 1, 0, 0U, 1, 0U, 2, 2,
                  8);
    for (c5_i15 = 0; c5_i15 < 16; c5_i15++) {
      c5_y[c5_i15] = c5_dv21[c5_i15];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_i_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_E1, const char_T *c5_identifier, real_T
  c5_y[2])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_E1), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_E1);
}

static void c5_j_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[2])
{
  real_T c5_dv22[2];
  int32_T c5_i16;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_E1_not_empty = FALSE;
  } else {
    chartInstance->c5_E1_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv22, 1, 0, 0U, 1, 0U, 1, 2);
    for (c5_i16 = 0; c5_i16 < 2; c5_i16++) {
      c5_y[c5_i16] = c5_dv22[c5_i16];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_k_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_LK, const char_T *c5_identifier, real_T
  c5_y[25000])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_LK), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_LK);
}

static void c5_l_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[25000])
{
  static real_T c5_dv23[25000];
  int32_T c5_i17;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_LK_not_empty = FALSE;
  } else {
    chartInstance->c5_LK_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv23, 1, 0, 0U, 1, 0U, 2,
                  100, 250);
    for (c5_i17 = 0; c5_i17 < 25000; c5_i17++) {
      c5_y[c5_i17] = c5_dv23[c5_i17];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_m_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Pabk, const char_T *c5_identifier, real_T
  c5_y[100])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_Pabk), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_Pabk);
}

static void c5_n_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[100])
{
  real_T c5_dv24[100];
  int32_T c5_i18;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_Pabk_not_empty = FALSE;
  } else {
    chartInstance->c5_Pabk_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv24, 1, 0, 0U, 1, 0U, 2, 10,
                  10);
    for (c5_i18 = 0; c5_i18 < 100; c5_i18++) {
      c5_y[c5_i18] = c5_dv24[c5_i18];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_o_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Pcd, const char_T *c5_identifier, real_T
  c5_y[64])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_Pcd), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_Pcd);
}

static void c5_p_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[64])
{
  real_T c5_dv25[64];
  int32_T c5_i19;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_Pcd_not_empty = FALSE;
  } else {
    chartInstance->c5_Pcd_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv25, 1, 0, 0U, 1, 0U, 2, 8,
                  8);
    for (c5_i19 = 0; c5_i19 < 64; c5_i19++) {
      c5_y[c5_i19] = c5_dv25[c5_i19];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_q_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Plk, const char_T *c5_identifier, real_T
  c5_y[63504])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_Plk), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_Plk);
}

static void c5_r_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[63504])
{
  static real_T c5_dv26[63504];
  int32_T c5_i20;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_Plk_not_empty = FALSE;
  } else {
    chartInstance->c5_Plk_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv26, 1, 0, 0U, 1, 0U, 2,
                  252, 252);
    for (c5_i20 = 0; c5_i20 < 63504; c5_i20++) {
      c5_y[c5_i20] = c5_dv26[c5_i20];
    }
  }

  sf_mex_destroy(&c5_u);
}

static real_T c5_s_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_U1, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_U1), &c5_thisId);
  sf_mex_destroy(&c5_b_U1);
  return c5_y;
}

static real_T c5_t_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d2;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_U1_not_empty = FALSE;
  } else {
    chartInstance->c5_U1_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d2, 1, 0, 0U, 0, 0U, 0);
    c5_y = c5_d2;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_u_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_VARX, const char_T *c5_identifier, real_T
  c5_y[504])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_v_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_VARX), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_VARX);
}

static void c5_v_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[504])
{
  real_T c5_dv27[504];
  int32_T c5_i21;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_VARX_not_empty = FALSE;
  } else {
    chartInstance->c5_VARX_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv27, 1, 0, 0U, 1, 0U, 2, 2,
                  252);
    for (c5_i21 = 0; c5_i21 < 504; c5_i21++) {
      c5_y[c5_i21] = c5_dv27[c5_i21];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_w_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_X, const char_T *c5_identifier, real_T
  c5_y[7])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_x_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_X), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_X);
}

static void c5_x_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7])
{
  real_T c5_dv28[7];
  int32_T c5_i22;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_X_not_empty = FALSE;
  } else {
    chartInstance->c5_X_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv28, 1, 0, 0U, 1, 0U, 1, 7);
    for (c5_i22 = 0; c5_i22 < 7; c5_i22++) {
      c5_y[c5_i22] = c5_dv28[c5_i22];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_y_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_X1, const char_T *c5_identifier, real_T
  c5_y[7])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_ab_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_X1), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_X1);
}

static void c5_ab_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[7])
{
  real_T c5_dv29[7];
  int32_T c5_i23;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_X1_not_empty = FALSE;
  } else {
    chartInstance->c5_X1_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv29, 1, 0, 0U, 1, 0U, 1, 7);
    for (c5_i23 = 0; c5_i23 < 7; c5_i23++) {
      c5_y[c5_i23] = c5_dv29[c5_i23];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_bb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Y1, const char_T *c5_identifier, real_T
  c5_y[2])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_cb_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_Y1), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_Y1);
}

static void c5_cb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[2])
{
  real_T c5_dv30[2];
  int32_T c5_i24;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_Y1_not_empty = FALSE;
  } else {
    chartInstance->c5_Y1_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv30, 1, 0, 0U, 1, 0U, 1, 2);
    for (c5_i24 = 0; c5_i24 < 2; c5_i24++) {
      c5_y[c5_i24] = c5_dv30[c5_i24];
    }
  }

  sf_mex_destroy(&c5_u);
}

static void c5_db_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_Zp, const char_T *c5_identifier, real_T
  c5_y[250])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_eb_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_Zp), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_Zp);
}

static void c5_eb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[250])
{
  real_T c5_dv31[250];
  int32_T c5_i25;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_Zp_not_empty = FALSE;
  } else {
    chartInstance->c5_Zp_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv31, 1, 0, 0U, 1, 0U, 1,
                  250);
    for (c5_i25 = 0; c5_i25 < 250; c5_i25++) {
      c5_y[c5_i25] = c5_dv31[c5_i25];
    }
  }

  sf_mex_destroy(&c5_u);
}

static real_T c5_fb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_k, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_gb_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_k), &c5_thisId);
  sf_mex_destroy(&c5_b_k);
  return c5_y;
}

static real_T c5_gb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d3;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_k_not_empty = FALSE;
  } else {
    chartInstance->c5_k_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d3, 1, 0, 0U, 0, 0U, 0);
    c5_y = c5_d3;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static real_T c5_hb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_saw, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_ib_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_saw), &c5_thisId);
  sf_mex_destroy(&c5_b_saw);
  return c5_y;
}

static real_T c5_ib_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d4;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_saw_not_empty = FALSE;
  } else {
    chartInstance->c5_saw_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d4, 1, 0, 0U, 0, 0U, 0);
    c5_y = c5_d4;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static real_T c5_jb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_start, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_kb_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_start),
    &c5_thisId);
  sf_mex_destroy(&c5_b_start);
  return c5_y;
}

static real_T c5_kb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d5;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_start_not_empty = FALSE;
  } else {
    chartInstance->c5_start_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d5, 1, 0, 0U, 0, 0U, 0);
    c5_y = c5_d5;
  }

  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_lb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_b_w, const char_T *c5_identifier, real_T
  c5_y[300])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_mb_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_w), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_b_w);
}

static void c5_mb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId,
  real_T c5_y[300])
{
  real_T c5_dv32[300];
  int32_T c5_i26;
  if (mxIsEmpty(c5_u)) {
    chartInstance->c5_w_not_empty = FALSE;
  } else {
    chartInstance->c5_w_not_empty = TRUE;
    sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv32, 1, 0, 0U, 1, 0U, 2, 1,
                  300);
    for (c5_i26 = 0; c5_i26 < 300; c5_i26++) {
      c5_y[c5_i26] = c5_dv32[c5_i26];
    }
  }

  sf_mex_destroy(&c5_u);
}

static uint8_T c5_nb_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct *
  chartInstance, const mxArray *c5_b_is_active_c5_sim1_rlti_varmax_gust2, const
  char_T *c5_identifier)
{
  uint8_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_ob_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c5_b_is_active_c5_sim1_rlti_varmax_gust2), &c5_thisId);
  sf_mex_destroy(&c5_b_is_active_c5_sim1_rlti_varmax_gust2);
  return c5_y;
}

static uint8_T c5_ob_emlrt_marshallIn(SFc5_sim1_rlti_varmax_gust2InstanceStruct *
  chartInstance, const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint8_T c5_y;
  uint8_T c5_u0;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u0, 1, 3, 0U, 0, 0U, 0);
  c5_y = c5_u0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_rls_ew_reg(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, real_T c5_z[252], real_T c5_y[2], real_T c5_theta[504], real_T
  c5_P[63504], real_T c5_lambda, real_T c5_reg, real_T *c5_b_saw)
{
  int32_T c5_i27;
  real_T c5_eta[252];
  real_T c5_b_y;
  int32_T c5_m;
  int32_T c5_n;
  static real_T c5_c_y[63504];
  int32_T c5_b_k;
  real_T c5_alpha1;
  int32_T c5_lda;
  int32_T c5_ldb;
  real_T c5_beta1;
  int32_T c5_ldc;
  char_T c5_TRANSA;
  char_T c5_TRANSB;
  static real_T c5_d_y[63504];
  real_T c5_a[252];
  real_T c5_e_y[252];
  real_T c5_f_y[2];
  real_T c5_g_y[504];
  real_T c5_h_y[504];
  for (c5_i27 = 0; c5_i27 < 252; c5_i27++) {
    c5_eta[c5_i27] = 0.0;
  }

  if (*c5_b_saw > 251.5) {
    *c5_b_saw = 1.0;
  } else {
    (*c5_b_saw)++;
  }

  c5_b_y = (1.0 - muDoubleScalarPower(c5_lambda, 252.0)) * c5_reg;
  if (c5_b_y < 0.0) {
    c5_eml_error(chartInstance);
  }

  c5_eta[sf_mex_lw_bounds_check((int32_T)*c5_b_saw, 1, 252) - 1] =
    muDoubleScalarSqrt(c5_b_y);
  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_c_y[c5_n + c5_i27] = c5_eta[c5_n] * c5_eta[c5_m];
    }

    c5_i27 += 252;
  }

  c5_m = 252;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 252;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 252;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 63504; c5_i27++) {
    c5_d_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_P[0],
          &c5_lda, &c5_c_y[0], &c5_ldb, &c5_beta1, &c5_d_y[0], &c5_ldc);
  c5_m = 252;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 252;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 252;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 63504; c5_i27++) {
    c5_c_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_d_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_c_y[0], &c5_ldc);
  c5_m = 1;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 1;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 1;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 252; c5_i27++) {
    c5_a[c5_i27] = c5_eta[c5_i27];
    c5_e_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_a[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_e_y[0], &c5_ldc);
  c5_b_y = 0.0;
  for (c5_b_k = 0; c5_b_k < 252; c5_b_k++) {
    c5_b_y += c5_e_y[c5_b_k] * c5_eta[c5_b_k];
  }

  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_P[c5_n + c5_i27] -= c5_c_y[c5_n + c5_i27] / (1.0 + c5_b_y);
    }

    c5_i27 += 252;
  }

  c5_b_y = 1.0 / c5_lambda;
  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_c_y[c5_n + c5_i27] = c5_z[c5_n] * c5_z[c5_m];
    }

    c5_i27 += 252;
  }

  c5_m = 252;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 252;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 252;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 63504; c5_i27++) {
    c5_d_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_P[0],
          &c5_lda, &c5_c_y[0], &c5_ldb, &c5_beta1, &c5_d_y[0], &c5_ldc);
  c5_m = 252;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 252;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 252;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 63504; c5_i27++) {
    c5_c_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_d_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_c_y[0], &c5_ldc);
  c5_m = 1;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 1;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 1;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 252; c5_i27++) {
    c5_a[c5_i27] = c5_z[c5_i27];
    c5_e_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_a[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_e_y[0], &c5_ldc);
  c5_alpha1 = 0.0;
  for (c5_b_k = 0; c5_b_k < 252; c5_b_k++) {
    c5_alpha1 += c5_e_y[c5_b_k] * c5_z[c5_b_k];
  }

  c5_alpha1 += c5_lambda;
  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_P[c5_n + c5_i27] = c5_b_y * (c5_P[c5_n + c5_i27] - c5_c_y[c5_n + c5_i27]
        / c5_alpha1);
    }

    c5_i27 += 252;
  }

  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    c5_n = 0;
    for (c5_b_k = 0; c5_b_k < 252; c5_b_k++) {
      c5_c_y[c5_b_k + c5_i27] = 0.5 * (c5_P[c5_b_k + c5_i27] + c5_P[c5_n + c5_m]);
      c5_n += 252;
    }

    c5_i27 += 252;
  }

  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_P[c5_n + c5_i27] = c5_c_y[c5_n + c5_i27];
    }

    c5_i27 += 252;
  }

  for (c5_i27 = 0; c5_i27 < 2; c5_i27++) {
    c5_alpha1 = 0.0;
    c5_m = 0;
    for (c5_n = 0; c5_n < 252; c5_n++) {
      c5_alpha1 += c5_theta[c5_m + c5_i27] * c5_z[c5_n];
      c5_m += 2;
    }

    c5_f_y[c5_i27] = c5_y[c5_i27] - c5_alpha1;
  }

  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 2; c5_n++) {
      c5_g_y[c5_n + c5_i27] = c5_f_y[c5_n] * c5_z[c5_m];
    }

    c5_i27 += 2;
  }

  c5_m = 2;
  c5_n = 252;
  c5_b_k = 252;
  c5_alpha1 = 1.0;
  c5_lda = 2;
  c5_ldb = 252;
  c5_beta1 = 0.0;
  c5_ldc = 2;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i27 = 0; c5_i27 < 504; c5_i27++) {
    c5_h_y[c5_i27] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_g_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_h_y[0], &c5_ldc);
  c5_i27 = 0;
  for (c5_m = 0; c5_m < 252; c5_m++) {
    for (c5_n = 0; c5_n < 2; c5_n++) {
      c5_theta[c5_n + c5_i27] += c5_h_y[c5_n + c5_i27];
    }

    c5_i27 += 2;
  }
}

static void c5_rls_ew(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                      real_T c5_z[8], real_T c5_y[2], real_T c5_theta[16],
                      real_T c5_P[64], real_T c5_lambda)
{
  real_T c5_b_y;
  int32_T c5_i28;
  int32_T c5_m;
  int32_T c5_n;
  real_T c5_c_y[64];
  int32_T c5_b_k;
  real_T c5_alpha1;
  int32_T c5_lda;
  int32_T c5_ldb;
  real_T c5_beta1;
  int32_T c5_ldc;
  char_T c5_TRANSA;
  char_T c5_TRANSB;
  real_T c5_d_y[64];
  real_T c5_e_y[8];
  real_T c5_f_y[2];
  real_T c5_g_y[16];
  real_T c5_b_theta[16];
  c5_b_y = 1.0 / c5_lambda;
  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    for (c5_n = 0; c5_n < 8; c5_n++) {
      c5_c_y[c5_n + c5_i28] = c5_z[c5_n] * c5_z[c5_m];
    }

    c5_i28 += 8;
  }

  c5_m = 8;
  c5_n = 8;
  c5_b_k = 8;
  c5_alpha1 = 1.0;
  c5_lda = 8;
  c5_ldb = 8;
  c5_beta1 = 0.0;
  c5_ldc = 8;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i28 = 0; c5_i28 < 64; c5_i28++) {
    c5_d_y[c5_i28] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_P[0],
          &c5_lda, &c5_c_y[0], &c5_ldb, &c5_beta1, &c5_d_y[0], &c5_ldc);
  c5_m = 8;
  c5_n = 8;
  c5_b_k = 8;
  c5_alpha1 = 1.0;
  c5_lda = 8;
  c5_ldb = 8;
  c5_beta1 = 0.0;
  c5_ldc = 8;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i28 = 0; c5_i28 < 64; c5_i28++) {
    c5_c_y[c5_i28] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_d_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_c_y[0], &c5_ldc);
  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    c5_e_y[c5_m] = 0.0;
    for (c5_n = 0; c5_n < 8; c5_n++) {
      c5_alpha1 = c5_e_y[c5_m] + c5_z[c5_n] * c5_P[c5_n + c5_i28];
      c5_e_y[c5_m] = c5_alpha1;
    }

    c5_i28 += 8;
  }

  c5_alpha1 = 0.0;
  for (c5_b_k = 0; c5_b_k < 8; c5_b_k++) {
    c5_alpha1 += c5_e_y[c5_b_k] * c5_z[c5_b_k];
  }

  c5_alpha1 += c5_lambda;
  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    for (c5_n = 0; c5_n < 8; c5_n++) {
      c5_P[c5_n + c5_i28] = c5_b_y * (c5_P[c5_n + c5_i28] - c5_c_y[c5_n + c5_i28]
        / c5_alpha1);
    }

    c5_i28 += 8;
  }

  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    c5_n = 0;
    for (c5_b_k = 0; c5_b_k < 8; c5_b_k++) {
      c5_c_y[c5_b_k + c5_i28] = 0.5 * (c5_P[c5_b_k + c5_i28] + c5_P[c5_n + c5_m]);
      c5_n += 8;
    }

    c5_i28 += 8;
  }

  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    for (c5_n = 0; c5_n < 8; c5_n++) {
      c5_P[c5_n + c5_i28] = c5_c_y[c5_n + c5_i28];
    }

    c5_i28 += 8;
  }

  for (c5_i28 = 0; c5_i28 < 2; c5_i28++) {
    c5_alpha1 = 0.0;
    c5_m = 0;
    for (c5_n = 0; c5_n < 8; c5_n++) {
      c5_alpha1 += c5_theta[c5_m + c5_i28] * c5_z[c5_n];
      c5_m += 2;
    }

    c5_f_y[c5_i28] = c5_y[c5_i28] - c5_alpha1;
  }

  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    for (c5_n = 0; c5_n < 2; c5_n++) {
      c5_g_y[c5_n + c5_i28] = c5_f_y[c5_n] * c5_z[c5_m];
    }

    c5_i28 += 2;
  }

  for (c5_i28 = 0; c5_i28 < 2; c5_i28++) {
    c5_m = 0;
    c5_n = 0;
    for (c5_b_k = 0; c5_b_k < 8; c5_b_k++) {
      c5_alpha1 = 0.0;
      c5_lda = 0;
      for (c5_ldb = 0; c5_ldb < 8; c5_ldb++) {
        c5_alpha1 += c5_g_y[c5_lda + c5_i28] * c5_P[c5_ldb + c5_n];
        c5_lda += 2;
      }

      c5_b_theta[c5_m + c5_i28] = c5_theta[c5_m + c5_i28] + c5_alpha1;
      c5_m += 2;
      c5_n += 8;
    }
  }

  c5_i28 = 0;
  for (c5_m = 0; c5_m < 8; c5_m++) {
    for (c5_n = 0; c5_n < 2; c5_n++) {
      c5_theta[c5_n + c5_i28] = c5_b_theta[c5_n + c5_i28];
    }

    c5_i28 += 2;
  }
}

static void c5_b_rls_ew(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
  real_T c5_z[10], real_T c5_y[7], real_T c5_theta[70], real_T c5_P[100], real_T
  c5_lambda)
{
  real_T c5_b_y;
  int32_T c5_i29;
  int32_T c5_m;
  int32_T c5_n;
  real_T c5_c_y[100];
  int32_T c5_b_k;
  real_T c5_alpha1;
  int32_T c5_lda;
  int32_T c5_ldb;
  real_T c5_beta1;
  int32_T c5_ldc;
  char_T c5_TRANSA;
  char_T c5_TRANSB;
  real_T c5_d_y[100];
  real_T c5_e_y[10];
  real_T c5_f_y[7];
  real_T c5_g_y[70];
  real_T c5_h_y[70];
  c5_b_y = 1.0 / c5_lambda;
  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    for (c5_n = 0; c5_n < 10; c5_n++) {
      c5_c_y[c5_n + c5_i29] = c5_z[c5_n] * c5_z[c5_m];
    }

    c5_i29 += 10;
  }

  c5_m = 10;
  c5_n = 10;
  c5_b_k = 10;
  c5_alpha1 = 1.0;
  c5_lda = 10;
  c5_ldb = 10;
  c5_beta1 = 0.0;
  c5_ldc = 10;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i29 = 0; c5_i29 < 100; c5_i29++) {
    c5_d_y[c5_i29] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_P[0],
          &c5_lda, &c5_c_y[0], &c5_ldb, &c5_beta1, &c5_d_y[0], &c5_ldc);
  c5_m = 10;
  c5_n = 10;
  c5_b_k = 10;
  c5_alpha1 = 1.0;
  c5_lda = 10;
  c5_ldb = 10;
  c5_beta1 = 0.0;
  c5_ldc = 10;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i29 = 0; c5_i29 < 100; c5_i29++) {
    c5_c_y[c5_i29] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_d_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_c_y[0], &c5_ldc);
  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    c5_e_y[c5_m] = 0.0;
    for (c5_n = 0; c5_n < 10; c5_n++) {
      c5_alpha1 = c5_e_y[c5_m] + c5_z[c5_n] * c5_P[c5_n + c5_i29];
      c5_e_y[c5_m] = c5_alpha1;
    }

    c5_i29 += 10;
  }

  c5_alpha1 = 0.0;
  for (c5_b_k = 0; c5_b_k < 10; c5_b_k++) {
    c5_alpha1 += c5_e_y[c5_b_k] * c5_z[c5_b_k];
  }

  c5_alpha1 += c5_lambda;
  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    for (c5_n = 0; c5_n < 10; c5_n++) {
      c5_P[c5_n + c5_i29] = c5_b_y * (c5_P[c5_n + c5_i29] - c5_c_y[c5_n + c5_i29]
        / c5_alpha1);
    }

    c5_i29 += 10;
  }

  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    c5_n = 0;
    for (c5_b_k = 0; c5_b_k < 10; c5_b_k++) {
      c5_c_y[c5_b_k + c5_i29] = 0.5 * (c5_P[c5_b_k + c5_i29] + c5_P[c5_n + c5_m]);
      c5_n += 10;
    }

    c5_i29 += 10;
  }

  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    for (c5_n = 0; c5_n < 10; c5_n++) {
      c5_P[c5_n + c5_i29] = c5_c_y[c5_n + c5_i29];
    }

    c5_i29 += 10;
  }

  for (c5_i29 = 0; c5_i29 < 7; c5_i29++) {
    c5_alpha1 = 0.0;
    c5_m = 0;
    for (c5_n = 0; c5_n < 10; c5_n++) {
      c5_alpha1 += c5_theta[c5_m + c5_i29] * c5_z[c5_n];
      c5_m += 7;
    }

    c5_f_y[c5_i29] = c5_y[c5_i29] - c5_alpha1;
  }

  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    for (c5_n = 0; c5_n < 7; c5_n++) {
      c5_g_y[c5_n + c5_i29] = c5_f_y[c5_n] * c5_z[c5_m];
    }

    c5_i29 += 7;
  }

  c5_m = 7;
  c5_n = 10;
  c5_b_k = 10;
  c5_alpha1 = 1.0;
  c5_lda = 7;
  c5_ldb = 10;
  c5_beta1 = 0.0;
  c5_ldc = 7;
  c5_TRANSA = 'N';
  c5_TRANSB = 'N';
  for (c5_i29 = 0; c5_i29 < 70; c5_i29++) {
    c5_h_y[c5_i29] = 0.0;
  }

  dgemm32(&c5_TRANSA, &c5_TRANSB, &c5_m, &c5_n, &c5_b_k, &c5_alpha1, &c5_g_y[0],
          &c5_lda, &c5_P[0], &c5_ldb, &c5_beta1, &c5_h_y[0], &c5_ldc);
  c5_i29 = 0;
  for (c5_m = 0; c5_m < 10; c5_m++) {
    for (c5_n = 0; c5_n < 7; c5_n++) {
      c5_theta[c5_n + c5_i29] += c5_h_y[c5_n + c5_i29];
    }

    c5_i29 += 7;
  }
}

static void c5_sqrt(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                    creal_T *c5_x)
{
  real_T c5_absxi;
  real_T c5_a;
  real_T c5_absxr;
  if (c5_x->im == 0.0) {
    if (c5_x->re < 0.0) {
      c5_absxi = 0.0;
      c5_a = muDoubleScalarSqrt(muDoubleScalarAbs(c5_x->re));
    } else {
      c5_absxi = muDoubleScalarSqrt(c5_x->re);
      c5_a = 0.0;
    }
  } else if (c5_x->re == 0.0) {
    if (c5_x->im < 0.0) {
      c5_absxi = muDoubleScalarSqrt(-c5_x->im / 2.0);
      c5_a = -c5_absxi;
    } else {
      c5_absxi = muDoubleScalarSqrt(c5_x->im / 2.0);
      c5_a = c5_absxi;
    }
  } else if (muDoubleScalarIsNaN(c5_x->re) || muDoubleScalarIsNaN(c5_x->im)) {
    c5_absxi = rtNaN;
    c5_a = rtNaN;
  } else if (muDoubleScalarIsInf(c5_x->im)) {
    c5_absxi = rtInf;
    c5_a = c5_x->im;
  } else if (muDoubleScalarIsInf(c5_x->re)) {
    if (c5_x->re < 0.0) {
      c5_absxi = 0.0;
      c5_a = rtInf;
    } else {
      c5_absxi = rtInf;
      c5_a = 0.0;
    }
  } else {
    c5_absxr = muDoubleScalarAbs(c5_x->re);
    c5_absxi = muDoubleScalarAbs(c5_x->im);
    if ((c5_absxr > 4.4942328371557893E+307) || (c5_absxi >
         4.4942328371557893E+307)) {
      c5_absxr *= 0.5;
      c5_absxi *= 0.5;
      if (c5_absxr < c5_absxi) {
        c5_a = c5_absxr / c5_absxi;
        c5_absxi *= muDoubleScalarSqrt(c5_a * c5_a + 1.0);
      } else if (c5_absxr > c5_absxi) {
        c5_absxi /= c5_absxr;
        c5_absxi = muDoubleScalarSqrt(c5_absxi * c5_absxi + 1.0) * c5_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c5_absxi)) {
          c5_absxi = c5_absxr * 1.4142135623730951;
        }
      }

      if (c5_absxi > c5_absxr) {
        c5_absxi = muDoubleScalarSqrt(c5_absxi) * muDoubleScalarSqrt(1.0 +
          c5_absxr / c5_absxi);
      } else {
        c5_absxi = muDoubleScalarSqrt(c5_absxi) * 1.4142135623730951;
      }
    } else {
      if (c5_absxr < c5_absxi) {
        c5_a = c5_absxr / c5_absxi;
        c5_absxi *= muDoubleScalarSqrt(c5_a * c5_a + 1.0);
      } else if (c5_absxr > c5_absxi) {
        c5_absxi /= c5_absxr;
        c5_absxi = muDoubleScalarSqrt(c5_absxi * c5_absxi + 1.0) * c5_absxr;
      } else {
        if (!muDoubleScalarIsNaN(c5_absxi)) {
          c5_absxi = c5_absxr * 1.4142135623730951;
        }
      }

      c5_absxi = muDoubleScalarSqrt((c5_absxi + c5_absxr) * 0.5);
    }

    if (c5_x->re > 0.0) {
      c5_a = 0.5 * (c5_x->im / c5_absxi);
    } else {
      if (c5_x->im < 0.0) {
        c5_a = -c5_absxi;
      } else {
        c5_a = c5_absxi;
      }

      c5_absxi = 0.5 * (c5_x->im / c5_a);
    }
  }

  c5_x->re = c5_absxi;
  c5_x->im = c5_a;
}

static void c5_exp(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                   creal_T *c5_x)
{
  real_T c5_r;
  creal_T c5_b_x;
  creal_T c5_c_x;
  c5_r = muDoubleScalarExp(c5_x->re / 2.0);
  c5_b_x = *c5_x;
  c5_c_x = *c5_x;
  c5_x->re = c5_r * (c5_r * muDoubleScalarCos(c5_b_x.im));
  c5_x->im = c5_r * (c5_r * muDoubleScalarSin(c5_c_x.im));
}

static void c5_log10(SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance,
                     real_T *c5_x)
{
  if (*c5_x < 0.0) {
    c5_b_eml_error(chartInstance);
  }

  *c5_x = muDoubleScalarLog10(*c5_x);
}

static int32_T c5_div_s32_floor(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance, int32_T c5_numerator, int32_T c5_denominator)
{
  int32_T c5_quotient;
  uint32_T c5_absNumerator;
  uint32_T c5_absDenominator;
  int32_T c5_quotientNeedsNegation;
  uint32_T c5_tempAbsQuotient;
  if (c5_denominator == 0) {
    c5_quotient = c5_numerator >= 0 ? MAX_int32_T : MIN_int32_T;
    sf_mex_dividebyzero_error();
  } else {
    c5_absNumerator = (uint32_T)(c5_numerator >= 0 ? c5_numerator :
      -c5_numerator);
    c5_absDenominator = (uint32_T)(c5_denominator >= 0 ? c5_denominator :
      -c5_denominator);
    c5_quotientNeedsNegation = (c5_numerator < 0 != c5_denominator < 0);
    c5_tempAbsQuotient = c5_absNumerator / c5_absDenominator;
    if ((uint32_T)c5_quotientNeedsNegation) {
      c5_absNumerator %= c5_absDenominator;
      if (c5_absNumerator > (uint32_T)0) {
        c5_tempAbsQuotient++;
      }
    }

    c5_quotient = (uint32_T)c5_quotientNeedsNegation ? -(int32_T)
      c5_tempAbsQuotient : (int32_T)c5_tempAbsQuotient;
  }

  return c5_quotient;
}

static void init_dsm_address_info(SFc5_sim1_rlti_varmax_gust2InstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
void sf_c5_sim1_rlti_varmax_gust2_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(425620736U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1059273904U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3747300354U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1930482429U);
}

mxArray *sf_c5_sim1_rlti_varmax_gust2_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("NBzCRZCbppGV1MgVolS4tG");
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

static const mxArray *sf_get_sim_state_info_c5_sim1_rlti_varmax_gust2(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x10'type','srcId','name','auxInfo'{{M[1],M[5],T\"Fs1\",},{M[1],M[6],T\"Fs2\",},{M[1],M[9],T\"Wn\",},{M[1],M[4],T\"Ws\",},{M[1],M[10],T\"Zn\",},{M[4],M[0],T\"ABK\",S'l','i','p'{{M1x2[1210 1213],M[0],}}},{M[4],M[0],T\"CD\",S'l','i','p'{{M1x2[1214 1216],M[0],}}},{M[4],M[0],T\"E1\",S'l','i','p'{{M1x2[1228 1230],M[0],}}},{M[4],M[0],T\"LK\",S'l','i','p'{{M1x2[1198 1200],M[0],}}},{M[4],M[0],T\"Pabk\",S'l','i','p'{{M1x2[1205 1209],M[0],}}}}",
    "100 S1x10'type','srcId','name','auxInfo'{{M[4],M[0],T\"Pcd\",S'l','i','p'{{M1x2[1201 1204],M[0],}}},{M[4],M[0],T\"Plk\",S'l','i','p'{{M1x2[1189 1192],M[0],}}},{M[4],M[0],T\"U1\",S'l','i','p'{{M1x2[1217 1219],M[0],}}},{M[4],M[0],T\"VARX\",S'l','i','p'{{M1x2[1193 1197],M[0],}}},{M[4],M[0],T\"X\",S'l','i','p'{{M1x2[1223 1224],M[0],}}},{M[4],M[0],T\"X1\",S'l','i','p'{{M1x2[1225 1227],M[0],}}},{M[4],M[0],T\"Y1\",S'l','i','p'{{M1x2[1220 1222],M[0],}}},{M[4],M[0],T\"Zp\",S'l','i','p'{{M1x2[1231 1233],M[0],}}},{M[4],M[0],T\"k\",S'l','i','p'{{M1x2[1236 1237],M[0],}}},{M[4],M[0],T\"saw\",S'l','i','p'{{M1x2[1244 1247],M[0],}}}}",
    "100 S1x3'type','srcId','name','auxInfo'{{M[4],M[0],T\"start\",S'l','i','p'{{M1x2[1238 1243],M[0],}}},{M[4],M[0],T\"w\",S'l','i','p'{{M1x2[1234 1235],M[0],}}},{M[8],M[0],T\"is_active_c5_sim1_rlti_varmax_gust2\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 23, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c5_sim1_rlti_varmax_gust2_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void sf_opaque_initialize_c5_sim1_rlti_varmax_gust2(void
  *chartInstanceVar)
{
  initialize_params_c5_sim1_rlti_varmax_gust2
    ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*) chartInstanceVar);
  initialize_c5_sim1_rlti_varmax_gust2
    ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c5_sim1_rlti_varmax_gust2(void *chartInstanceVar)
{
  enable_c5_sim1_rlti_varmax_gust2((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c5_sim1_rlti_varmax_gust2(void *chartInstanceVar)
{
  disable_c5_sim1_rlti_varmax_gust2((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c5_sim1_rlti_varmax_gust2(void *chartInstanceVar)
{
  sf_c5_sim1_rlti_varmax_gust2((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c5_sim1_rlti_varmax_gust2
  (SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c5_sim1_rlti_varmax_gust2
    ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_sim1_rlti_varmax_gust2();/* state var info */
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

extern void sf_internal_set_sim_state_c5_sim1_rlti_varmax_gust2(SimStruct* S,
  const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_sim1_rlti_varmax_gust2();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c5_sim1_rlti_varmax_gust2
    ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)chartInfo->chartInstance,
     mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c5_sim1_rlti_varmax_gust2
  (SimStruct* S)
{
  return sf_internal_get_sim_state_c5_sim1_rlti_varmax_gust2(S);
}

static void sf_opaque_set_sim_state_c5_sim1_rlti_varmax_gust2(SimStruct* S,
  const mxArray *st)
{
  sf_internal_set_sim_state_c5_sim1_rlti_varmax_gust2(S, st);
}

static void sf_opaque_terminate_c5_sim1_rlti_varmax_gust2(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)
                    chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
    }

    finalize_c5_sim1_rlti_varmax_gust2
      ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*) chartInstanceVar);
    free((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }

  unload_sim1_rlti_varmax_gust2_optimization_info();
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc5_sim1_rlti_varmax_gust2
    ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c5_sim1_rlti_varmax_gust2(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c5_sim1_rlti_varmax_gust2
      ((SFc5_sim1_rlti_varmax_gust2InstanceStruct*)(((ChartInfoStruct *)
         ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c5_sim1_rlti_varmax_gust2(SimStruct *S)
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
    mxArray *infoStruct = load_sim1_rlti_varmax_gust2_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,infoStruct,5);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,infoStruct,5,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,infoStruct,5,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,infoStruct,5,4);
      sf_mark_chart_reusable_outputs(S,infoStruct,5,5);
    }

    sf_set_rtw_dwork_info(S,infoStruct,5);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3237911301U));
  ssSetChecksum1(S,(1643827518U));
  ssSetChecksum2(S,(2208136295U));
  ssSetChecksum3(S,(1907997904U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
}

static void mdlRTW_c5_sim1_rlti_varmax_gust2(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c5_sim1_rlti_varmax_gust2(SimStruct *S)
{
  SFc5_sim1_rlti_varmax_gust2InstanceStruct *chartInstance;
  chartInstance = (SFc5_sim1_rlti_varmax_gust2InstanceStruct *)malloc(sizeof
    (SFc5_sim1_rlti_varmax_gust2InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc5_sim1_rlti_varmax_gust2InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.mdlStart = mdlStart_c5_sim1_rlti_varmax_gust2;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c5_sim1_rlti_varmax_gust2;
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

void c5_sim1_rlti_varmax_gust2_method_dispatcher(SimStruct *S, int_T method,
  void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c5_sim1_rlti_varmax_gust2(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c5_sim1_rlti_varmax_gust2(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c5_sim1_rlti_varmax_gust2(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c5_sim1_rlti_varmax_gust2_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
