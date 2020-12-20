#ifndef __c1_sim1_rlti_varx_fast_gust2_h__
#define __c1_sim1_rlti_varx_fast_gust2_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"

/* Type Definitions */
typedef struct {
  const char * context;
  const char * name;
  const char * dominantType;
  const char * resolved;
  uint32_T fileTimeLo;
  uint32_T fileTimeHi;
  uint32_T mFileTimeLo;
  uint32_T mFileTimeHi;
} c1_ResolvedFunctionInfo;

typedef struct {
  SimStruct *S;
  real_T c1_ABK[70];
  real_T c1_CD[16];
  real_T c1_Du;
  real_T c1_Dy[4];
  real_T c1_Glk[459];
  real_T c1_LK[15000];
  real_T c1_Pabk[100];
  real_T c1_Pcd[64];
  real_T c1_Plk[1099];
  real_T c1_U1;
  real_T c1_VARX[306];
  real_T c1_VARX1[300];
  real_T c1_W[700];
  real_T c1_X[7];
  real_T c1_X1[7];
  real_T c1_Y1[2];
  real_T c1_Ylk[3];
  real_T c1_Z[156];
  real_T c1_eta[468];
  real_T c1_k;
  real_T c1_saw;
  real_T c1_start;
  real_T c1_w[300];
  boolean_T c1_ABK_not_empty;
  boolean_T c1_CD_not_empty;
  boolean_T c1_Glk_not_empty;
  boolean_T c1_LK_not_empty;
  boolean_T c1_Pabk_not_empty;
  boolean_T c1_Pcd_not_empty;
  boolean_T c1_Plk_not_empty;
  boolean_T c1_U1_not_empty;
  boolean_T c1_VARX1_not_empty;
  boolean_T c1_VARX_not_empty;
  boolean_T c1_X1_not_empty;
  boolean_T c1_X_not_empty;
  boolean_T c1_Y1_not_empty;
  boolean_T c1_Ylk_not_empty;
  boolean_T c1_Z_not_empty;
  boolean_T c1_eta_not_empty;
  boolean_T c1_k_not_empty;
  boolean_T c1_saw_not_empty;
  boolean_T c1_start_not_empty;
  boolean_T c1_w_not_empty;
  uint8_T c1_is_active_c1_sim1_rlti_varx_fast_gust2;
  ChartInfoStruct chartInfo;
} SFc1_sim1_rlti_varx_fast_gust2InstanceStruct;

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c1_sim1_rlti_varx_fast_gust2_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_sim1_rlti_varx_fast_gust2_get_check_sum(mxArray *plhs[]);
extern void c1_sim1_rlti_varx_fast_gust2_method_dispatcher(SimStruct *S, int_T
  method, void *data);

#endif
