#ifndef __c5_sim1_rlti_varmax_gust2_h__
#define __c5_sim1_rlti_varmax_gust2_h__

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
} c5_ResolvedFunctionInfo;

typedef struct {
  SimStruct *S;
  real_T c5_ABK[70];
  real_T c5_CD[16];
  real_T c5_Du;
  real_T c5_Dy[4];
  real_T c5_E1[2];
  real_T c5_LK[25000];
  real_T c5_Pabk[100];
  real_T c5_Pcd[64];
  real_T c5_Plk[63504];
  real_T c5_U1;
  real_T c5_VARX[504];
  real_T c5_W[700];
  real_T c5_X[7];
  real_T c5_X1[7];
  real_T c5_Y1[2];
  real_T c5_Zp[250];
  real_T c5_k;
  real_T c5_saw;
  real_T c5_start;
  real_T c5_w[300];
  boolean_T c5_ABK_not_empty;
  boolean_T c5_CD_not_empty;
  boolean_T c5_E1_not_empty;
  boolean_T c5_LK_not_empty;
  boolean_T c5_Pabk_not_empty;
  boolean_T c5_Pcd_not_empty;
  boolean_T c5_Plk_not_empty;
  boolean_T c5_U1_not_empty;
  boolean_T c5_VARX_not_empty;
  boolean_T c5_X1_not_empty;
  boolean_T c5_X_not_empty;
  boolean_T c5_Y1_not_empty;
  boolean_T c5_Zp_not_empty;
  boolean_T c5_k_not_empty;
  boolean_T c5_saw_not_empty;
  boolean_T c5_start_not_empty;
  boolean_T c5_w_not_empty;
  uint8_T c5_is_active_c5_sim1_rlti_varmax_gust2;
  ChartInfoStruct chartInfo;
} SFc5_sim1_rlti_varmax_gust2InstanceStruct;

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c5_sim1_rlti_varmax_gust2_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c5_sim1_rlti_varmax_gust2_get_check_sum(mxArray *plhs[]);
extern void c5_sim1_rlti_varmax_gust2_method_dispatcher(SimStruct *S, int_T
  method, void *data);

#endif
