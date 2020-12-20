#ifndef __c2_sim1_rlti_varx_gust2_h__
#define __c2_sim1_rlti_varx_gust2_h__

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
} c2_ResolvedFunctionInfo;

typedef struct {
  SimStruct *S;
  real_T c2_ABK[70];
  real_T c2_CD[16];
  real_T c2_Du;
  real_T c2_Dy[4];
  real_T c2_LK[15000];
  real_T c2_Pabk[100];
  real_T c2_Pcd[64];
  real_T c2_Plk[23104];
  real_T c2_U1;
  real_T c2_VARX[304];
  real_T c2_W[700];
  real_T c2_X[7];
  real_T c2_X1[7];
  real_T c2_Y1[2];
  real_T c2_Zp[150];
  real_T c2_k;
  real_T c2_saw;
  real_T c2_start;
  real_T c2_w[300];
  boolean_T c2_ABK_not_empty;
  boolean_T c2_CD_not_empty;
  boolean_T c2_LK_not_empty;
  boolean_T c2_Pabk_not_empty;
  boolean_T c2_Pcd_not_empty;
  boolean_T c2_Plk_not_empty;
  boolean_T c2_U1_not_empty;
  boolean_T c2_VARX_not_empty;
  boolean_T c2_X1_not_empty;
  boolean_T c2_X_not_empty;
  boolean_T c2_Y1_not_empty;
  boolean_T c2_Zp_not_empty;
  boolean_T c2_k_not_empty;
  boolean_T c2_saw_not_empty;
  boolean_T c2_start_not_empty;
  boolean_T c2_w_not_empty;
  uint8_T c2_is_active_c2_sim1_rlti_varx_gust2;
  ChartInfoStruct chartInfo;
} SFc2_sim1_rlti_varx_gust2InstanceStruct;

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_sim1_rlti_varx_gust2_get_eml_resolved_functions_info
  (void);

/* Function Definitions */
extern void sf_c2_sim1_rlti_varx_gust2_get_check_sum(mxArray *plhs[]);
extern void c2_sim1_rlti_varx_gust2_method_dispatcher(SimStruct *S, int_T method,
  void *data);

#endif
