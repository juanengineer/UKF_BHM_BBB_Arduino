/*
 * File: _coder_UKF_bhm3_api.h
 *
 * MATLAB Coder version            : 3.0
 * C/C++ source code generated on  : 05-Apr-2016 17:59:50
 */

#ifndef ___CODER_UKF_BHM3_API_H__
#define ___CODER_UKF_BHM3_API_H__

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_UKF_bhm3_api.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T qmax;
  real_T Cmax;
  real_T Ccb0;
  real_T Ccb1;
  real_T Ccb2;
  real_T Ccb3;
  real_T Rs;
  real_T Cs;
  real_T Rcp0;
  real_T Rcp1;
  real_T Rcp2;
  real_T Ccp0;
  real_T Ccp1;
  real_T Ccp2;
  real_T Rp;
  real_T Jt;
  real_T hcp;
  real_T hcs;
  real_T ha;
  real_T Ta;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  real_T P[9];
  real_T state_ukf[6];
  real_T est_ukf[6];
  real_T SOC;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
void UKF_bhm3(real_T i, struct0_T *Param, real_T z_ukf, real_T state_ukf[6],
              real_T w_ukf[3], real_T est_ukf[6], real_T P[9], real_T Q[9],
              real_T R, real_T dt, real_T ndim, real_T zdim, struct1_T *ukfout);
void UKF_bhm3_api(const mxArray *prhs[12], const mxArray *plhs[1]);
void UKF_bhm3_atexit(void);
void UKF_bhm3_initialize(void);
void UKF_bhm3_terminate(void);
void UKF_bhm3_xil_terminate(void);

#endif

/*
 * File trailer for _coder_UKF_bhm3_api.h
 *
 * [EOF]
 */
