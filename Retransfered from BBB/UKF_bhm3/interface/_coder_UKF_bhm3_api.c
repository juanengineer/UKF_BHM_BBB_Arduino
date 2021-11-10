/*
 * File: _coder_UKF_bhm3_api.c
 *
 * MATLAB Coder version            : 3.0
 * C/C++ source code generated on  : 05-Apr-2016 17:59:50
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_UKF_bhm3_api.h"
#include "_coder_UKF_bhm3_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true, false, 131419U, NULL, "UKF_bhm3", NULL,
  false, { 2045744189U, 2170104910U, 2743257031U, 4284093946U }, NULL };

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Param, const
  char_T *identifier, struct0_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *state_ukf, const char_T *identifier))[6];
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *i, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const struct1_T *u);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6];
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *w_ukf,
  const char_T *identifier))[3];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3];
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *P, const
  char_T *identifier))[9];
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[9];
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6];
static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3];
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[9];

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *Param
 *                const char_T *identifier
 *                struct0_T *y
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Param, const
  char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(Param), &thisId, y);
  emlrtDestroyArray(&Param);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct0_T *y
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const int32_T dims = 0;
  static const char * fieldNames[20] = { "qmax", "Cmax", "Ccb0", "Ccb1", "Ccb2",
    "Ccb3", "Rs", "Cs", "Rcp0", "Rcp1", "Rcp2", "Ccp0", "Ccp1", "Ccp2", "Rp",
    "Jt", "hcp", "hcs", "ha", "Ta" };

  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 20, fieldNames, 0U, &dims);
  thisId.fIdentifier = "qmax";
  y->qmax = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "qmax")), &thisId);
  thisId.fIdentifier = "Cmax";
  y->Cmax = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Cmax")), &thisId);
  thisId.fIdentifier = "Ccb0";
  y->Ccb0 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccb0")), &thisId);
  thisId.fIdentifier = "Ccb1";
  y->Ccb1 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccb1")), &thisId);
  thisId.fIdentifier = "Ccb2";
  y->Ccb2 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccb2")), &thisId);
  thisId.fIdentifier = "Ccb3";
  y->Ccb3 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccb3")), &thisId);
  thisId.fIdentifier = "Rs";
  y->Rs = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Rs")),
    &thisId);
  thisId.fIdentifier = "Cs";
  y->Cs = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Cs")),
    &thisId);
  thisId.fIdentifier = "Rcp0";
  y->Rcp0 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Rcp0")), &thisId);
  thisId.fIdentifier = "Rcp1";
  y->Rcp1 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Rcp1")), &thisId);
  thisId.fIdentifier = "Rcp2";
  y->Rcp2 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Rcp2")), &thisId);
  thisId.fIdentifier = "Ccp0";
  y->Ccp0 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccp0")), &thisId);
  thisId.fIdentifier = "Ccp1";
  y->Ccp1 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccp1")), &thisId);
  thisId.fIdentifier = "Ccp2";
  y->Ccp2 = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "Ccp2")), &thisId);
  thisId.fIdentifier = "Rp";
  y->Rp = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Rp")),
    &thisId);
  thisId.fIdentifier = "Jt";
  y->Jt = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Jt")),
    &thisId);
  thisId.fIdentifier = "hcp";
  y->hcp = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "hcp")),
    &thisId);
  thisId.fIdentifier = "hcs";
  y->hcs = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "hcs")),
    &thisId);
  thisId.fIdentifier = "ha";
  y->ha = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "ha")),
    &thisId);
  thisId.fIdentifier = "Ta";
  y->Ta = b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Ta")),
    &thisId);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *state_ukf
 *                const char_T *identifier
 * Return Type  : real_T (*)[6]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *state_ukf, const char_T *identifier))[6]
{
  real_T (*y)[6];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(state_ukf), &thisId);
  emlrtDestroyArray(&state_ukf);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *i
 *                const char_T *identifier
 * Return Type  : real_T
 */
  static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *i, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(i), &thisId);
  emlrtDestroyArray(&i);
  return y;
}

/*
 * Arguments    : const struct1_T *u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const struct1_T *u)
{
  const mxArray *y;
  const mxArray *b_y;
  static const int32_T iv0[2] = { 3, 3 };

  const mxArray *m0;
  real_T *pData;
  int32_T i;
  const mxArray *c_y;
  static const int32_T iv1[2] = { 3, 2 };

  const mxArray *m1;
  real_T *b_pData;
  int32_T b_i;
  const mxArray *d_y;
  static const int32_T iv2[2] = { 3, 2 };

  const mxArray *m2;
  real_T *c_pData;
  int32_T c_i;
  const mxArray *e_y;
  const mxArray *m3;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 0, NULL));
  b_y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  pData = (real_T *)mxGetPr(m0);
  for (i = 0; i < 9; i++) {
    pData[i] = u->P[i];
  }

  emlrtAssign(&b_y, m0);
  emlrtAddField(y, b_y, "P", 0);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  b_pData = (real_T *)mxGetPr(m1);
  for (b_i = 0; b_i < 6; b_i++) {
    b_pData[b_i] = u->state_ukf[b_i];
  }

  emlrtAssign(&c_y, m1);
  emlrtAddField(y, c_y, "state_ukf", 0);
  d_y = NULL;
  m2 = emlrtCreateNumericArray(2, iv2, mxDOUBLE_CLASS, mxREAL);
  c_pData = (real_T *)mxGetPr(m2);
  for (c_i = 0; c_i < 6; c_i++) {
    c_pData[c_i] = u->est_ukf[c_i];
  }

  emlrtAssign(&d_y, m2);
  emlrtAddField(y, d_y, "est_ukf", 0);
  e_y = NULL;
  m3 = emlrtCreateDoubleScalar(u->SOC);
  emlrtAssign(&e_y, m3);
  emlrtAddField(y, e_y, "SOC", 0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[6]
 */
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6]
{
  real_T (*y)[6];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *w_ukf
 *                const char_T *identifier
 * Return Type  : real_T (*)[3]
 */
  static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *w_ukf,
  const char_T *identifier))[3]
{
  real_T (*y)[3];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(w_ukf), &thisId);
  emlrtDestroyArray(&w_ukf);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[3]
 */
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3]
{
  real_T (*y)[3];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *P
 *                const char_T *identifier
 * Return Type  : real_T (*)[9]
 */
  static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *P,
  const char_T *identifier))[9]
{
  real_T (*y)[9];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(P), &thisId);
  emlrtDestroyArray(&P);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[9]
 */
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[9]
{
  real_T (*y)[9];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
  static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[6]
 */
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6]
{
  real_T (*ret)[6];
  static const int32_T dims[2] = { 3, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[6])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[3]
 */
  static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3]
{
  real_T (*ret)[3];
  static const int32_T dims[1] = { 3 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[3])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[9]
 */
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[9]
{
  real_T (*ret)[9];
  static const int32_T dims[2] = { 3, 3 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[9])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const mxArray *prhs[12]
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
  void UKF_bhm3_api(const mxArray *prhs[12], const mxArray *plhs[1])
{
  real_T i;
  struct0_T Param;
  real_T z_ukf;
  real_T (*state_ukf)[6];
  real_T (*w_ukf)[3];
  real_T (*est_ukf)[6];
  real_T (*P)[9];
  real_T (*Q)[9];
  real_T R;
  real_T dt;
  real_T ndim;
  real_T zdim;
  struct1_T ukfout;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  prhs[3] = emlrtProtectR2012b(prhs[3], 3, false, -1);
  prhs[4] = emlrtProtectR2012b(prhs[4], 4, false, -1);
  prhs[5] = emlrtProtectR2012b(prhs[5], 5, false, -1);
  prhs[6] = emlrtProtectR2012b(prhs[6], 6, false, -1);
  prhs[7] = emlrtProtectR2012b(prhs[7], 7, false, -1);

  /* Marshall function inputs */
  i = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "i");
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Param", &Param);
  z_ukf = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "z_ukf");
  state_ukf = e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "state_ukf");
  w_ukf = g_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "w_ukf");
  est_ukf = e_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "est_ukf");
  P = i_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "P");
  Q = i_emlrt_marshallIn(&st, emlrtAlias(prhs[7]), "Q");
  R = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "R");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "dt");
  ndim = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "ndim");
  zdim = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "zdim");

  /* Invoke the target function */
  UKF_bhm3(i, &Param, z_ukf, *state_ukf, *w_ukf, *est_ukf, *P, *Q, R, dt, ndim,
           zdim, &ukfout);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&ukfout);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void UKF_bhm3_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  UKF_bhm3_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void UKF_bhm3_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void UKF_bhm3_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_UKF_bhm3_api.c
 *
 * [EOF]
 */
