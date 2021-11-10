/*
 * File: UKF_bhm3.c
 *
 * MATLAB Coder version            : 3.0
 * C/C++ source code generated on  : 05-Apr-2016 17:59:50
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "UKF_bhm3.h"

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray__common*/

#ifndef typedef_emxArray__common
#define typedef_emxArray__common

typedef struct emxArray__common emxArray__common;

#endif                                 /*typedef_emxArray__common*/

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_int32_T*/

#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T

typedef struct emxArray_int32_T emxArray_int32_T;

#endif                                 /*typedef_emxArray_int32_T*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Function Declarations */
static void b_sqrt(creal_T *x);
static void b_xrot(int n, double x[9], int ix0, int iy0, double c, double s);
static void c_xrot(double x[9], int ix0, int iy0, double c, double s);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static void invNxN(const emxArray_real_T *x, emxArray_real_T *y);
static void mpower(const emxArray_real_T *a, emxArray_real_T *c);
static double rt_hypotd_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);
static void schur(const double A[9], creal_T V[9], creal_T T[9]);
static void sqrtm(const double A[9], creal_T X[9]);
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn);
static void xgehrd(double a[9], double tau[2]);
static void xgerc(int m, int n, double alpha1, int ix0, const double y[3],
                  double A[9], int ia0);
static int xhseqr(double h[9], double z[9]);
static double xnrm2(int n, const double x[3]);
static void xrot(int n, double x[9], int ix0, int iy0, double c, double s);
static double xzlarfg(int n, double *alpha1, double x[3]);

/* Function Definitions */

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
static void b_sqrt(creal_T *x)
{
  double yr;
  double yi;
  double absxr;
  double absxi;
  double absxd2;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      yr = 0.0;
      yi = sqrt(fabs(x->re));
    } else {
      yr = sqrt(x->re);
      yi = 0.0;
    }
  } else if (x->re == 0.0) {
    if (x->im < 0.0) {
      yr = sqrt(-x->im / 2.0);
      yi = -yr;
    } else {
      yr = sqrt(x->im / 2.0);
      yi = yr;
    }
  } else if (rtIsNaN(x->re) || rtIsNaN(x->im)) {
    yr = rtNaN;
    yi = rtNaN;
  } else if (rtIsInf(x->im)) {
    yr = rtInf;
    yi = x->im;
  } else if (rtIsInf(x->re)) {
    if (x->re < 0.0) {
      yr = 0.0;
      yi = rtInf;
    } else {
      yr = rtInf;
      yi = 0.0;
    }
  } else {
    absxr = fabs(x->re);
    absxi = fabs(x->im);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxd2 = rt_hypotd_snf(absxr, absxi);
      if (absxd2 > absxr) {
        yr = sqrt(absxd2) * sqrt(1.0 + absxr / absxd2);
      } else {
        yr = sqrt(absxd2) * 1.4142135623730951;
      }
    } else {
      yr = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (x->re > 0.0) {
      yi = 0.5 * (x->im / yr);
    } else {
      if (x->im < 0.0) {
        yi = -yr;
      } else {
        yi = yr;
      }

      yr = 0.5 * (x->im / yi);
    }
  }

  x->re = yr;
  x->im = yi;
}

/*
 * Arguments    : int n
 *                double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
static void b_xrot(int n, double x[9], int ix0, int iy0, double c, double s)
{
  double temp;
  if (n < 1) {
  } else {
    temp = c * x[ix0 - 1] + s * x[iy0 - 1];
    x[iy0 - 1] = c * x[iy0 - 1] - s * x[ix0 - 1];
    x[ix0 - 1] = temp;
  }
}

/*
 * Arguments    : double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
static void c_xrot(double x[9], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 3; k++) {
    temp = c * x[ix] + s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = temp;
    iy++;
    ix++;
  }
}

/*
 * Arguments    : emxArray__common *emxArray
 *                int oldNumel
 *                int elementSize
 * Return Type  : void
 */
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize)
{
  int newNumel;
  int i;
  int newCapacity;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    newCapacity = emxArray->allocatedSize;
    if (newCapacity < 16) {
      newCapacity = 16;
    }

    while (newCapacity < newNumel) {
      newCapacity = newCapacity * 2;
    }

    newData = calloc((unsigned int)newCapacity, (unsigned int)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = newCapacity;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void invNxN(const emxArray_real_T *x, emxArray_real_T *y)
{
  int n;
  int i129;
  int loop_ub;
  int i130;
  emxArray_real_T *A;
  int i131;
  int b_loop_ub;
  int i132;
  int minval;
  int b_n;
  emxArray_int32_T *ipiv;
  int i133;
  int yk;
  int k;
  int i134;
  int j;
  int mmj;
  int c;
  int idxmax;
  int ix;
  double smax;
  int b_k;
  double s;
  int b_ix;
  int iy;
  int c_k;
  double temp;
  int i135;
  int i;
  int b_c;
  int jA;
  int jy;
  int b_j;
  double yjy;
  int c_ix;
  int i136;
  int ijA;
  int c_n;
  emxArray_int32_T *p;
  int i137;
  int b_yk;
  int d_k;
  int e_k;
  int pipk;
  int f_k;
  int c_c;
  int c_j;
  int b_i;
  int d_j;
  int jBcol;
  int g_k;
  int kAcol;
  double b_x;
  double b_y;
  int c_i;
  n = x->size[0];
  i129 = y->size[0] * y->size[1];
  y->size[0] = x->size[0];
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i129, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i130 = 0; i130 < loop_ub; i130++) {
    y->data[i130] = 0.0;
  }

  emxInit_real_T1(&A, 2);
  i131 = A->size[0] * A->size[1];
  A->size[0] = x->size[0];
  A->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)A, i131, (int)sizeof(double));
  b_loop_ub = x->size[0] * x->size[1];
  for (i132 = 0; i132 < b_loop_ub; i132++) {
    A->data[i132] = x->data[i132];
  }

  minval = x->size[0];
  if (minval < 1) {
    b_n = 0;
  } else {
    b_n = minval;
  }

  emxInit_int32_T1(&ipiv, 2);
  i133 = ipiv->size[0] * ipiv->size[1];
  ipiv->size[0] = 1;
  ipiv->size[1] = b_n;
  emxEnsureCapacity((emxArray__common *)ipiv, i133, (int)sizeof(int));
  if (b_n > 0) {
    ipiv->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      ipiv->data[k - 1] = yk;
    }
  }

  if (x->size[0] < 1) {
  } else {
    if (x->size[0] - 1 <= x->size[0]) {
      i134 = x->size[0] - 1;
    } else {
      i134 = x->size[0];
    }

    for (j = 0; j + 1 <= i134; j++) {
      mmj = n - j;
      c = j * (n + 1);
      if (mmj < 1) {
        idxmax = -1;
      } else {
        idxmax = 0;
        if (mmj > 1) {
          ix = c;
          smax = fabs(A->data[c]);
          for (b_k = 1; b_k + 1 <= mmj; b_k++) {
            ix++;
            s = fabs(A->data[ix]);
            if (s > smax) {
              idxmax = b_k;
              smax = s;
            }
          }
        }
      }

      if (A->data[c + idxmax] != 0.0) {
        if (idxmax != 0) {
          ipiv->data[j] = (j + idxmax) + 1;
          b_ix = j;
          iy = j + idxmax;
          for (c_k = 1; c_k <= n; c_k++) {
            temp = A->data[b_ix];
            A->data[b_ix] = A->data[iy];
            A->data[iy] = temp;
            b_ix += n;
            iy += n;
          }
        }

        i135 = c + mmj;
        for (i = c + 1; i + 1 <= i135; i++) {
          A->data[i] /= A->data[c];
        }
      }

      b_c = (n - j) - 1;
      jA = c + n;
      jy = c + n;
      for (b_j = 1; b_j <= b_c; b_j++) {
        yjy = A->data[jy];
        if (A->data[jy] != 0.0) {
          c_ix = c + 1;
          i136 = mmj + jA;
          for (ijA = 1 + jA; ijA + 1 <= i136; ijA++) {
            A->data[ijA] += A->data[c_ix] * -yjy;
            c_ix++;
          }
        }

        jy += n;
        jA += n;
      }
    }
  }

  if (x->size[0] < 1) {
    c_n = 0;
  } else {
    c_n = x->size[0];
  }

  emxInit_int32_T1(&p, 2);
  i137 = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = c_n;
  emxEnsureCapacity((emxArray__common *)p, i137, (int)sizeof(int));
  if (c_n > 0) {
    p->data[0] = 1;
    b_yk = 1;
    for (d_k = 2; d_k <= c_n; d_k++) {
      b_yk++;
      p->data[d_k - 1] = b_yk;
    }
  }

  for (e_k = 0; e_k < ipiv->size[1]; e_k++) {
    if (ipiv->data[e_k] > 1 + e_k) {
      pipk = p->data[ipiv->data[e_k] - 1];
      p->data[ipiv->data[e_k] - 1] = p->data[e_k];
      p->data[e_k] = pipk;
    }
  }

  emxFree_int32_T(&ipiv);
  for (f_k = 0; f_k + 1 <= n; f_k++) {
    c_c = p->data[f_k] - 1;
    y->data[f_k + y->size[0] * (p->data[f_k] - 1)] = 1.0;
    for (c_j = f_k; c_j + 1 <= n; c_j++) {
      if (y->data[c_j + y->size[0] * c_c] != 0.0) {
        for (b_i = c_j + 1; b_i + 1 <= n; b_i++) {
          y->data[b_i + y->size[0] * c_c] -= y->data[c_j + y->size[0] * c_c] *
            A->data[b_i + A->size[0] * c_j];
        }
      }
    }
  }

  emxFree_int32_T(&p);
  if ((x->size[0] == 0) || ((y->size[0] == 0) || (y->size[1] == 0))) {
  } else {
    for (d_j = 1; d_j <= n; d_j++) {
      jBcol = n * (d_j - 1);
      for (g_k = n - 1; g_k + 1 > 0; g_k--) {
        kAcol = n * g_k;
        if (y->data[g_k + jBcol] != 0.0) {
          b_x = y->data[g_k + jBcol];
          b_y = A->data[g_k + kAcol];
          y->data[g_k + jBcol] = b_x / b_y;
          for (c_i = 0; c_i + 1 <= g_k; c_i++) {
            y->data[c_i + jBcol] -= y->data[g_k + jBcol] * A->data[c_i + kAcol];
          }
        }
      }
    }
  }

  emxFree_real_T(&A);
}

/*
 * Arguments    : const emxArray_real_T *a
 *                emxArray_real_T *c
 * Return Type  : void
 */
static void mpower(const emxArray_real_T *a, emxArray_real_T *c)
{
  int i127;
  int loop_ub;
  int i128;
  if ((a->size[0] == 0) || (a->size[1] == 0)) {
    i127 = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)c, i127, (int)sizeof(double));
    loop_ub = a->size[0] * a->size[1];
    for (i128 = 0; i128 < loop_ub; i128++) {
      c->data[i128] = a->data[i128];
    }
  } else {
    invNxN(a, c);
  }
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d10;
  double d11;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d10 = fabs(u0);
    d11 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d10 == 1.0) {
        y = rtNaN;
      } else if (d10 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d11 == 0.0) {
      y = 1.0;
    } else if (d11 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : const double A[9]
 *                creal_T V[9]
 *                creal_T T[9]
 * Return Type  : void
 */
static void schur(const double A[9], creal_T V[9], creal_T T[9])
{
  double b_A[9];
  double tau[2];
  double Vr[9];
  int j;
  int c;
  int i;
  int b_i;
  double work[3];
  int c_i;
  int itau;
  int d_i;
  int iaii;
  int lastv;
  int e_i;
  int lastc;
  boolean_T exitg2;
  int ia;
  int exitg1;
  int b_lastc;
  int iy;
  int iac;
  int ix;
  double b_c;
  int i125;
  int b_ia;
  int k;
  int b_j;
  int i126;
  int m;
  double a;
  double b;
  double c_c;
  double d;
  double sn;
  double cs;
  double rt2i;
  double rt2r;
  double rt1i;
  double rt1r;
  double mu1_re;
  double r;
  double mu1_im;
  double s;
  int c_j;
  double t1_re;
  double t1_im;
  double T_re;
  double b_mu1_im;
  int f_i;
  double c_mu1_im;
  double b_T_re;
  int g_i;
  double d_mu1_im;
  double V_re;
  memcpy(&b_A[0], &A[0], 9U * sizeof(double));
  xgehrd(b_A, tau);
  memcpy(&Vr[0], &b_A[0], 9U * sizeof(double));
  for (j = 1; j >= 0; j += -1) {
    c = (j + 1) * 3 - 1;
    for (i = 1; i <= j + 1; i++) {
      Vr[c + i] = 0.0;
    }

    b_i = j + 3;
    while (b_i <= 3) {
      Vr[c + 3] = Vr[c];
      b_i = 4;
    }
  }

  for (c_i = 0; c_i < 3; c_i++) {
    Vr[c_i] = 0.0;
    work[c_i] = 0.0;
  }

  Vr[0] = 1.0;
  itau = 1;
  for (d_i = 1; d_i >= 0; d_i += -1) {
    iaii = (d_i + d_i * 3) + 5;
    if (d_i + 1 < 2) {
      Vr[iaii - 1] = 1.0;
      if (tau[itau] != 0.0) {
        lastv = 2;
        e_i = iaii;
        while ((lastv > 0) && (Vr[e_i] == 0.0)) {
          lastv--;
          e_i--;
        }

        lastc = 1;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          ia = iaii + 2;
          do {
            exitg1 = 0;
            if (ia + 1 <= (iaii + lastv) + 2) {
              if (Vr[ia] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc = 0;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        b_lastc = lastc;
      } else {
        lastv = 0;
        b_lastc = 0;
      }

      if (lastv > 0) {
        if (b_lastc == 0) {
        } else {
          work[0] = 0.0;
          iy = 0;
          for (iac = iaii + 3; iac <= iaii + 3; iac += 3) {
            ix = iaii;
            b_c = 0.0;
            i125 = (iac + lastv) - 1;
            for (b_ia = iac; b_ia <= i125; b_ia++) {
              b_c += Vr[b_ia - 1] * Vr[ix - 1];
              ix++;
            }

            work[iy] += b_c;
            iy++;
          }
        }

        xgerc(lastv, b_lastc, -tau[itau], iaii, work, Vr, iaii + 3);
      }
    }

    if (d_i + 1 < 2) {
      for (k = iaii; k + 1 <= iaii + 1; k++) {
        Vr[k] *= -tau[itau];
      }
    }

    Vr[iaii - 1] = 1.0 - tau[itau];
    b_j = 1;
    while (b_j <= d_i) {
      Vr[iaii - 2] = 0.0;
      b_j = 2;
    }

    itau--;
  }

  xhseqr(b_A, Vr);
  for (i126 = 0; i126 < 9; i126++) {
    T[i126].re = b_A[i126];
    T[i126].im = 0.0;
    V[i126].re = Vr[i126];
    V[i126].im = 0.0;
  }

  for (m = 1; m >= 0; m += -1) {
    if (b_A[(m + 3 * m) + 1] != 0.0) {
      a = b_A[m + 3 * m];
      b = b_A[m + 3 * (m + 1)];
      c_c = b_A[(m + 3 * m) + 1];
      d = b_A[(m + 3 * (m + 1)) + 1];
      xdlanv2(&a, &b, &c_c, &d, &rt1r, &rt1i, &rt2r, &rt2i, &cs, &sn);
      mu1_re = rt1r - b_A[(m + 3 * (m + 1)) + 1];
      r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), b_A[(m + 3 * m) + 1]);
      if (rt1i == 0.0) {
        mu1_re /= r;
        mu1_im = 0.0;
      } else if (mu1_re == 0.0) {
        mu1_re = 0.0;
        mu1_im = rt1i / r;
      } else {
        mu1_re /= r;
        mu1_im = rt1i / r;
      }

      s = b_A[(m + 3 * m) + 1] / r;
      for (c_j = m; c_j + 1 < 4; c_j++) {
        t1_re = T[m + 3 * c_j].re;
        t1_im = T[m + 3 * c_j].im;
        T_re = T[m + 3 * c_j].re;
        T[m + 3 * c_j].re = (mu1_re * T[m + 3 * c_j].re + mu1_im * T[m + 3 * c_j]
                             .im) + s * T[(m + 3 * c_j) + 1].re;
        T[m + 3 * c_j].im = (mu1_re * T[m + 3 * c_j].im - mu1_im * T_re) + s *
          T[(m + 3 * c_j) + 1].im;
        b_mu1_im = mu1_re * T[(m + 3 * c_j) + 1].im + mu1_im * T[(m + 3 * c_j) +
          1].re;
        T[(m + 3 * c_j) + 1].re = (mu1_re * T[(m + 3 * c_j) + 1].re - mu1_im *
          T[(m + 3 * c_j) + 1].im) - s * t1_re;
        T[(m + 3 * c_j) + 1].im = b_mu1_im - s * t1_im;
      }

      for (f_i = 0; f_i + 1 <= m + 2; f_i++) {
        t1_re = T[f_i + 3 * m].re;
        t1_im = T[f_i + 3 * m].im;
        c_mu1_im = mu1_re * T[f_i + 3 * m].im + mu1_im * T[f_i + 3 * m].re;
        T[f_i + 3 * m].re = (mu1_re * T[f_i + 3 * m].re - mu1_im * T[f_i + 3 * m]
                             .im) + s * T[f_i + 3 * (m + 1)].re;
        T[f_i + 3 * m].im = c_mu1_im + s * T[f_i + 3 * (m + 1)].im;
        b_T_re = T[f_i + 3 * (m + 1)].re;
        T[f_i + 3 * (m + 1)].re = (mu1_re * T[f_i + 3 * (m + 1)].re + mu1_im *
          T[f_i + 3 * (m + 1)].im) - s * t1_re;
        T[f_i + 3 * (m + 1)].im = (mu1_re * T[f_i + 3 * (m + 1)].im - mu1_im *
          b_T_re) - s * t1_im;
      }

      for (g_i = 0; g_i < 3; g_i++) {
        t1_re = V[g_i + 3 * m].re;
        t1_im = V[g_i + 3 * m].im;
        d_mu1_im = mu1_re * V[g_i + 3 * m].im + mu1_im * V[g_i + 3 * m].re;
        V[g_i + 3 * m].re = (mu1_re * V[g_i + 3 * m].re - mu1_im * V[g_i + 3 * m]
                             .im) + s * V[g_i + 3 * (m + 1)].re;
        V[g_i + 3 * m].im = d_mu1_im + s * V[g_i + 3 * (m + 1)].im;
        V_re = V[g_i + 3 * (m + 1)].re;
        V[g_i + 3 * (m + 1)].re = (mu1_re * V[g_i + 3 * (m + 1)].re + mu1_im *
          V[g_i + 3 * (m + 1)].im) - s * t1_re;
        V[g_i + 3 * (m + 1)].im = (mu1_re * V[g_i + 3 * (m + 1)].im - mu1_im *
          V_re) - s * t1_im;
      }

      T[(m + 3 * m) + 1].re = 0.0;
      T[(m + 3 * m) + 1].im = 0.0;
    }
  }
}

/*
 * Arguments    : const double A[9]
 *                creal_T X[9]
 * Return Type  : void
 */
static void sqrtm(const double A[9], creal_T X[9])
{
  creal_T T[9];
  creal_T Q[9];
  creal_T R[9];
  int i118;
  int j;
  int exitg4;
  int i;
  int exitg3;
  boolean_T p;
  int b_j;
  int c_j;
  int b_i;
  double s_re;
  double s_im;
  int k;
  double R_re;
  double R_im;
  double ar;
  double ai;
  double brm;
  double bim;
  double s;
  double d;
  double sgnbr;
  double sgnbi;
  creal_T b_Q[9];
  int i119;
  int i120;
  int i121;
  int i122;
  int i123;
  double Q_re;
  double Q_im;
  double x[9];
  int i124;
  double y;
  int d_j;
  boolean_T exitg2;
  double b_s;
  int c_i;
  double b_y;
  int e_j;
  boolean_T exitg1;
  double c_s;
  int d_i;
  int f_j;
  int e_i;
  schur(A, Q, T);
  for (i118 = 0; i118 < 9; i118++) {
    R[i118].re = 0.0;
    R[i118].im = 0.0;
  }

  j = 0;
  do {
    exitg4 = 0;
    if (j + 1 < 4) {
      i = 1;
      do {
        exitg3 = 0;
        if (i <= j) {
          if ((T[(i + 3 * j) - 1].re != 0.0) || (T[(i + 3 * j) - 1].im != 0.0))
          {
            p = false;
            exitg3 = 1;
          } else {
            i++;
          }
        } else {
          j++;
          exitg3 = 2;
        }
      } while (exitg3 == 0);

      if (exitg3 == 1) {
        exitg4 = 1;
      }
    } else {
      p = true;
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  if (p) {
    for (b_j = 0; b_j < 3; b_j++) {
      R[b_j + 3 * b_j] = T[b_j + 3 * b_j];
      b_sqrt(&R[b_j + 3 * b_j]);
    }
  } else {
    for (c_j = 0; c_j < 3; c_j++) {
      R[c_j + 3 * c_j] = T[c_j + 3 * c_j];
      b_sqrt(&R[c_j + 3 * c_j]);
      for (b_i = c_j - 1; b_i + 1 > 0; b_i--) {
        s_re = 0.0;
        s_im = 0.0;
        k = b_i + 2;
        while (k <= c_j) {
          s_re += R[3 + b_i].re * R[1 + 3 * c_j].re - R[3 + b_i].im * R[1 + 3 *
            c_j].im;
          s_im += R[3 + b_i].re * R[1 + 3 * c_j].im + R[3 + b_i].im * R[1 + 3 *
            c_j].re;
          k = 3;
        }

        R_re = R[b_i + 3 * b_i].re + R[c_j + 3 * c_j].re;
        R_im = R[b_i + 3 * b_i].im + R[c_j + 3 * c_j].im;
        ar = T[b_i + 3 * c_j].re - s_re;
        ai = T[b_i + 3 * c_j].im - s_im;
        if (R_im == 0.0) {
          if (ai == 0.0) {
            R[b_i + 3 * c_j].re = ar / R_re;
            R[b_i + 3 * c_j].im = 0.0;
          } else if (ar == 0.0) {
            R[b_i + 3 * c_j].re = 0.0;
            R[b_i + 3 * c_j].im = ai / R_re;
          } else {
            R[b_i + 3 * c_j].re = ar / R_re;
            R[b_i + 3 * c_j].im = ai / R_re;
          }
        } else if (R_re == 0.0) {
          if (ar == 0.0) {
            R[b_i + 3 * c_j].re = ai / R_im;
            R[b_i + 3 * c_j].im = 0.0;
          } else if (ai == 0.0) {
            R[b_i + 3 * c_j].re = 0.0;
            R[b_i + 3 * c_j].im = -(ar / R_im);
          } else {
            R[b_i + 3 * c_j].re = ai / R_im;
            R[b_i + 3 * c_j].im = -(ar / R_im);
          }
        } else {
          brm = fabs(R_re);
          bim = fabs(R_im);
          if (brm > bim) {
            s = R_im / R_re;
            d = R_re + s * R_im;
            R[b_i + 3 * c_j].re = (ar + s * ai) / d;
            R[b_i + 3 * c_j].im = (ai - s * ar) / d;
          } else if (bim == brm) {
            if (R_re > 0.0) {
              sgnbr = 0.5;
            } else {
              sgnbr = -0.5;
            }

            if (R_im > 0.0) {
              sgnbi = 0.5;
            } else {
              sgnbi = -0.5;
            }

            R[b_i + 3 * c_j].re = (ar * sgnbr + ai * sgnbi) / brm;
            R[b_i + 3 * c_j].im = (ai * sgnbr - ar * sgnbi) / brm;
          } else {
            s = R_re / R_im;
            d = R_im + s * R_re;
            R[b_i + 3 * c_j].re = (s * ar + ai) / d;
            R[b_i + 3 * c_j].im = (s * ai - ar) / d;
          }
        }
      }
    }
  }

  for (i119 = 0; i119 < 3; i119++) {
    for (i120 = 0; i120 < 3; i120++) {
      b_Q[i119 + 3 * i120].re = 0.0;
      b_Q[i119 + 3 * i120].im = 0.0;
      for (i121 = 0; i121 < 3; i121++) {
        b_Q[i119 + 3 * i120].re += Q[i119 + 3 * i121].re * R[i121 + 3 * i120].re
          - Q[i119 + 3 * i121].im * R[i121 + 3 * i120].im;
        b_Q[i119 + 3 * i120].im += Q[i119 + 3 * i121].re * R[i121 + 3 * i120].im
          + Q[i119 + 3 * i121].im * R[i121 + 3 * i120].re;
      }
    }

    for (i122 = 0; i122 < 3; i122++) {
      X[i119 + 3 * i122].re = 0.0;
      X[i119 + 3 * i122].im = 0.0;
      for (i123 = 0; i123 < 3; i123++) {
        Q_re = Q[i122 + 3 * i123].re;
        Q_im = -Q[i122 + 3 * i123].im;
        X[i119 + 3 * i122].re += b_Q[i119 + 3 * i123].re * Q_re - b_Q[i119 + 3 *
          i123].im * Q_im;
        X[i119 + 3 * i122].im += b_Q[i119 + 3 * i123].re * Q_im + b_Q[i119 + 3 *
          i123].im * Q_re;
      }
    }
  }

  for (i124 = 0; i124 < 9; i124++) {
    x[i124] = X[i124].im;
  }

  y = 0.0;
  d_j = 0;
  exitg2 = false;
  while ((!exitg2) && (d_j < 3)) {
    b_s = 0.0;
    for (c_i = 0; c_i < 3; c_i++) {
      b_s += fabs(x[c_i + 3 * d_j]);
    }

    if (rtIsNaN(b_s)) {
      y = rtNaN;
      exitg2 = true;
    } else {
      if (b_s > y) {
        y = b_s;
      }

      d_j++;
    }
  }

  b_y = 0.0;
  e_j = 0;
  exitg1 = false;
  while ((!exitg1) && (e_j < 3)) {
    c_s = 0.0;
    for (d_i = 0; d_i < 3; d_i++) {
      c_s += rt_hypotd_snf(X[d_i + 3 * e_j].re, X[d_i + 3 * e_j].im);
    }

    if (rtIsNaN(c_s)) {
      b_y = rtNaN;
      exitg1 = true;
    } else {
      if (c_s > b_y) {
        b_y = c_s;
      }

      e_j++;
    }
  }

  if (y <= 6.6613381477509392E-15 * b_y) {
    for (f_j = 0; f_j < 3; f_j++) {
      for (e_i = 0; e_i < 3; e_i++) {
        X[e_i + 3 * f_j].im = 0.0;
      }
    }
  }
}

/*
 * Arguments    : double *a
 *                double *b
 *                double *c
 *                double *d
 *                double *rt1r
 *                double *rt1i
 *                double *rt2r
 *                double *rt2i
 *                double *cs
 *                double *sn
 * Return Type  : void
 */
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn)
{
  double temp;
  double p;
  double y;
  double b_y;
  double bcmax;
  double c_y;
  double d_y;
  double e_y;
  int b_b;
  int b_c;
  double bcmis;
  double f_y;
  double scale;
  double z;
  double b_a;
  double b_p;
  double tau;
  double sigma;
  int b_sigma;
  double aa;
  double bb;
  double cc;
  double dd;
  double sab;
  double sac;
  double c_a;
  double cs1;
  double sn1;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    temp = *d;
    *d = *a;
    *a = temp;
    *b = -*c;
    *c = 0.0;
  } else if ((*a - *d == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
    *cs = 1.0;
    *sn = 0.0;
  } else {
    temp = *a - *d;
    p = 0.5 * temp;
    y = fabs(*b);
    b_y = fabs(*c);
    if ((y >= b_y) || rtIsNaN(b_y)) {
      bcmax = y;
    } else {
      bcmax = b_y;
    }

    c_y = fabs(*b);
    d_y = fabs(*c);
    if ((c_y <= d_y) || rtIsNaN(d_y)) {
      e_y = c_y;
    } else {
      e_y = d_y;
    }

    if (!(*b < 0.0)) {
      b_b = 1;
    } else {
      b_b = -1;
    }

    if (!(*c < 0.0)) {
      b_c = 1;
    } else {
      b_c = -1;
    }

    bcmis = e_y * (double)b_b * (double)b_c;
    f_y = fabs(p);
    if ((f_y >= bcmax) || rtIsNaN(bcmax)) {
      scale = f_y;
    } else {
      scale = bcmax;
    }

    z = p / scale * p + bcmax / scale * bcmis;
    if (z >= 8.8817841970012523E-16) {
      b_a = sqrt(scale) * sqrt(z);
      if (!(p < 0.0)) {
        b_p = b_a;
      } else {
        b_p = -b_a;
      }

      z = p + b_p;
      *a = *d + z;
      *d -= bcmax / z * bcmis;
      tau = rt_hypotd_snf(*c, z);
      *cs = z / tau;
      *sn = *c / tau;
      *b -= *c;
      *c = 0.0;
    } else {
      sigma = *b + *c;
      tau = rt_hypotd_snf(sigma, temp);
      *cs = sqrt(0.5 * (1.0 + fabs(sigma) / tau));
      if (!(sigma < 0.0)) {
        b_sigma = 1;
      } else {
        b_sigma = -1;
      }

      *sn = -(p / (tau * *cs)) * (double)b_sigma;
      aa = *a * *cs + *b * *sn;
      bb = -*a * *sn + *b * *cs;
      cc = *c * *cs + *d * *sn;
      dd = -*c * *sn + *d * *cs;
      *b = bb * *cs + dd * *sn;
      *c = -aa * *sn + cc * *cs;
      temp = 0.5 * ((aa * *cs + cc * *sn) + (-bb * *sn + dd * *cs));
      *a = temp;
      *d = temp;
      if (*c != 0.0) {
        if (*b != 0.0) {
          if ((*b < 0.0) == (*c < 0.0)) {
            sab = sqrt(fabs(*b));
            sac = sqrt(fabs(*c));
            c_a = sab * sac;
            if (!(*c < 0.0)) {
              p = c_a;
            } else {
              p = -c_a;
            }

            tau = 1.0 / sqrt(fabs(*b + *c));
            *a = temp + p;
            *d = temp - p;
            *b -= *c;
            *c = 0.0;
            cs1 = sab * tau;
            sn1 = sac * tau;
            temp = *cs * cs1 - *sn * sn1;
            *sn = *cs * sn1 + *sn * cs1;
            *cs = temp;
          }
        } else {
          *b = -*c;
          *c = 0.0;
          temp = *cs;
          *cs = -*sn;
          *sn = temp;
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
}

/*
 * Arguments    : double a[9]
 *                double tau[2]
 * Return Type  : void
 */
static void xgehrd(double a[9], double tau[2])
{
  double work[3];
  int i;
  int b_i;
  int im1n;
  int in;
  int c;
  double alpha1;
  double d12;
  double xnorm;
  double beta1;
  int knt;
  int i138;
  int k;
  int i139;
  int b_k;
  int c_k;
  int i140;
  int d_k;
  int b_c;
  int lastv;
  int c_i;
  int lastc;
  boolean_T exitg4;
  int rowleft;
  int ia;
  int exitg3;
  int b_lastc;
  int iy;
  int ix;
  int i141;
  int iac;
  int b_iy;
  int i142;
  int b_ia;
  int jA;
  int jy;
  int j;
  double temp;
  int b_ix;
  int i143;
  int ijA;
  int c_c;
  int d_c;
  int b_lastv;
  int d_i;
  int c_lastc;
  boolean_T exitg2;
  int coltop;
  int c_ia;
  int exitg1;
  int d_lastc;
  int c_iy;
  int d_iy;
  int i144;
  int b_iac;
  int c_ix;
  double e_c;
  int i145;
  int d_ia;
  for (i = 0; i < 3; i++) {
    work[i] = 0.0;
  }

  for (b_i = 0; b_i < 2; b_i++) {
    im1n = b_i * 3 + 2;
    in = (b_i + 1) * 3;
    c = b_i * 3 + 2;
    alpha1 = a[(b_i + 3 * b_i) + 1];
    d12 = 0.0;
    xnorm = 0.0;
    if (1 - b_i < 1) {
    } else {
      xnorm = fabs(a[c]);
    }

    if (xnorm != 0.0) {
      beta1 = rt_hypotd_snf(a[(b_i + 3 * b_i) + 1], xnorm);
      if (a[(b_i + 3 * b_i) + 1] >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          i138 = c - b_i;
          for (k = c; k + 1 <= i138 + 1; k++) {
            a[k] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(beta1) >= 1.0020841800044864E-292));

        xnorm = 0.0;
        if (1 - b_i < 1) {
        } else {
          xnorm = fabs(a[c]);
        }

        beta1 = rt_hypotd_snf(alpha1, xnorm);
        if (alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        d12 = (beta1 - alpha1) / beta1;
        alpha1 = 1.0 / (alpha1 - beta1);
        i139 = c - b_i;
        for (b_k = c; b_k + 1 <= i139 + 1; b_k++) {
          a[b_k] *= alpha1;
        }

        for (c_k = 1; c_k <= knt; c_k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1 = beta1;
      } else {
        d12 = (beta1 - a[(b_i + 3 * b_i) + 1]) / beta1;
        alpha1 = 1.0 / (a[(b_i + 3 * b_i) + 1] - beta1);
        i140 = c - b_i;
        for (d_k = c; d_k + 1 <= i140 + 1; d_k++) {
          a[d_k] *= alpha1;
        }

        alpha1 = beta1;
      }
    }

    tau[b_i] = d12;
    a[(b_i + 3 * b_i) + 1] = 1.0;
    b_c = b_i + im1n;
    if (tau[b_i] != 0.0) {
      lastv = 2 - b_i;
      c_i = b_c - b_i;
      while ((lastv > 0) && (a[c_i] == 0.0)) {
        lastv--;
        c_i--;
      }

      lastc = 3;
      exitg4 = false;
      while ((!exitg4) && (lastc > 0)) {
        rowleft = in + lastc;
        ia = rowleft;
        do {
          exitg3 = 0;
          if (ia <= rowleft + (lastv - 1) * 3) {
            if (a[ia - 1] != 0.0) {
              exitg3 = 1;
            } else {
              ia += 3;
            }
          } else {
            lastc--;
            exitg3 = 2;
          }
        } while (exitg3 == 0);

        if (exitg3 == 1) {
          exitg4 = true;
        }
      }

      b_lastc = lastc;
    } else {
      lastv = 0;
      b_lastc = 0;
    }

    if (lastv > 0) {
      if (b_lastc == 0) {
      } else {
        for (iy = 1; iy <= b_lastc; iy++) {
          work[iy - 1] = 0.0;
        }

        ix = b_c;
        i141 = (in + 3 * (lastv - 1)) + 1;
        for (iac = in + 1; iac <= i141; iac += 3) {
          b_iy = 0;
          i142 = (iac + b_lastc) - 1;
          for (b_ia = iac; b_ia <= i142; b_ia++) {
            work[b_iy] += a[b_ia - 1] * a[ix - 1];
            b_iy++;
          }

          ix++;
        }
      }

      if (-tau[b_i] == 0.0) {
      } else {
        jA = in;
        jy = b_c - 1;
        for (j = 1; j <= lastv; j++) {
          if (a[jy] != 0.0) {
            temp = a[jy] * -tau[b_i];
            b_ix = 0;
            i143 = b_lastc + jA;
            for (ijA = jA; ijA + 1 <= i143; ijA++) {
              a[ijA] += work[b_ix] * temp;
              b_ix++;
            }
          }

          jy++;
          jA += 3;
        }
      }
    }

    c_c = b_i + im1n;
    d_c = (b_i + in) + 2;
    if (tau[b_i] != 0.0) {
      b_lastv = 2 - b_i;
      d_i = c_c - b_i;
      while ((b_lastv > 0) && (a[d_i] == 0.0)) {
        b_lastv--;
        d_i--;
      }

      c_lastc = 2 - b_i;
      exitg2 = false;
      while ((!exitg2) && (c_lastc > 0)) {
        coltop = d_c + (c_lastc - 1) * 3;
        c_ia = coltop;
        do {
          exitg1 = 0;
          if (c_ia <= (coltop + b_lastv) - 1) {
            if (a[c_ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              c_ia++;
            }
          } else {
            c_lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }

      d_lastc = c_lastc;
    } else {
      b_lastv = 0;
      d_lastc = 0;
    }

    if (b_lastv > 0) {
      if (d_lastc == 0) {
      } else {
        for (c_iy = 1; c_iy <= d_lastc; c_iy++) {
          work[c_iy - 1] = 0.0;
        }

        d_iy = 0;
        i144 = d_c + 3 * (d_lastc - 1);
        for (b_iac = d_c; b_iac <= i144; b_iac += 3) {
          c_ix = c_c;
          e_c = 0.0;
          i145 = (b_iac + b_lastv) - 1;
          for (d_ia = b_iac; d_ia <= i145; d_ia++) {
            e_c += a[d_ia - 1] * a[c_ix - 1];
            c_ix++;
          }

          work[d_iy] += e_c;
          d_iy++;
        }
      }

      xgerc(b_lastv, d_lastc, -tau[b_i], c_c, work, a, d_c);
    }

    a[(b_i + 3 * b_i) + 1] = alpha1;
  }
}

/*
 * Arguments    : int m
 *                int n
 *                double alpha1
 *                int ix0
 *                const double y[3]
 *                double A[9]
 *                int ia0
 * Return Type  : void
 */
static void xgerc(int m, int n, double alpha1, int ix0, const double y[3],
                  double A[9], int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i146;
  int ijA;
  if (alpha1 == 0.0) {
  } else {
    jA = ia0 - 1;
    jy = 0;
    for (j = 1; j <= n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i146 = m + jA;
        for (ijA = jA; ijA + 1 <= i146; ijA++) {
          A[ijA] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

/*
 * Arguments    : double h[9]
 *                double z[9]
 * Return Type  : int
 */
static int xhseqr(double h[9], double z[9])
{
  int b_info;
  double v[3];
  int i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  double tst;
  boolean_T guard1 = false;
  double htmp1;
  double htmp2;
  double ab;
  double ba;
  double aa;
  double bb;
  double s;
  double varargin_2;
  double d13;
  double h11;
  double h12;
  double h21;
  double h22;
  double rt1r;
  double rt1i;
  double rt2r;
  double rt2i;
  double tr;
  double det;
  double rtdisc;
  double h21s;
  int b_k;
  int c;
  int nr;
  int hoffset;
  int j;
  double alpha;
  double t1;
  double v2;
  double t2;
  double v3;
  double t3;
  int b_j;
  double sum1;
  int i147;
  int c_j;
  int d_j;
  int e_j;
  int f_j;
  int g_j;
  double d14;
  double d15;
  double d16;
  double sn;
  double cs;
  double unusedU3;
  double unusedU2;
  double unusedU1;
  double unusedU0;
  b_info = -1;
  h[2] = 0.0;
  i = 2;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 31)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && ((k + 1 > 1) && (!(fabs(h[k + 3 * (k - 1)]) <=
                3.0062525400134592E-292)))) {
        tst = fabs(h[(k + 3 * (k - 1)) - 1]) + fabs(h[k + 3 * k]);
        if (tst == 0.0) {
          if (k - 1 >= 1) {
            tst = fabs(h[(k + 3 * (k - 2)) - 1]);
          }

          if (k + 2 <= 3) {
            tst += fabs(h[(k + 3 * k) + 1]);
          }
        }

        guard1 = false;
        if (fabs(h[k + 3 * (k - 1)]) <= 2.2204460492503131E-16 * tst) {
          htmp1 = fabs(h[k + 3 * (k - 1)]);
          htmp2 = fabs(h[(k + 3 * k) - 1]);
          if (htmp1 > htmp2) {
            ab = htmp1;
            ba = htmp2;
          } else {
            ab = htmp2;
            ba = htmp1;
          }

          htmp1 = fabs(h[k + 3 * k]);
          htmp2 = fabs(h[(k + 3 * (k - 1)) - 1] - h[k + 3 * k]);
          if (htmp1 > htmp2) {
            aa = htmp1;
            bb = htmp2;
          } else {
            aa = htmp2;
            bb = htmp1;
          }

          s = aa + ab;
          varargin_2 = 2.2204460492503131E-16 * (bb * (aa / s));
          if ((3.0062525400134592E-292 >= varargin_2) || rtIsNaN(varargin_2)) {
            d13 = 3.0062525400134592E-292;
          } else {
            d13 = varargin_2;
          }

          if (ba * (ab / s) <= d13) {
            exitg3 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          k--;
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 3 * (k - 1)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          s = fabs(h[1]) + fabs(h[5]);
          h11 = 0.75 * s + h[0];
          h12 = -0.4375 * s;
          h21 = s;
          h22 = h11;
        } else if (its == 20) {
          s = fabs(h[i + 3 * (i - 1)]) + fabs(h[i - 1]);
          h11 = 0.75 * s + h[i + 3 * i];
          h12 = -0.4375 * s;
          h21 = s;
          h22 = h11;
        } else {
          h11 = h[(i + 3 * (i - 1)) - 1];
          h21 = h[i + 3 * (i - 1)];
          h12 = h[(i + 3 * i) - 1];
          h22 = h[i + 3 * i];
        }

        s = ((fabs(h11) + fabs(h12)) + fabs(h21)) + fabs(h22);
        if (s == 0.0) {
          rt1r = 0.0;
          rt1i = 0.0;
          rt2r = 0.0;
          rt2i = 0.0;
        } else {
          h11 /= s;
          h21 /= s;
          h12 /= s;
          h22 /= s;
          tr = (h11 + h22) / 2.0;
          det = (h11 - tr) * (h22 - tr) - h12 * h21;
          rtdisc = sqrt(fabs(det));
          if (det >= 0.0) {
            rt1r = tr * s;
            rt2r = rt1r;
            rt1i = rtdisc * s;
            rt2i = -rt1i;
          } else {
            rt1r = tr + rtdisc;
            rt2r = tr - rtdisc;
            if (fabs(rt1r - h22) <= fabs(rt2r - h22)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }

            rt1i = 0.0;
            rt2i = 0.0;
          }
        }

        if (i - 1 >= 1) {
          s = (fabs(h[0] - rt2r) + fabs(rt2i)) + fabs(h[1]);
          h21s = h[1] / s;
          v[0] = (h21s * h[3] + (h[0] - rt1r) * ((h[0] - rt2r) / s)) - rt1i *
            (rt2i / s);
          v[1] = h21s * (((h[0] + h[4]) - rt1r) - rt2r);
          v[2] = h21s * h[5];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
        }

        for (b_k = i - 1; b_k <= i; b_k++) {
          c = (i - b_k) + 2;
          if (3 <= c) {
            nr = 3;
          } else {
            nr = c;
          }

          if (b_k > i - 1) {
            hoffset = b_k + 3 * (b_k - 2);
            for (j = 1; j <= nr; j++) {
              v[j - 1] = h[(j + hoffset) - 2];
            }
          }

          alpha = v[0];
          t1 = xzlarfg(nr, &alpha, v);
          v[0] = alpha;
          if (b_k > i - 1) {
            h[b_k - 1] = alpha;
            h[b_k] = 0.0;
            if (b_k < i) {
              h[2] = 0.0;
            }
          }

          v2 = v[1];
          t2 = t1 * v[1];
          if (nr == 3) {
            v3 = v[2];
            t3 = t1 * v[2];
            for (b_j = b_k - 1; b_j + 1 < 4; b_j++) {
              sum1 = (h[(b_k + 3 * b_j) - 1] + v2 * h[b_k + 3 * b_j]) + v3 * h
                [(b_k + 3 * b_j) + 1];
              h[(b_k + 3 * b_j) - 1] -= sum1 * t1;
              h[b_k + 3 * b_j] -= sum1 * t2;
              h[2 + 3 * b_j] -= sum1 * t3;
            }

            if (b_k + 3 <= i + 1) {
              i147 = b_k;
            } else {
              i147 = i - 2;
            }

            for (c_j = 0; c_j + 1 <= i147 + 3; c_j++) {
              sum1 = (h[c_j + 3 * (b_k - 1)] + v2 * h[c_j + 3 * b_k]) + v3 *
                h[c_j + 3 * (b_k + 1)];
              h[c_j + 3 * (b_k - 1)] -= sum1 * t1;
              h[c_j + 3 * b_k] -= sum1 * t2;
              h[6 + c_j] -= sum1 * t3;
            }

            for (d_j = 0; d_j < 3; d_j++) {
              sum1 = (z[d_j + 3 * (b_k - 1)] + v2 * z[d_j + 3 * b_k]) + v3 *
                z[d_j + 3 * (b_k + 1)];
              z[d_j + 3 * (b_k - 1)] -= sum1 * t1;
              z[d_j + 3 * b_k] -= sum1 * t2;
              z[6 + d_j] -= sum1 * t3;
            }
          } else {
            if (nr == 2) {
              for (e_j = b_k - 1; e_j + 1 < 4; e_j++) {
                sum1 = h[(b_k + 3 * e_j) - 1] + v2 * h[b_k + 3 * e_j];
                h[(b_k + 3 * e_j) - 1] -= sum1 * t1;
                h[b_k + 3 * e_j] -= sum1 * t2;
              }

              for (f_j = 0; f_j + 1 <= i + 1; f_j++) {
                sum1 = h[f_j + 3 * (b_k - 1)] + v2 * h[f_j + 3 * b_k];
                h[f_j + 3 * (b_k - 1)] -= sum1 * t1;
                h[f_j + 3 * b_k] -= sum1 * t2;
              }

              for (g_j = 0; g_j < 3; g_j++) {
                sum1 = z[g_j + 3 * (b_k - 1)] + v2 * z[g_j + 3 * b_k];
                z[g_j + 3 * (b_k - 1)] -= sum1 * t1;
                z[g_j + 3 * b_k] -= sum1 * t2;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      b_info = i;
      exitg1 = true;
    } else {
      if ((L == i + 1) || (!(L == i))) {
      } else {
        d14 = h[(i + 3 * i) - 1];
        d15 = h[i + 3 * (i - 1)];
        d16 = h[i + 3 * i];
        xdlanv2(&h[(i + 3 * (i - 1)) - 1], &d14, &d15, &d16, &unusedU0,
                &unusedU1, &unusedU2, &unusedU3, &cs, &sn);
        h[(i + 3 * i) - 1] = d14;
        h[i + 3 * (i - 1)] = d15;
        h[i + 3 * i] = d16;
        if (3 > i + 1) {
          xrot(2 - i, h, i + (i + 1) * 3, (i + (i + 1) * 3) + 1, cs, sn);
        }

        b_xrot(i - 1, h, 1 + (i - 1) * 3, 1 + i * 3, cs, sn);
        c_xrot(z, 1 + (i - 1) * 3, 1 + i * 3, cs, sn);
      }

      i = L - 2;
    }
  }

  return b_info + 1;
}

/*
 * Arguments    : int n
 *                const double x[3]
 * Return Type  : double
 */
static double xnrm2(int n, const double x[3])
{
  double y;
  double scale;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x[1]);
  } else {
    scale = 2.2250738585072014E-308;
    for (k = 0; k < 2; k++) {
      absxk = fabs(x[k + 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

/*
 * Arguments    : int n
 *                double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
static void xrot(int n, double x[9], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = c * x[ix] + s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = temp;
    iy += 3;
    ix += 3;
  }
}

/*
 * Arguments    : int n
 *                double *alpha1
 *                double x[3]
 * Return Type  : double
 */
static double xzlarfg(int n, double *alpha1, double x[3])
{
  double tau;
  double xnorm;
  double beta1;
  int knt;
  int k;
  int b_k;
  int c_k;
  int d_k;
  tau = 0.0;
  if (n <= 0) {
  } else {
    xnorm = xnrm2(n - 1, x);
    if (xnorm != 0.0) {
      beta1 = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          for (k = 1; k + 1 <= n; k++) {
            x[k] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(beta1) >= 1.0020841800044864E-292));

        beta1 = rt_hypotd_snf(*alpha1, xnrm2(n - 1, x));
        if (*alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        tau = (beta1 - *alpha1) / beta1;
        *alpha1 = 1.0 / (*alpha1 - beta1);
        for (b_k = 1; b_k + 1 <= n; b_k++) {
          x[b_k] *= *alpha1;
        }

        for (c_k = 1; c_k <= knt; c_k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        *alpha1 = beta1;
      } else {
        tau = (beta1 - *alpha1) / beta1;
        *alpha1 = 1.0 / (*alpha1 - beta1);
        for (d_k = 1; d_k + 1 <= n; d_k++) {
          x[d_k] *= *alpha1;
        }

        *alpha1 = beta1;
      }
    }
  }

  return tau;
}

/*
 * Arguments    : double i
 *                const struct0_T *Param
 *                double z_ukf
 *                double state_ukf[6]
 *                const double w_ukf[3]
 *                double est_ukf[6]
 *                const double P[9]
 *                const double Q[9]
 *                double R
 *                double dt
 *                double ndim
 *                double zdim
 *                struct1_T *ukfout
 * Return Type  : void
 */
void UKF_bhm3(double i, const struct0_T *Param, double z_ukf, double state_ukf[6],
              const double w_ukf[3], double est_ukf[6], const double P[9], const
              double Q[9], double R, double dt, double ndim, double zdim,
              struct1_T *ukfout)
{
  double A;
  double y;
  double SOC;
  double Cb;
  double Ccp;
  double Vcp;
  double Vp;
  double ib;
  double icp;
  double b_state_ukf[3];
  int i0;
  emxArray_real_T *W;
  int i1;
  int loop_ub;
  int i2;
  double d0;
  double dv0[9];
  int i3;
  creal_T dcv0[9];
  double SQRTM_P[9];
  int i4;
  emxArray_real_T *Chi;
  int i5;
  int b_loop_ub;
  int i6;
  emxArray_int32_T *r0;
  int c_loop_ub;
  int i7;
  int i8;
  double b_est_ukf[3];
  int i9;
  int unnamed_idx_0;
  int i10;
  int k;
  int d_loop_ub;
  int i11;
  int i12;
  double c_est_ukf[3];
  int i13;
  int b_unnamed_idx_0;
  int i14;
  int e_loop_ub;
  int i15;
  int i16;
  int i17;
  double d_est_ukf[3];
  int i18;
  int c_unnamed_idx_0;
  int i19;
  emxArray_real_T *ChiF;
  int i20;
  int f_loop_ub;
  int i21;
  double d1;
  int b_k;
  int g_loop_ub;
  int i22;
  int i23;
  double b_Chi[3];
  int d_unnamed_idx_0;
  int i24;
  emxArray_real_T *ChiF_mean;
  int i25;
  int h_loop_ub;
  int i26;
  double d2;
  int c_k;
  double b_W;
  int i27;
  int i_loop_ub;
  int i28;
  emxArray_real_T *PX;
  int i29;
  int j_loop_ub;
  int i30;
  double d3;
  int d_k;
  emxArray_real_T *c_W;
  emxArray_real_T *b_ChiF;
  double d_W;
  int k_loop_ub;
  int l_loop_ub;
  int i31;
  int i32;
  int i33;
  int i34;
  int i35;
  int m_loop_ub;
  int i36;
  int n_loop_ub;
  int i37;
  double d4;
  double a[3];
  emxArray_real_T *ZF;
  int i38;
  int o_loop_ub;
  int i39;
  int i40;
  unsigned int unnamed_idx_1;
  int i41;
  int n;
  int i42;
  int p_loop_ub;
  int i43;
  int cr;
  int ic;
  int br;
  int b_cr;
  int ar;
  int b_ib;
  int ia;
  int b_ic;
  emxArray_real_T *ZF_mean;
  int i44;
  int q_loop_ub;
  int i45;
  double d5;
  int e_k;
  double e_W;
  int i46;
  int r_loop_ub;
  int i47;
  emxArray_real_T *PZ;
  int i48;
  int s_loop_ub;
  int i49;
  double d6;
  int f_k;
  emxArray_real_T *f_W;
  emxArray_real_T *b_ZF;
  double g_W;
  double c_ZF;
  double d_ZF;
  int i50;
  int t_loop_ub;
  int i51;
  int i52;
  int u_loop_ub;
  int i53;
  int i54;
  int v_loop_ub;
  int i55;
  int w_loop_ub;
  int i56;
  double d7;
  int i57;
  int b_PZ;
  int c_PZ;
  int x_loop_ub;
  int i58;
  emxArray_real_T *PXZ;
  int i59;
  int y_loop_ub;
  int i60;
  double d8;
  int g_k;
  emxArray_real_T *h_W;
  emxArray_real_T *e_ZF;
  double i_W;
  int ab_loop_ub;
  double f_ZF;
  int i61;
  int i62;
  int i63;
  int bb_loop_ub;
  int i64;
  int i65;
  int cb_loop_ub;
  int i66;
  int db_loop_ub;
  int i67;
  double d9;
  emxArray_real_T *K;
  int i68;
  int eb_loop_ub;
  int i69;
  int fb_loop_ub;
  int i70;
  int gb_loop_ub;
  int i71;
  int h_k;
  unsigned int e_unnamed_idx_0;
  unsigned int b_unnamed_idx_1;
  int i72;
  int m;
  int i73;
  int hb_loop_ub;
  int i74;
  int ib_loop_ub;
  int i75;
  int c;
  int c_cr;
  int i76;
  int c_ic;
  int b_br;
  int d_cr;
  int b_ar;
  int i77;
  int c_ib;
  int b_ia;
  int i78;
  int d_ic;
  int i79;
  int jb_loop_ub;
  int i80;
  emxArray_real_T *C;
  int i81;
  int kb_loop_ub;
  int i82;
  int lb_loop_ub;
  int i83;
  int i_k;
  unsigned int K_idx_0;
  int i84;
  int b_m;
  int b_C;
  int i85;
  int i86;
  int e_cr;
  int e_ic;
  int c_br;
  int f_cr;
  int c_ar;
  int i87;
  int d_ib;
  int c_ia;
  int f_ic;
  emxArray_real_T *b_y;
  int i88;
  int mb_loop_ub;
  int i89;
  int nb_loop_ub;
  int i90;
  int ob_loop_ub;
  int i91;
  int j_k;
  int i92;
  int c_m;
  int i93;
  int pb_loop_ub;
  int i94;
  int qb_loop_ub;
  int i95;
  int b_c;
  int g_cr;
  int i96;
  int g_ic;
  int d_br;
  int h_cr;
  int d_ar;
  int i97;
  int e_ib;
  int d_ia;
  int i98;
  int h_ic;
  int i99;
  int rb_loop_ub;
  int i100;
  int sb_loop_ub;
  int i101;
  emxArray_real_T *c_C;
  int i102;
  int tb_loop_ub;
  int i103;
  int ub_loop_ub;
  int i104;
  int vb_loop_ub;
  int i105;
  int k_k;
  int i106;
  int d_m;
  int i107;
  int wb_loop_ub;
  int i108;
  int xb_loop_ub;
  int i109;
  int c_c;
  int i_cr;
  int i110;
  int i_ic;
  int e_br;
  int j_cr;
  int e_ar;
  int i111;
  int f_ib;
  int e_ia;
  int i112;
  int j_ic;
  emxArray_real_T *b_ChiF_mean;
  int i113;
  int yb_loop_ub;
  int i114;
  int i115;
  int i116;
  int i117;

  /* function [P,state_ukf,est_ukf, SOC] = UKF_bhm3(i,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim) */
  /* function [t_predict,P,state_ukf,est_ukf] = UKF_bhm(z_ukf,tau_model,state_ukf,w_ukf,est_ukf,P,Q,R,dt,xa0,ndim) */
  /*  Function expanded to N dimensional state vector */
  /*  Array indices must begin at 1 and not 0 */
  /*  Sigma points are denoted by the "Chi" variables */
  /* %SOC = ones(1,N);  */
  /* %state_meas = ones(2,1); */
  /* heat transfer coeff */
  /* ambient temp degC */
  /* SOC(1) = 100; */
  /* %x_total(:,1) = x; */
  /* Tb = state_ukf(4,1); %Removed Tb from state and meas vectors */
  A = Param->qmax - state_ukf[0];
  y = A / Param->Cmax;
  SOC = 1.0 - A / Param->Cmax;
  Cb = ((Param->Ccb0 + Param->Ccb1 * (1.0 - y)) + Param->Ccb2 * (SOC * SOC)) +
    Param->Ccb3 * rt_powd_snf(SOC, 3.0);
  Ccp = Param->Ccp0 + Param->Ccp1 * exp(Param->Ccp2 * (1.0 - y));

  /* qhat = qmax - (qmax*0.99 - Cmax)*(1-SOC); */
  /* Vb = qhat/Cb; %Diverges */
  /* Works with this! Only diverges at very end. */
  Vcp = state_ukf[1] / Ccp;
  Vp = (state_ukf[0] / Cb - Vcp) - state_ukf[2] / Param->Cs;
  ib = Vp / Param->Rp + i;
  icp = ib - Vcp / (Param->Rcp0 + Param->Rcp1 * exp(Param->Rcp2 * (1.0 - y)));

  /* Tbdot = 1/Jt*(hcp*Vcp^2/Rcp + hcs*Vcs^2/Rs - ha*(Tb-Ta)); */
  /* {   */
  /* x = [ qb + qbdot*dt; */
  /*         qcp + qcpdot*dt; */
  /*         qcsnew; */
  /*         Tb + Tbdot*dt];     */
  /* } */
  /* voltage and tau */
  b_state_ukf[0] = state_ukf[0] + -ib * dt;
  b_state_ukf[1] = state_ukf[1] + icp * dt;
  b_state_ukf[2] = ib * Param->Rs * Param->Cs;
  for (i0 = 0; i0 < 3; i0++) {
    state_ukf[3 + i0] = b_state_ukf[i0] + w_ukf[i0];
  }

  emxInit_real_T(&W, 1);

  /* % Selection of Sigma Points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* W0 = -0.8;  % -1 < W0 < 1 */
  /* W0 = 0.99;              %Weights for tau are different than for voltage */
  /* W1 = (1-W0)/(2*2); */
  i1 = W->size[0];
  W->size[0] = (int)(2.0 * ndim + 1.0);
  emxEnsureCapacity((emxArray__common *)W, i1, (int)sizeof(double));
  loop_ub = (int)(2.0 * ndim + 1.0);
  for (i2 = 0; i2 < loop_ub; i2++) {
    W->data[i2] = 0.0;
  }

  W->data[0] = 0.99;

  /* Sigma point array begins at 1 and not 0 */
  d0 = 1.0 - W->data[0];
  for (i3 = 0; i3 < 9; i3++) {
    dv0[i3] = 2.0 * P[i3] / d0;
  }

  sqrtm(dv0, dcv0);
  for (i4 = 0; i4 < 9; i4++) {
    SQRTM_P[i4] = dcv0[i4].re;
  }

  emxInit_real_T1(&Chi, 2);

  /* Square root of matrix is converted to real type */
  i5 = Chi->size[0] * Chi->size[1];
  Chi->size[0] = (int)ndim;
  Chi->size[1] = (int)(2.0 * ndim + 1.0);
  emxEnsureCapacity((emxArray__common *)Chi, i5, (int)sizeof(double));
  b_loop_ub = (int)ndim * (int)(2.0 * ndim + 1.0);
  for (i6 = 0; i6 < b_loop_ub; i6++) {
    Chi->data[i6] = 0.0;
  }

  emxInit_int32_T(&r0, 1);
  c_loop_ub = (int)ndim;
  i7 = r0->size[0];
  r0->size[0] = (int)ndim;
  emxEnsureCapacity((emxArray__common *)r0, i7, (int)sizeof(int));
  for (i8 = 0; i8 < c_loop_ub; i8++) {
    r0->data[i8] = i8;
  }

  for (i9 = 0; i9 < 3; i9++) {
    b_est_ukf[i9] = est_ukf[i9];
  }

  unnamed_idx_0 = r0->size[0];
  for (i10 = 0; i10 < unnamed_idx_0; i10++) {
    Chi->data[r0->data[i10]] = b_est_ukf[i10];
  }

  /*  xa is old state estimate feed back from algorithm */
  for (k = 0; k < (int)((ndim + 1.0) + -1.0); k++) {
    /* Generate sigma points with +sqrt. Increase indices by 1 since arrays don't start from 0 */
    d_loop_ub = Chi->size[0];
    i11 = r0->size[0];
    r0->size[0] = d_loop_ub;
    emxEnsureCapacity((emxArray__common *)r0, i11, (int)sizeof(int));
    for (i12 = 0; i12 < d_loop_ub; i12++) {
      r0->data[i12] = i12;
    }

    for (i13 = 0; i13 < 3; i13++) {
      c_est_ukf[i13] = est_ukf[i13] + SQRTM_P[i13 + 3 * ((int)((2.0 + (double)k)
        - 1.0) - 1)];
    }

    b_unnamed_idx_0 = r0->size[0];
    for (i14 = 0; i14 < b_unnamed_idx_0; i14++) {
      Chi->data[r0->data[i14] + Chi->size[0] * (k + 1)] = c_est_ukf[i14];
    }

    W->data[k + 1] = (1.0 - W->data[0]) / 4.0;
    e_loop_ub = Chi->size[0];
    i15 = r0->size[0];
    r0->size[0] = e_loop_ub;
    emxEnsureCapacity((emxArray__common *)r0, i15, (int)sizeof(int));
    for (i16 = 0; i16 < e_loop_ub; i16++) {
      r0->data[i16] = i16;
    }

    i17 = (int)(ndim + (2.0 + (double)k)) - 1;
    for (i18 = 0; i18 < 3; i18++) {
      d_est_ukf[i18] = est_ukf[i18] - SQRTM_P[i18 + 3 * ((int)((2.0 + (double)k)
        - 1.0) - 1)];
    }

    c_unnamed_idx_0 = r0->size[0];
    for (i19 = 0; i19 < c_unnamed_idx_0; i19++) {
      Chi->data[r0->data[i19] + Chi->size[0] * i17] = d_est_ukf[i19];
    }

    W->data[(int)(ndim + (2.0 + (double)k)) - 1] = (1.0 - W->data[0]) / 4.0;
  }

  emxInit_real_T1(&ChiF, 2);

  /* % Model Forecast Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Propagate Sigma Points through non-linear transform */
  i20 = ChiF->size[0] * ChiF->size[1];
  ChiF->size[0] = (int)ndim;
  ChiF->size[1] = (int)(2.0 * ndim + 1.0);
  emxEnsureCapacity((emxArray__common *)ChiF, i20, (int)sizeof(double));
  f_loop_ub = (int)ndim * (int)(2.0 * ndim + 1.0);
  for (i21 = 0; i21 < f_loop_ub; i21++) {
    ChiF->data[i21] = 0.0;
  }

  d1 = 2.0 * ndim + 1.0;
  for (b_k = 0; b_k < (int)d1; b_k++) {
    /* Forecast sigma points */
    /* %Consider finding qcsdot */
    /* Symbol for sigma pt is chi */
    g_loop_ub = ChiF->size[0];
    i22 = r0->size[0];
    r0->size[0] = g_loop_ub;
    emxEnsureCapacity((emxArray__common *)r0, i22, (int)sizeof(int));
    for (i23 = 0; i23 < g_loop_ub; i23++) {
      r0->data[i23] = i23;
    }

    b_Chi[0] = Chi->data[Chi->size[0] * b_k] + -ib * dt;
    b_Chi[1] = Chi->data[1 + Chi->size[0] * b_k] + icp * dt;
    b_Chi[2] = (Vp / Param->Rp + i) * Param->Rs * Param->Cs;
    d_unnamed_idx_0 = r0->size[0];
    for (i24 = 0; i24 < d_unnamed_idx_0; i24++) {
      ChiF->data[r0->data[i24] + ChiF->size[0] * b_k] = b_Chi[i24];
    }

    /* %% Not sure about this!!! */
    /* { */
    /* ChiF(:,k) = [ xa(1) + qbdot*dt;    %Symbol for sigma pt is chi */
    /*         xa(2) + qcpdot*dt; */
    /*         Chi_qcsnew;            %%% Not sure about this!!! */
    /*         xa(4) + Tbdot*dt]; */
    /* } */
    /* a = (1/(1+dt/tau_model))*Chi(1,k); */
    /* tau_meas = (1/(-1+Chi(1,1)/a))*dt;  %There is an issue forecasting tau.  To work better, need 2 points in history */
    /* ChiF(:,k) = [a; tau_model]+w_ukf; */
    /* ChiF(:,k) = [a; tau_meas]; */
  }

  emxFree_int32_T(&r0);
  emxInit_real_T(&ChiF_mean, 1);

  /* % Compute the mean and covariance of forecast */
  i25 = ChiF_mean->size[0];
  ChiF_mean->size[0] = (int)ndim;
  emxEnsureCapacity((emxArray__common *)ChiF_mean, i25, (int)sizeof(double));
  h_loop_ub = (int)ndim;
  for (i26 = 0; i26 < h_loop_ub; i26++) {
    ChiF_mean->data[i26] = 0.0;
  }

  d2 = 2.0 * ndim + 1.0;
  for (c_k = 0; c_k < (int)d2; c_k++) {
    /* Calculate mean of sigma points */
    b_W = W->data[c_k];
    i27 = ChiF_mean->size[0];
    emxEnsureCapacity((emxArray__common *)ChiF_mean, i27, (int)sizeof(double));
    i_loop_ub = ChiF_mean->size[0];
    for (i28 = 0; i28 < i_loop_ub; i28++) {
      ChiF_mean->data[i28] += b_W * ChiF->data[i28 + ChiF->size[0] * c_k];
    }
  }

  emxInit_real_T1(&PX, 2);
  i29 = PX->size[0] * PX->size[1];
  PX->size[0] = (int)ndim;
  PX->size[1] = (int)ndim;
  emxEnsureCapacity((emxArray__common *)PX, i29, (int)sizeof(double));
  j_loop_ub = (int)ndim * (int)ndim;
  for (i30 = 0; i30 < j_loop_ub; i30++) {
    PX->data[i30] = 0.0;
  }

  /* PX = zeros(ndim,ndim);       */
  d3 = 2.0 * ndim + 1.0;
  d_k = 0;
  emxInit_real_T(&c_W, 1);
  emxInit_real_T1(&b_ChiF, 2);
  while (d_k <= (int)d3 - 1) {
    /* Calculate covariance matrix of sigma pts */
    d_W = W->data[d_k];
    k_loop_ub = ChiF->size[0];
    l_loop_ub = ChiF->size[0];
    i31 = c_W->size[0];
    c_W->size[0] = k_loop_ub;
    emxEnsureCapacity((emxArray__common *)c_W, i31, (int)sizeof(double));
    for (i32 = 0; i32 < k_loop_ub; i32++) {
      c_W->data[i32] = d_W * (ChiF->data[i32 + ChiF->size[0] * d_k] -
        ChiF_mean->data[i32]);
    }

    i33 = b_ChiF->size[0] * b_ChiF->size[1];
    b_ChiF->size[0] = 1;
    b_ChiF->size[1] = l_loop_ub;
    emxEnsureCapacity((emxArray__common *)b_ChiF, i33, (int)sizeof(double));
    for (i34 = 0; i34 < l_loop_ub; i34++) {
      b_ChiF->data[b_ChiF->size[0] * i34] = ChiF->data[i34 + ChiF->size[0] * d_k]
        - ChiF_mean->data[i34];
    }

    i35 = PX->size[0] * PX->size[1];
    PX->size[0] = c_W->size[0];
    PX->size[1] = b_ChiF->size[1];
    emxEnsureCapacity((emxArray__common *)PX, i35, (int)sizeof(double));
    m_loop_ub = c_W->size[0];
    for (i36 = 0; i36 < m_loop_ub; i36++) {
      n_loop_ub = b_ChiF->size[1];
      for (i37 = 0; i37 < n_loop_ub; i37++) {
        d4 = c_W->data[i36] * b_ChiF->data[b_ChiF->size[0] * i37];
        PX->data[i36 + PX->size[0] * i37] += d4;
      }
    }

    d_k++;
  }

  emxFree_real_T(&b_ChiF);
  emxFree_real_T(&c_W);

  /* % Propagate the sigma points through the observation model */
  /* ZF = ChiF;                  %Expected measurements based on forecasted sigma points */
  a[0] = 1.0 / Cb;
  a[1] = -1.0 / Ccp;
  a[2] = -1.0 / Param->Cs;
  emxInit_real_T1(&ZF, 2);
  if (ChiF->size[0] == 1) {
    i38 = ZF->size[0] * ZF->size[1];
    ZF->size[0] = 1;
    ZF->size[1] = ChiF->size[1];
    emxEnsureCapacity((emxArray__common *)ZF, i38, (int)sizeof(double));
    o_loop_ub = ChiF->size[1];
    for (i39 = 0; i39 < o_loop_ub; i39++) {
      ZF->data[ZF->size[0] * i39] = 0.0;
      for (i40 = 0; i40 < 3; i40++) {
        ZF->data[ZF->size[0] * i39] += a[i40] * ChiF->data[i40 + ChiF->size[0] *
          i39];
      }
    }
  } else {
    unnamed_idx_1 = (unsigned int)ChiF->size[1];
    i41 = ZF->size[0] * ZF->size[1];
    ZF->size[0] = 1;
    ZF->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)ZF, i41, (int)sizeof(double));
    n = ChiF->size[1] - 1;
    i42 = ZF->size[0] * ZF->size[1];
    ZF->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)ZF, i42, (int)sizeof(double));
    p_loop_ub = ZF->size[1];
    for (i43 = 0; i43 < p_loop_ub; i43++) {
      ZF->data[ZF->size[0] * i43] = 0.0;
    }

    for (cr = 0; cr <= n; cr++) {
      for (ic = cr; ic + 1 <= cr + 1; ic++) {
        ZF->data[ic] = 0.0;
      }
    }

    br = 0;
    for (b_cr = 0; b_cr <= n; b_cr++) {
      ar = 0;
      for (b_ib = br; b_ib + 1 <= br + 3; b_ib++) {
        if (ChiF->data[b_ib] != 0.0) {
          ia = ar;
          for (b_ic = b_cr; b_ic + 1 <= b_cr + 1; b_ic++) {
            ia++;
            ZF->data[b_ic] += ChiF->data[b_ib] * a[ia - 1];
          }
        }

        ar++;
      }

      br += 3;
    }
  }

  emxInit_real_T(&ZF_mean, 1);

  /* Expected measurements based on forecasted sigma points */
  /* % Compute mean and covariance of transformed observations */
  i44 = ZF_mean->size[0];
  ZF_mean->size[0] = (int)zdim;
  emxEnsureCapacity((emxArray__common *)ZF_mean, i44, (int)sizeof(double));
  q_loop_ub = (int)zdim;
  for (i45 = 0; i45 < q_loop_ub; i45++) {
    ZF_mean->data[i45] = 0.0;
  }

  /* ZF_mean = zeros(ndim,1); */
  d5 = 2.0 * ndim + 1.0;
  for (e_k = 0; e_k < (int)d5; e_k++) {
    /* Calculate mean of expected measurements */
    e_W = W->data[e_k] * ZF->data[ZF->size[0] * e_k];
    i46 = ZF_mean->size[0];
    emxEnsureCapacity((emxArray__common *)ZF_mean, i46, (int)sizeof(double));
    r_loop_ub = ZF_mean->size[0];
    for (i47 = 0; i47 < r_loop_ub; i47++) {
      ZF_mean->data[i47] += e_W;
    }
  }

  emxInit_real_T1(&PZ, 2);
  i48 = PZ->size[0] * PZ->size[1];
  PZ->size[0] = (int)zdim;
  PZ->size[1] = (int)zdim;
  emxEnsureCapacity((emxArray__common *)PZ, i48, (int)sizeof(double));
  s_loop_ub = (int)zdim * (int)zdim;
  for (i49 = 0; i49 < s_loop_ub; i49++) {
    PZ->data[i49] = 0.0;
  }

  /* PZ = zeros(ndim,ndim); */
  d6 = 2.0 * ndim + 1.0;
  f_k = 0;
  emxInit_real_T(&f_W, 1);
  emxInit_real_T1(&b_ZF, 2);
  while (f_k <= (int)d6 - 1) {
    /* Covar of expected measurements */
    g_W = W->data[f_k];
    c_ZF = ZF->data[ZF->size[0] * f_k];
    d_ZF = ZF->data[ZF->size[0] * f_k];
    i50 = f_W->size[0];
    f_W->size[0] = ZF_mean->size[0];
    emxEnsureCapacity((emxArray__common *)f_W, i50, (int)sizeof(double));
    t_loop_ub = ZF_mean->size[0];
    for (i51 = 0; i51 < t_loop_ub; i51++) {
      f_W->data[i51] = g_W * (c_ZF - ZF_mean->data[i51]);
    }

    i52 = b_ZF->size[0] * b_ZF->size[1];
    b_ZF->size[0] = 1;
    b_ZF->size[1] = ZF_mean->size[0];
    emxEnsureCapacity((emxArray__common *)b_ZF, i52, (int)sizeof(double));
    u_loop_ub = ZF_mean->size[0];
    for (i53 = 0; i53 < u_loop_ub; i53++) {
      b_ZF->data[b_ZF->size[0] * i53] = d_ZF - ZF_mean->data[i53];
    }

    i54 = PZ->size[0] * PZ->size[1];
    PZ->size[0] = f_W->size[0];
    PZ->size[1] = b_ZF->size[1];
    emxEnsureCapacity((emxArray__common *)PZ, i54, (int)sizeof(double));
    v_loop_ub = f_W->size[0];
    for (i55 = 0; i55 < v_loop_ub; i55++) {
      w_loop_ub = b_ZF->size[1];
      for (i56 = 0; i56 < w_loop_ub; i56++) {
        d7 = f_W->data[i55] * b_ZF->data[b_ZF->size[0] * i56];
        PZ->data[i55 + PZ->size[0] * i56] += d7;
      }
    }

    f_k++;
  }

  emxFree_real_T(&b_ZF);
  emxFree_real_T(&f_W);
  i57 = PZ->size[0] * PZ->size[1];
  emxEnsureCapacity((emxArray__common *)PZ, i57, (int)sizeof(double));
  b_PZ = PZ->size[0];
  c_PZ = PZ->size[1];
  x_loop_ub = b_PZ * c_PZ;
  for (i58 = 0; i58 < x_loop_ub; i58++) {
    PZ->data[i58] += R;
  }

  emxInit_real_T1(&PXZ, 2);

  /* % Compute the cross covariance between XF and Zf */
  i59 = PXZ->size[0] * PXZ->size[1];
  PXZ->size[0] = (int)ndim;
  PXZ->size[1] = (int)zdim;
  emxEnsureCapacity((emxArray__common *)PXZ, i59, (int)sizeof(double));
  y_loop_ub = (int)ndim * (int)zdim;
  for (i60 = 0; i60 < y_loop_ub; i60++) {
    PXZ->data[i60] = 0.0;
  }

  /* PXZ = zeros(ndim,ndim); */
  d8 = 2.0 * ndim + 1.0;
  g_k = 0;
  emxInit_real_T(&h_W, 1);
  emxInit_real_T1(&e_ZF, 2);
  while (g_k <= (int)d8 - 1) {
    /* Cross covariance of sigma points and expected measurements */
    i_W = W->data[g_k];
    ab_loop_ub = ChiF->size[0];
    f_ZF = ZF->data[ZF->size[0] * g_k];
    i61 = h_W->size[0];
    h_W->size[0] = ab_loop_ub;
    emxEnsureCapacity((emxArray__common *)h_W, i61, (int)sizeof(double));
    for (i62 = 0; i62 < ab_loop_ub; i62++) {
      h_W->data[i62] = i_W * (ChiF->data[i62 + ChiF->size[0] * g_k] -
        ChiF_mean->data[i62]);
    }

    i63 = e_ZF->size[0] * e_ZF->size[1];
    e_ZF->size[0] = 1;
    e_ZF->size[1] = ZF_mean->size[0];
    emxEnsureCapacity((emxArray__common *)e_ZF, i63, (int)sizeof(double));
    bb_loop_ub = ZF_mean->size[0];
    for (i64 = 0; i64 < bb_loop_ub; i64++) {
      e_ZF->data[e_ZF->size[0] * i64] = f_ZF - ZF_mean->data[i64];
    }

    i65 = PXZ->size[0] * PXZ->size[1];
    PXZ->size[0] = h_W->size[0];
    PXZ->size[1] = e_ZF->size[1];
    emxEnsureCapacity((emxArray__common *)PXZ, i65, (int)sizeof(double));
    cb_loop_ub = h_W->size[0];
    for (i66 = 0; i66 < cb_loop_ub; i66++) {
      db_loop_ub = e_ZF->size[1];
      for (i67 = 0; i67 < db_loop_ub; i67++) {
        d9 = h_W->data[i66] * e_ZF->data[e_ZF->size[0] * i67];
        PXZ->data[i66 + PXZ->size[0] * i67] += d9;
      }
    }

    g_k++;
  }

  emxFree_real_T(&e_ZF);
  emxFree_real_T(&h_W);
  emxFree_real_T(&ZF);
  emxFree_real_T(&ChiF);
  emxFree_real_T(&W);

  /* % Data Assimilation Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* K = PXZ*pinv(PZ); */
  mpower(PZ, Chi);
  emxInit_real_T1(&K, 2);
  if ((PXZ->size[1] == 1) || (Chi->size[0] == 1)) {
    i68 = K->size[0] * K->size[1];
    K->size[0] = PXZ->size[0];
    K->size[1] = Chi->size[1];
    emxEnsureCapacity((emxArray__common *)K, i68, (int)sizeof(double));
    eb_loop_ub = PXZ->size[0];
    for (i69 = 0; i69 < eb_loop_ub; i69++) {
      fb_loop_ub = Chi->size[1];
      for (i70 = 0; i70 < fb_loop_ub; i70++) {
        K->data[i69 + K->size[0] * i70] = 0.0;
        gb_loop_ub = PXZ->size[1];
        for (i71 = 0; i71 < gb_loop_ub; i71++) {
          K->data[i69 + K->size[0] * i70] += PXZ->data[i69 + PXZ->size[0] * i71]
            * Chi->data[i71 + Chi->size[0] * i70];
        }
      }
    }
  } else {
    h_k = PXZ->size[1];
    e_unnamed_idx_0 = (unsigned int)PXZ->size[0];
    b_unnamed_idx_1 = (unsigned int)Chi->size[1];
    i72 = K->size[0] * K->size[1];
    K->size[0] = (int)e_unnamed_idx_0;
    K->size[1] = (int)b_unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)K, i72, (int)sizeof(double));
    m = PXZ->size[0];
    i73 = K->size[0] * K->size[1];
    emxEnsureCapacity((emxArray__common *)K, i73, (int)sizeof(double));
    hb_loop_ub = K->size[1];
    for (i74 = 0; i74 < hb_loop_ub; i74++) {
      ib_loop_ub = K->size[0];
      for (i75 = 0; i75 < ib_loop_ub; i75++) {
        K->data[i75 + K->size[0] * i74] = 0.0;
      }
    }

    if ((PXZ->size[0] == 0) || (Chi->size[1] == 0)) {
    } else {
      c = PXZ->size[0] * (Chi->size[1] - 1);
      c_cr = 0;
      while ((m > 0) && (c_cr <= c)) {
        i76 = c_cr + m;
        for (c_ic = c_cr; c_ic + 1 <= i76; c_ic++) {
          K->data[c_ic] = 0.0;
        }

        c_cr += m;
      }

      b_br = 0;
      d_cr = 0;
      while ((m > 0) && (d_cr <= c)) {
        b_ar = 0;
        i77 = b_br + h_k;
        for (c_ib = b_br; c_ib + 1 <= i77; c_ib++) {
          if (Chi->data[c_ib] != 0.0) {
            b_ia = b_ar;
            i78 = d_cr + m;
            for (d_ic = d_cr; d_ic + 1 <= i78; d_ic++) {
              b_ia++;
              K->data[d_ic] += Chi->data[c_ib] * PXZ->data[b_ia - 1];
            }
          }

          b_ar += m;
        }

        b_br += h_k;
        d_cr += m;
      }
    }
  }

  emxFree_real_T(&PXZ);
  i79 = ZF_mean->size[0];
  emxEnsureCapacity((emxArray__common *)ZF_mean, i79, (int)sizeof(double));
  jb_loop_ub = ZF_mean->size[0];
  for (i80 = 0; i80 < jb_loop_ub; i80++) {
    ZF_mean->data[i80] = z_ukf - ZF_mean->data[i80];
  }

  emxInit_real_T(&C, 1);
  if ((K->size[1] == 1) || (ZF_mean->size[0] == 1)) {
    i81 = C->size[0];
    C->size[0] = K->size[0];
    emxEnsureCapacity((emxArray__common *)C, i81, (int)sizeof(double));
    kb_loop_ub = K->size[0];
    for (i82 = 0; i82 < kb_loop_ub; i82++) {
      C->data[i82] = 0.0;
      lb_loop_ub = K->size[1];
      for (i83 = 0; i83 < lb_loop_ub; i83++) {
        C->data[i82] += K->data[i82 + K->size[0] * i83] * ZF_mean->data[i83];
      }
    }
  } else {
    i_k = K->size[1];
    K_idx_0 = (unsigned int)K->size[0];
    i84 = C->size[0];
    C->size[0] = (int)K_idx_0;
    emxEnsureCapacity((emxArray__common *)C, i84, (int)sizeof(double));
    b_m = K->size[0];
    b_C = C->size[0];
    i85 = C->size[0];
    C->size[0] = b_C;
    emxEnsureCapacity((emxArray__common *)C, i85, (int)sizeof(double));
    for (i86 = 0; i86 < b_C; i86++) {
      C->data[i86] = 0.0;
    }

    if (K->size[0] == 0) {
    } else {
      e_cr = 0;
      while ((b_m > 0) && (e_cr <= 0)) {
        for (e_ic = 1; e_ic <= b_m; e_ic++) {
          C->data[e_ic - 1] = 0.0;
        }

        e_cr = b_m;
      }

      c_br = 0;
      f_cr = 0;
      while ((b_m > 0) && (f_cr <= 0)) {
        c_ar = 0;
        i87 = c_br + i_k;
        for (d_ib = c_br; d_ib + 1 <= i87; d_ib++) {
          if (ZF_mean->data[d_ib] != 0.0) {
            c_ia = c_ar;
            for (f_ic = 0; f_ic + 1 <= b_m; f_ic++) {
              c_ia++;
              C->data[f_ic] += ZF_mean->data[d_ib] * K->data[c_ia - 1];
            }
          }

          c_ar += b_m;
        }

        c_br += i_k;
        f_cr = b_m;
      }
    }
  }

  emxFree_real_T(&ZF_mean);
  emxInit_real_T1(&b_y, 2);
  if ((K->size[1] == 1) || (PZ->size[0] == 1)) {
    i88 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = K->size[0];
    b_y->size[1] = PZ->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i88, (int)sizeof(double));
    mb_loop_ub = K->size[0];
    for (i89 = 0; i89 < mb_loop_ub; i89++) {
      nb_loop_ub = PZ->size[1];
      for (i90 = 0; i90 < nb_loop_ub; i90++) {
        b_y->data[i89 + b_y->size[0] * i90] = 0.0;
        ob_loop_ub = K->size[1];
        for (i91 = 0; i91 < ob_loop_ub; i91++) {
          b_y->data[i89 + b_y->size[0] * i90] += K->data[i89 + K->size[0] * i91]
            * PZ->data[i91 + PZ->size[0] * i90];
        }
      }
    }
  } else {
    j_k = K->size[1];
    e_unnamed_idx_0 = (unsigned int)K->size[0];
    b_unnamed_idx_1 = (unsigned int)PZ->size[1];
    i92 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = (int)e_unnamed_idx_0;
    b_y->size[1] = (int)b_unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)b_y, i92, (int)sizeof(double));
    c_m = K->size[0];
    i93 = b_y->size[0] * b_y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i93, (int)sizeof(double));
    pb_loop_ub = b_y->size[1];
    for (i94 = 0; i94 < pb_loop_ub; i94++) {
      qb_loop_ub = b_y->size[0];
      for (i95 = 0; i95 < qb_loop_ub; i95++) {
        b_y->data[i95 + b_y->size[0] * i94] = 0.0;
      }
    }

    if ((K->size[0] == 0) || (PZ->size[1] == 0)) {
    } else {
      b_c = K->size[0] * (PZ->size[1] - 1);
      g_cr = 0;
      while ((c_m > 0) && (g_cr <= b_c)) {
        i96 = g_cr + c_m;
        for (g_ic = g_cr; g_ic + 1 <= i96; g_ic++) {
          b_y->data[g_ic] = 0.0;
        }

        g_cr += c_m;
      }

      d_br = 0;
      h_cr = 0;
      while ((c_m > 0) && (h_cr <= b_c)) {
        d_ar = 0;
        i97 = d_br + j_k;
        for (e_ib = d_br; e_ib + 1 <= i97; e_ib++) {
          if (PZ->data[e_ib] != 0.0) {
            d_ia = d_ar;
            i98 = h_cr + c_m;
            for (h_ic = h_cr; h_ic + 1 <= i98; h_ic++) {
              d_ia++;
              b_y->data[h_ic] += PZ->data[e_ib] * K->data[d_ia - 1];
            }
          }

          d_ar += c_m;
        }

        d_br += j_k;
        h_cr += c_m;
      }
    }
  }

  emxFree_real_T(&PZ);
  i99 = Chi->size[0] * Chi->size[1];
  Chi->size[0] = K->size[1];
  Chi->size[1] = K->size[0];
  emxEnsureCapacity((emxArray__common *)Chi, i99, (int)sizeof(double));
  rb_loop_ub = K->size[0];
  for (i100 = 0; i100 < rb_loop_ub; i100++) {
    sb_loop_ub = K->size[1];
    for (i101 = 0; i101 < sb_loop_ub; i101++) {
      Chi->data[i101 + Chi->size[0] * i100] = K->data[i100 + K->size[0] * i101];
    }
  }

  emxFree_real_T(&K);
  emxInit_real_T1(&c_C, 2);
  if ((b_y->size[1] == 1) || (Chi->size[0] == 1)) {
    i102 = c_C->size[0] * c_C->size[1];
    c_C->size[0] = b_y->size[0];
    c_C->size[1] = Chi->size[1];
    emxEnsureCapacity((emxArray__common *)c_C, i102, (int)sizeof(double));
    tb_loop_ub = b_y->size[0];
    for (i103 = 0; i103 < tb_loop_ub; i103++) {
      ub_loop_ub = Chi->size[1];
      for (i104 = 0; i104 < ub_loop_ub; i104++) {
        c_C->data[i103 + c_C->size[0] * i104] = 0.0;
        vb_loop_ub = b_y->size[1];
        for (i105 = 0; i105 < vb_loop_ub; i105++) {
          c_C->data[i103 + c_C->size[0] * i104] += b_y->data[i103 + b_y->size[0]
            * i105] * Chi->data[i105 + Chi->size[0] * i104];
        }
      }
    }
  } else {
    k_k = b_y->size[1];
    e_unnamed_idx_0 = (unsigned int)b_y->size[0];
    b_unnamed_idx_1 = (unsigned int)Chi->size[1];
    i106 = c_C->size[0] * c_C->size[1];
    c_C->size[0] = (int)e_unnamed_idx_0;
    c_C->size[1] = (int)b_unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)c_C, i106, (int)sizeof(double));
    d_m = b_y->size[0];
    i107 = c_C->size[0] * c_C->size[1];
    emxEnsureCapacity((emxArray__common *)c_C, i107, (int)sizeof(double));
    wb_loop_ub = c_C->size[1];
    for (i108 = 0; i108 < wb_loop_ub; i108++) {
      xb_loop_ub = c_C->size[0];
      for (i109 = 0; i109 < xb_loop_ub; i109++) {
        c_C->data[i109 + c_C->size[0] * i108] = 0.0;
      }
    }

    if ((b_y->size[0] == 0) || (Chi->size[1] == 0)) {
    } else {
      c_c = b_y->size[0] * (Chi->size[1] - 1);
      i_cr = 0;
      while ((d_m > 0) && (i_cr <= c_c)) {
        i110 = i_cr + d_m;
        for (i_ic = i_cr; i_ic + 1 <= i110; i_ic++) {
          c_C->data[i_ic] = 0.0;
        }

        i_cr += d_m;
      }

      e_br = 0;
      j_cr = 0;
      while ((d_m > 0) && (j_cr <= c_c)) {
        e_ar = 0;
        i111 = e_br + k_k;
        for (f_ib = e_br; f_ib + 1 <= i111; f_ib++) {
          if (Chi->data[f_ib] != 0.0) {
            e_ia = e_ar;
            i112 = j_cr + d_m;
            for (j_ic = j_cr; j_ic + 1 <= i112; j_ic++) {
              e_ia++;
              c_C->data[j_ic] += Chi->data[f_ib] * b_y->data[e_ia - 1];
            }
          }

          e_ar += d_m;
        }

        e_br += k_k;
        j_cr += d_m;
      }
    }
  }

  emxFree_real_T(&b_y);
  emxFree_real_T(&Chi);
  emxInit_real_T(&b_ChiF_mean, 1);
  i113 = b_ChiF_mean->size[0];
  b_ChiF_mean->size[0] = ChiF_mean->size[0];
  emxEnsureCapacity((emxArray__common *)b_ChiF_mean, i113, (int)sizeof(double));
  yb_loop_ub = ChiF_mean->size[0];
  for (i114 = 0; i114 < yb_loop_ub; i114++) {
    b_ChiF_mean->data[i114] = ChiF_mean->data[i114] + C->data[i114];
  }

  emxFree_real_T(&C);
  emxFree_real_T(&ChiF_mean);
  for (i115 = 0; i115 < 3; i115++) {
    est_ukf[3 + i115] = b_ChiF_mean->data[i115];
  }

  emxFree_real_T(&b_ChiF_mean);

  /* t_predict = -xa(2,1)*log(0.5/xa0); */
  for (i116 = 0; i116 < 9; i116++) {
    ukfout->P[i116] = (PX->data[i116] + Q[i116]) - c_C->data[i116];
  }

  emxFree_real_T(&c_C);
  emxFree_real_T(&PX);
  for (i117 = 0; i117 < 6; i117++) {
    ukfout->state_ukf[i117] = state_ukf[i117];
    ukfout->est_ukf[i117] = est_ukf[i117];
  }

  ukfout->SOC = 1.0 - y;
 // ukfout.SOC = 0.1234;   //Test successful JCS 4.16.16
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void UKF_bhm3_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void UKF_bhm3_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for UKF_bhm3.c
 *
 * [EOF]
 */
