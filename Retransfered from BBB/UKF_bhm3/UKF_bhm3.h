/*
 * File: UKF_bhm3.h
 *
 * MATLAB Coder version            : 3.0
 * C/C++ source code generated on  : 05-Apr-2016 17:59:50
 */

#ifndef __UKF_BHM3_H__
#define __UKF_BHM3_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "UKF_bhm3_types.h"

/* Function Declarations */
void UKF_bhm3(double i, const struct0_T *Param, double z_ukf, double state_ukf[6],
              const double w_ukf[3], double est_ukf[6], const double P[9], const
              double Q[9], double R, double dt, double ndim, double zdim,
              struct1_T *ukfout);
void UKF_bhm3_initialize(void);
void UKF_bhm3_terminate(void);

#endif

/*
 * File trailer for UKF_bhm3.h
 *
 * [EOF]
 */
