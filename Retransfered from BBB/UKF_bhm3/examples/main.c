/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.0
 * C/C++ source code generated on  : 05-Apr-2016 17:59:50
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "UKF_bhm3.h"
#include "main.h"

/* Function Declarations */
static void argInit_3x1_real_T(double result[3]);
static void argInit_3x2_real_T(double result[6]);
static void argInit_3x3_real_T(double result[9]);
static double argInit_real_T(void);
static void argInit_struct0_T(struct0_T *result);
static void main_UKF_bhm3(void);

/* Function Definitions */

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[6]
 * Return Type  : void
 */
static void argInit_3x2_real_T(double result[6])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 2; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_3x3_real_T(double result[9])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->qmax = argInit_real_T();
  result->Cmax = argInit_real_T();
  result->Ccb0 = argInit_real_T();
  result->Ccb1 = argInit_real_T();
  result->Ccb2 = argInit_real_T();
  result->Ccb3 = argInit_real_T();
  result->Rs = argInit_real_T();
  result->Cs = argInit_real_T();
  result->Rcp0 = argInit_real_T();
  result->Rcp1 = argInit_real_T();
  result->Rcp2 = argInit_real_T();
  result->Ccp0 = argInit_real_T();
  result->Ccp1 = argInit_real_T();
  result->Ccp2 = argInit_real_T();
  result->Rp = argInit_real_T();
  result->Jt = argInit_real_T();
  result->hcp = argInit_real_T();
  result->hcs = argInit_real_T();
  result->ha = argInit_real_T();
  result->Ta = argInit_real_T();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_UKF_bhm3(void)
{
  double i;
  struct0_T Param;
  double dv1[6];
  double dv2[3];
  double dv3[6];
  double dv4[9];
  double dv5[9];
  struct1_T ukfout;

  /* Initialize function 'UKF_bhm3' input arguments. */
  i = argInit_real_T();

  /* Initialize function input argument 'Param'. */
  argInit_struct0_T(&Param);

  /* Initialize function input argument 'state_ukf'. */
  /* Initialize function input argument 'w_ukf'. */
  /* Initialize function input argument 'est_ukf'. */
  /* Initialize function input argument 'P'. */
  /* Initialize function input argument 'Q'. */
  /* Call the entry-point 'UKF_bhm3'. */
  argInit_3x2_real_T(dv1);
  argInit_3x1_real_T(dv2);
  argInit_3x2_real_T(dv3);
  argInit_3x3_real_T(dv4);
  argInit_3x3_real_T(dv5);
  UKF_bhm3(i, &Param, argInit_real_T(), dv1, dv2, dv3, dv4, dv5, argInit_real_T(),
           argInit_real_T(), argInit_real_T(), argInit_real_T(), &ukfout);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  UKF_bhm3_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_UKF_bhm3();

  /* Terminate the application.
     You do not need to do this more than one time. */
  UKF_bhm3_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
