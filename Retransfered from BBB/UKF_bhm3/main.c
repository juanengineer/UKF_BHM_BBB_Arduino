/* This is a copy of the original main.c file in the Examples folder //JCS 4.5.16 
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
#include <math.h>   //JCS 4.5.16 
#include <stdio.h>  //JCS 4.5.16 	

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
  result->qmax = 2.88*pow(10,4);
  result->Cmax = 2.85*pow(10,4);
  result->Ccb0 = 19.4;
  result->Ccb1 = 1576;
  result->Ccb2 = 41.7;
  result->Ccb3 = -203;
  result->Rs = 2.77*pow(10,-2);
  result->Cs = 89.3;
  result->Rcp0 = 1.6*pow(10,-3);
  result->Rcp1 = 8.45;
  result->Rcp2 = -61.9;
  result->Ccp0 = 2689;
  result->Ccp1 = -2285;
  result->Ccp2 = -0.73;
  result->Rp = 1*pow(10,5);
  result->Jt = 800;
  result->hcp = 19;
  result->hcs = 1;
  result->ha = 0.5;
  result->Ta = 25;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_UKF_bhm3(void)
{
/*
  
  double i;
  struct0_T Param;
  double dv1[6];
  double dv2[3];
  double dv3[6];
  double dv4[9];
  double dv5[9];
  struct1_T ukfout;
*/   
//Original JCS 4.5.16
  int incrementvariable;
  double i_testarray[5]={60,60,60,60,60};
  double z_testarray[5]={0, -320.5640, 7.3615, 11.6589, 14.2643};
  double i = 1;
  struct0_T Param;
  double z_ukf = 1;
  double state_ukf[6]={1,1,1,1,1,1};
  double w_ukf[3]={0,0,0};
  double est_ukf[6]={1,1,1,1,1,1};
  double P[9]={1,1,1,1,1,1,1,1,1};
  double Q[9]={1,1,1,1,1,1,1,1,1};
  double R=1;
  double dt=1;
  double ndim=3;
  double zdim=1;
  struct1_T ukfout;      //Modified from argument declarations from UKF_bhm3.c
  double *Ptest = ukfout.P;
  double *state_ukf_test = ukfout.state_ukf;
  double *est_ukf_test = ukfout.est_ukf;


  Param.qmax = 2.88*pow(10,4);
  Param.Cmax = 2.85*pow(10,4);
  Param.Ccb0 = 19.4;
  Param.Ccb1 = 1576;
  Param.Ccb2 = 41.7;
  Param.Ccb3 = -203;
  Param.Rs = 2.77*pow(10,-2);
  Param.Cs = 89.3;
  Param.Rcp0 = 1.6*pow(10,-3);
  Param.Rcp1 = 8.45;
  Param.Rcp2 = -61.9;
  Param.Ccp0 = 2689;
  Param.Ccp1 = -2285;
  Param.Ccp2 = -0.73;
  Param.Rp = 1*pow(10,5);
  Param.Jt = 800;
  Param.hcp = 19;
  Param.hcs = 1;
  Param.ha = 0.5;
  Param.Ta = 25;

  ukfout.P[0] = 1;
  ukfout.P[1] = 1;
  ukfout.P[2] = 1;
  ukfout.P[3] = 1;
  ukfout.P[4] = 1;
  ukfout.P[5] = 1;
  ukfout.P[6] = 1;
  ukfout.P[7] = 1;
  ukfout.P[8] = 1;
  ukfout.P[9] = 1;

  ukfout.state_ukf[0] = 1;
  ukfout.state_ukf[1] = 1;
  ukfout.state_ukf[2] = 1;
  ukfout.state_ukf[3] = 1;
  ukfout.state_ukf[4] = 1;
  ukfout.state_ukf[5] = 1;

  ukfout.est_ukf[0] = 1;
  ukfout.est_ukf[1] = 1;
  ukfout.est_ukf[2] = 1;
  ukfout.est_ukf[3] = 1;
  ukfout.est_ukf[4] = 1;
  ukfout.est_ukf[5] = 1;
 
  ukfout.SOC = 0.6969;
  
  /* Initialize function 'UKF_bhm3' input arguments. */
  i = argInit_real_T();

  /* Initialize function input argument 'Param'. */
  //argInit_struct0_T(&Param);

  /* Initialize function input argument 'state_ukf'. */
  /* Initialize function input argument 'w_ukf'. */
  /* Initialize function input argument 'est_ukf'. */
  /* Initialize function input argument 'P'. */
  /* Initialize function input argument 'Q'. */
  /* Call the entry-point 'UKF_bhm3'. */
//  argInit_3x2_real_T(dv1);
//  argInit_3x1_real_T(dv2);
//  argInit_3x2_real_T(dv3);
// argInit_3x3_real_T(dv4);
//  argInit_3x3_real_T(dv5);
/*
printf("The qmax value is %f\n",Param.qmax);		//JCS 4.5.16
printf("The Cmax value is %f\n",Param.Cmax);		//JCS 4.5.16
printf("The Ccb0 value is %f\n",Param.Ccb0);		//JCS 4.5.16
printf("The Ccb1 value is %f\n",Param.Ccb1);		//JCS 4.5.16
printf("The Ccb2 value is %f\n",Param.Ccb2);		//JCS 4.5.16
printf("The Ccb3 value is %f\n",Param.Ccb3);		//JCS 4.5.16
printf("The Rs value is %f\n",Param.Rs);		//JCS 4.5.16
printf("The Cs value is %f\n",Param.Cs);		//JCS 4.5.16
printf("The Rcp0 value is %f\n",Param.Rcp0);		//JCS 4.5.16
printf("The Rcp1 value is %f\n",Param.Rcp1);		//JCS 4.5.16
printf("The Rcp2 value is %f\n",Param.Rcp2);		//JCS 4.5.16
printf("The Ccp0 value is %f\n",Param.Ccp0);		//JCS 4.5.16
printf("The Ccp1 value is %f\n",Param.Ccp1);		//JCS 4.5.16
printf("The Ccp2 value is %f\n",Param.Ccp2);		//JCS 4.5.16
printf("The Rp value is %f\n",Param.Rp);		//JCS 4.5.16
printf("The Jt value is %f\n",Param.Jt);		//JCS 4.5.16
printf("The hcp  value is %f\n",Param.hcp );		//JCS 4.5.16
printf("The hcs value is %f\n",Param.hcs);		//JCS 4.5.16
printf("The ha value is %f\n",Param.ha);		//JCS 4.5.16
printf("The Ta value is %f\n",Param.Ta);		//JCS 4.5.16
*/
//test print successful! JCS 4.5.16

/*
UKF_bhm3(i, &Param, argInit_real_T(), dv1, dv2, dv3, dv4, dv5, argInit_real_T(),
           argInit_real_T(), argInit_real_T(), argInit_real_T(), &ukfout);*/ 
// Original function call JCS 4.5.16

for (incrementvariable=0; incrementvariable<5; incrementvariable++) {

i = i_testarray[incrementvariable];
z_ukf = z_testarray[incrementvariable];

UKF_bhm3(i, &Param, z_ukf, state_ukf_test,
              w_ukf, est_ukf_test, Ptest, 
               Q,  R,  dt,  ndim,  zdim,
              &ukfout); 	//JCS 4.5.16 modified from UKF_bhm3.c

Ptest = ukfout.P;     
state_ukf_test = ukfout.state_ukf;
est_ukf_test = ukfout.est_ukf;

//state_ukf = ukfout.state_ukf;
//est_ukf = ukfout.est_ukf;

printf("The SOC value is %f",ukfout.SOC);		//JCS 4.5.16
//printf("The P[0] value is %f ",ukfout.P[4]);		//JCS 4.5.16
printf(". The Ptest[4] value is %f",Ptest[4]);		//JCS 4.5.16
printf(". The state_ukf_test[0] value is %f",state_ukf_test[0]);		//JCS 4.5.16
printf(". The est_ukf_test[0] value is %f\n",est_ukf_test[0]);		//JCS 4.5.16
}

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
