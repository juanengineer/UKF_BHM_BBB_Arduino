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
#include<fcntl.h>   //added from BBB_UKF_4-12.c
#include<unistd.h>  //added from BBB_UKF_4-12.c
#include<termios.h>  //added from BBB_UKF_4-12.c  
#include<string.h>  //added from BBB_UKF_4-12.c
#include<stdlib.h>  //added from BBB_UKF_4-12.c
#include <math.h>   //added from BBB_UKF_4-12.c
#include <time.h>
#include <sys/time.h>
	
#include "rtGetInf.c"
#include "rtGetNaN.c"
#include "rt_nonfinite.c"
#include "UKF_bhm3.c"
//above 4 lines may break in the future, but must investigate. 4-26-2016
//To compile: gcc main.c -lm

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
 			 //Linduino sends data at approx 20 Hz or Tsig= 50000 us
			 //Nyquist sampling at twice the signal frequency or greater
			 // (which is also half the signal period or less)

 //int Linduino_period = 50000; //in microseconds (us)
 int Linduino_period = 754000; //in microseconds (us)
 //int Nyquist_period =  Linduino_period/2.6;   //Division of doubles and then decimal result is truncated.
 int Nyquist_period =  Linduino_period/2.0;   //Seems to work better when the BBB delay period matches the Arduino delay period.
 int Nyquist_frequency = 1.0/(Nyquist_period/1000000.0); 
			
 struct timeval tv;
 double tminus1, t, tzero; 
 gettimeofday(&tv, NULL);
        
 tzero = tv.tv_sec+tv.tv_usec/1000000.0;   //seconds
 tminus1 = 0;   //seconds
 //printf("the time is %f/n",tzero);
 FILE *fp;
 fp = fopen("/media/JUAN_SD/testfolder/testdatafile.csv","w+");
 fprintf(fp,"time       ,dt,       i,         V,           SOC,      t_EOD\n");

//fputs("This is testing for fputs...\n",fp);

 int file, count;  //added from BBB_UKF_4-12.c
 int res; 
 char str4[25];    //added from BBB_UKF_4-12.c
 float val;        //added from BBB_UKF_4-12.c

 //const int N = (Nyquist_frequency*60*20); //Test.  N=72,000 is for 60 Hz for 1200 seconds (20 minutes)
 const int N = 72000; //Test.  N=72,000 is for 60 Hz for 1200 seconds (20 minutes)
 const int ndim=3;
 const int zdim=1;
 int k;
 double i_testarray[]={49,39,29,19,9};
 double z_testarray[]={0, -320.5640, 7.3615, 11.6589, 14.2643};
 double i = 1;  //current and NOT INDEX VARIABLE
 struct0_T Param;
 double z_ukf = 1;
 double w_ukf[N];
 // double w_ukf[]={0.02,0.02,0.02,0.02,0.02};
 double x[ndim*N];
 double x_est[ndim*N];
  
 const double Q[9]={1,0,0,0,1,0,0,0,1}; //process noise covar
 const double R=1;  //Measurement noise covariance
 double dt=1;
  
 struct1_T ukfout;      //Modified from argument declarations from UKF_bhm3.c
 double SOC_array[N]; 
 double *P = ukfout.P;    //state estimate covariance
 double *state_ukf = ukfout.state_ukf;
 double *est_ukf = ukfout.est_ukf;
 double SOC = ukfout.SOC;  //Does not require pointer type because SOC is single value, maybe???
 double qb;
  
 double SOC_EOD = 0.999;
 double voltage_knee = 16.7;
 double newton_error = 1;
 double Cb;
 double Cb_prime;
 double qb_knee;
 double f_of_SOC;
 double f_prime;
 double q_EOD;
 double t_present;
 double t_EOD;
 char buffer[30];
   
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
 
 double fully_charged_battery_voltage = 21.0;
 double Cb_initial = Param.Ccb0 + Param.Ccb1 + Param.Ccb2 + Param.Ccb3;
 double Ccp_initial = Param.Ccp0 + Param.Ccp1*exp(Param.Ccp2);
 double Rcp_initial = Param.Rcp0 + Param.Rcp1*exp(Param.Rcp2);

 double qb_initial = Param.qmax;
 double qcp_initial = fully_charged_battery_voltage/Param.Rp*Ccp_initial*Rcp_initial;
 double qcs_initial = Param.Cs*(qb_initial/Cb_initial-qcp_initial/Ccp_initial - fully_charged_battery_voltage);

 x[0] = qb_initial;  //qb Initial conditions
 x[1] = qcp_initial;  //qcp Initial conditions
 x[2] = qcs_initial;  //qcs Initial conditions
   
 //printf("Initial conditions are %f  %f  %f\n",x[0],x[1],x[2]);

 x_est[0] = 1; 
 x_est[1] = 1; 
 x_est[2] = 1; 

 ukfout.P[0] = 1; ukfout.P[1] = 0; ukfout.P[2] = 0;
 ukfout.P[3] = 0; ukfout.P[4] = 1; ukfout.P[5] = 0;
 ukfout.P[6] = 0; ukfout.P[7] = 0; ukfout.P[8] = 1;
 
 ukfout.state_ukf[0] = 1; ukfout.state_ukf[3] = 1; 
 ukfout.state_ukf[1] = 1; ukfout.state_ukf[4] = 1; 
 ukfout.state_ukf[2] = 1; ukfout.state_ukf[5] = 1; 

 ukfout.est_ukf[0] = 1; ukfout.est_ukf[3] = 1; 
 ukfout.est_ukf[1] = 1; ukfout.est_ukf[4] = 1; 
 ukfout.est_ukf[2] = 1; ukfout.est_ukf[5] = 1;
  
 //Section below has to do with End of Discharge (EOD) time prediction
 //Below is Newton method
 //counter = 0;
 voltage_knee = 16.7;

 SOC_EOD = 0.999; //State of Charge at EOD
 newton_error = 1;

   while (newton_error > 0.01) { //This is fractional error
     Cb = Param.Ccb0 + Param.Ccb1*SOC_EOD+Param.Ccb2*pow(SOC_EOD,2)+Param.Ccb3*pow(SOC_EOD,3);
     Cb_prime = Param.Ccb1 + 2*Param.Ccb2*SOC_EOD + 3*Param.Ccb3*pow(SOC_EOD,2);
     qb_knee = voltage_knee*Cb;
     f_of_SOC = 1-(Param.qmax-qb_knee)/Param.Cmax-SOC_EOD;
     f_prime = -voltage_knee/Param.Cmax*Cb_prime;
     SOC_EOD = SOC_EOD - f_of_SOC/f_prime;
     newton_error = fabs((f_of_SOC/f_prime)/SOC_EOD);
   }//end while

 q_EOD = Param.qmax - Param.Cmax*(1-SOC_EOD); //q_EOD is of qb
 //printf("The q_EOD is %f\n",q_EOD);
 printf("The SOC_EOD is %lf\n",SOC_EOD);
 //printf("Cb_knee is %lf\n",Cb);

//************code inserted below***********************************
 
   if ((file = open("/dev/ttyUSB0", O_RDWR | O_NOCTTY | O_NDELAY))<0){
     perror("USB: Failed to open the file.\n");
     //return -1;
   }//end if

     for (k=1; k<N+1; k++){

	printf("The increment var is %d\n",k);
	i = 0.002;
	z_ukf = 0.2; //For testing value will be changed when read correctly
			// Above 0.1 was commented out to allow previous z-ukf to be used until next update
	
	state_ukf[0]= x[(ndim)*(k-1)];
	state_ukf[1]= x[(ndim)*(k-1)+1];
	state_ukf[2]= x[(ndim)*(k-1)+2];

	est_ukf[0]= x_est[(ndim)*(k-1)];
	est_ukf[1]= x_est[(ndim)*(k-1)+1];
	est_ukf[2]= x_est[(ndim)*(k-1)+2];

        int char_array_length = 22;	      //When battery is connected, 20 chars expected
        struct termios options;               //The termios structure is vital
        tcgetattr(file, &options);            //Sets the parameters associated with file
        //Set up the communications options:
        //9600 baud, 8-bit, enable receiver, no modem control lines
        options.c_cflag = B9600 | CS8 | CREAD | CLOCAL;
	options.c_iflag = IGNPAR | ICRNL;    //ignore partity errors, CR -> newline
   	tcflush(file, TCIFLUSH);             //discard file information not transmitted
   	tcsetattr(file, TCSANOW, &options);  //changes occur immmediately
   	//usleep(55000);                  //give the Arduino a chance to respond in microseconds
	
	usleep(Nyquist_period);		//give the Arduino a chance to respond in microseconds
	
	//unsigned char receivetest[char_array_length];      //declare a buffer for receiving data
	unsigned char receive[char_array_length];      //declare a buffer for receiving data
	

	//if ((count = read(file, (void*)receivetest, char_array_length))<0){ 
   	if ((count = read(file, (void*)receive, char_array_length))<0){   //receive the data
      	  perror("Failed to read from the input\n");
           //return -1;
		//continue; //continue to next iteration in for loop
   	}//end if
	
	//count = 1; //test
	//unsigned char *receive= "i= 60.0000 V=3.14159";
	
	if (count==0) {
		printf("There was no data available to read!\n");
		//continue; //continue to next ieration in for loop
	}
   	
	//if checksum before is equal to checksum after then use string
/*
	unsigned char *ans1, *ans2;
 	unsigned char *receive_pointer = receive;
	int strpos1,strpos2;
	
	ans1 = strstr(receive,"i="); //search for prefix in the receive string
	strpos1 = ans1 - receive_pointer;     //subtract pointers
	
	ans2 = strstr(receive,"V="); 
	strpos2 = ans2 - receive_pointer;     //subtract pointers

	if ((strpos1 < 0)||(strpos1 >= char_array_length)) {
		//continue; //continue to next iteration in for loop
	}

	if ((strpos2 < 0)||(strpos2 >= char_array_length)) {
		//continue; //continue to next iteration in for loop
	}
	
*/
	else {
        	printf("The following was read in [%d]: %s\n",count,receive);
 		printf("The size of the string is %d: \n",strlen(receive));
		
		unsigned char *ans1, *ans2;
 		unsigned char *receive_pointer = receive;
		int strpos1,strpos2;
 	
 		ans1 = strstr(receive,"i="); //search for prefix in the receive string
		strpos1 = ans1 - receive_pointer;     //subtract pointers
      		printf("strpos1 = %d\n",strpos1);
	  
	  	if (strpos1 >= char_array_length) printf("Prefix (i=) not found\n");
	  
	  	if (strpos1 < char_array_length){	 
          	    str4[0] = receive[strpos1+2];
	   	    str4[1] = receive[strpos1+3];
	    	    str4[2] = receive[strpos1+4];
		    str4[3] = receive[strpos1+5];
		    str4[4] = receive[strpos1+6];
		    str4[5] = receive[strpos1+7];
		    str4[6] = receive[strpos1+8];
		    //str4[7] = receive[strpos+9];
        	    val = atof(str4);  //Convert string to floating point and add previous value
		    printf("Battery current is: %f\n",val);
		    i = val;
	  	}//end if

 	  	ans2 = strstr(receive,"V="); //search for prefix in the receive string
      	  	strpos2 = ans2 - receive_pointer;      //subtract pointers
      	  	printf("strpos2 = %d\n",strpos2);

		  if (strpos2 >= char_array_length) printf("Prefix (V=) not found\n");
	  
		  if (strpos2 < char_array_length){
		    str4[0] = receive[strpos2+2];
	   	    str4[1] = receive[strpos2+3];
	    	    str4[2] = receive[strpos2+4];
		    str4[3] = receive[strpos2+5];
		    str4[4] = receive[strpos2+6];
		    str4[5] = receive[strpos2+7];
		    str4[6] = receive[strpos2+8];
		    //str4[7] = receive[strpos+9];
		    val = atof(str4);
		    printf("Voltage is: %f\n",val);
		    z_ukf = val;
	  	}//end if
   	
		//************code inserted above********************************

		//usleep(0);   //delay in microseconds
		gettimeofday(&tv, NULL);
		t = tv.tv_sec + tv.tv_usec/1000000.0 - tzero;   //seconds
		dt = t-tminus1;
		//        what happens when time goes to 0.999 and then restarts at 0.000?
		//printf("The dt is %lf\n",dt);
		
		UKF_bhm3(i, &Param, z_ukf, state_ukf,
  	            w_ukf, est_ukf, P, 
   	            Q,  R,  dt,  ndim,  zdim,
   	           &ukfout); 	//JCS 4.5.16 modified from UKF_bhm3.c

		P = ukfout.P;     
		state_ukf = ukfout.state_ukf;   //first 3 elements = "before" column vector, next 3 are "after"
		est_ukf = ukfout.est_ukf;
		SOC = ukfout.SOC;
		printf("SOC is %f\n",SOC);
		//qb = state_ukf[3];   //SOC and t_EOD directly affected by errors in "i"
		qb = est_ukf[3];
 
		t_EOD = (q_EOD-qb)/(-i);
		t_EOD = t_EOD+t;

		x[(ndim)*(k)] = state_ukf[3];
		x[(ndim)*(k)+1] = state_ukf[4];
		x[(ndim)*(k)+2] = state_ukf[5];

		x_est[(ndim)*(k)] = est_ukf[3];
		x_est[(ndim)*(k)+1] = est_ukf[4];
		x_est[(ndim)*(k)+2] = est_ukf[5];

		fprintf(fp,"%f, %f, %f, %f, %f, %f\n",t,dt, i, z_ukf, SOC, t_EOD);
		printf("Battery voltage is %f\n",z_ukf);
		
		tminus1 = t; //Move current time to past time slot

	}//end else
     }//end for	
   
   close(file);
 
 fclose(fp);
 
}//end of main function

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
