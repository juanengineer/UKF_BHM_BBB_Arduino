/** Serial communication test program that sends data to the Arduino echo
* program and reads the response on UART4. 
* 
* Written by Derek Molloy for the book "Exploring BeagleBone: Tools and 
* Techniques for Building with Embedded Linux" by John Wiley & Sons, 2014
* ISBN 9781118935125. Please see the file README.md in the repository root 
* directory for copyright and GNU GPLv3 license information.            */


#include<stdio.h>
#include<fcntl.h>
#include<unistd.h>
#include<termios.h>   // using the termios.h library
#include<string.h>
#include<stdlib.h>
#include <math.h>

// function declarations
//float UKF_state_estimate(float x, float x_est, float z, float P);
//float P_covar_calc(float PX, float K, float PZ);
//float x_est = 0;  //Be careful wth global variables
//float P = 0;	  //Be careful wth global variables
//float time = 0;	  //Be careful wth global variables

struct updates{
	float x;
	float x_est;
	float z;
	float time;
	float P;
};

struct updates UKF_state_estimate(struct updates new, struct updates old);
//updates UKF_state_estimate(UKF_input);

int main(){
   struct updates old;
   struct updates new;
   struct updates ukf_output;
   int file, count;
   int res, strpos;
   
   char str4[25];
   float val;

 for( ; ; ){
   
   if ((file = open("/dev/ttyUSB0", O_RDWR | O_NOCTTY | O_NDELAY))<0){
      perror("USB: Failed to open the file.\n");
      return -1;
   }
   struct termios options;               //The termios structure is vital
   tcgetattr(file, &options);            //Sets the parameters associated with file

   // Set up the communications options:
   //   9600 baud, 8-bit, enable receiver, no modem control lines
   options.c_cflag = B115200 | CS8 | CREAD | CLOCAL;
   options.c_iflag = IGNPAR | ICRNL;    //ignore partity errors, CR -> newline
   tcflush(file, TCIFLUSH);             //discard file information not transmitted
   tcsetattr(file, TCSANOW, &options);  //changes occur immmediately


   usleep(100000);                  //give the Arduino a chance to respond

	old = ukf_output;
		

   unsigned char receive[400];      //declare a buffer for receiving data
	
   if ((count = read(file, (void*)receive, 400))<0){   //receive the data
      perror("Failed to read from the input\n");
      return -1;
   }
   if (count==0) printf("There was no data available to read!\n");
   else {
      printf("The following was read in [%d]: %s\n",count,receive);
      printf("The size of the string is %d: \n",strlen(receive));
      

	unsigned char *ans1;
	unsigned char *receive_pointer = receive;
	  ans1 = strstr(receive,"Total Battery Voltage = "); //search for prefix in the receive string
      strpos = ans1 - receive_pointer;      //subtract pointers
      if (strpos >= 400) printf("Prefix not found\n");
	 
	  if (strpos < 400){
	  printf("strpos = %d\n",strpos);
 /*     
          str4[0] = receive[strpos+24];
	  str4[1] = receive[strpos+25];
	  str4[2] = receive[strpos+26];
	  str4[3] = receive[strpos+27];
	  str4[4] = receive[strpos+28];
	  str4[5] = receive[strpos+29];
	  str4[6] = receive[strpos+30];
	  str4[7] = receive[strpos+31];

	 
      new.time = atof(str4);  //Convert string to floating point and add previous value
 //     printf("Convert to float and add2: %f\n",val+2);
*/	 

	  str4[0] = receive[strpos+24];
	  str4[1] = receive[strpos+25];
	  str4[2] = receive[strpos+26];
	  str4[3] = receive[strpos+27];
	  str4[4] = receive[strpos+28];
	  str4[5] = receive[strpos+29];
	  str4[6] = receive[strpos+30];
	  str4[7] = receive[strpos+31];
	  val = atof(str4);
	  printf("Voltage is: %f\n",val);

/*	  new.z = atof(str4); 

	  ukf_output = UKF_state_estimate(new,old);
	  printf("The new voltage estimate is: %f\n",ukf_output.x_est);
*/
	  }

	
}
   close(file);
 }
   return 0;
}

// UKF function that estimates state

struct updates UKF_state_estimate(struct updates new, struct updates old){

	struct updates ukf_output;
	const float tau_model = 3300*pow(10,-6)*10000*7; 
	const float W0 = 0.999;
	const float W1 = (1-W0)/(2*2);
	const float W2 = W1;
//	float xa = x_est;
	float xa = old.x_est;
	float X0; float X1; float X2;
	float X0F; float X1F; float X2F; float XF_mean;
	float Z0F; float Z1F; float Z2F; float ZF_mean; 
	float PX; float PZ; float PXZ;
	float K;
	float Q = 1; float R = 1;
//	float dt = 0.04;
	float dt = new.time - old.time;
	
	float a; 
	float SQRTM_P;


	a = old.x/(1 + dt/tau_model);  //Physics model	
	ukf_output.x = a;	

	SQRTM_P = sqrt(2*old.P/(1-W0));
	
	X0 = xa;
	X1 = xa + SQRTM_P;
	X2 = xa - SQRTM_P;

	a = X0/(1 + dt/tau_model);
	X0F = a;

	a = X1/(1 + dt/tau_model);
	X1F = a;

	a = X2/(1 + dt/tau_model);
	X2F = a;

	XF_mean = W0*X0F + W1*X1F + W2*X2F;
	PX =  W0*(X0F - XF_mean)*(X0F - XF_mean) + W1*(X1F - XF_mean)*(X1F - XF_mean) +  W2*(X2F - XF_mean)*(X2F - XF_mean) + Q;
	
	Z0F = X0F;
	Z1F = X1F;
	Z2F = X2F;
	ZF_mean = W0*Z0F + W1*Z1F + W2*Z2F;
	PZ = W0*(Z0F - ZF_mean)*(Z0F - ZF_mean) + W1*(Z1F - ZF_mean)*(Z1F - ZF_mean) +  W2*(Z2F - ZF_mean)*(Z2F - ZF_mean) + R;
	PXZ = W0*(X0F - XF_mean)*(Z0F - ZF_mean) + W1*(X1F - XF_mean)*(Z1F - ZF_mean) + W2*(X2F - XF_mean)*(Z2F - ZF_mean);
	K = PXZ/PZ;
	xa = XF_mean + K*(new.z - ZF_mean);  
	ukf_output.x_est = xa;
	ukf_output.P = PX - K*PZ*K;

	
	return ukf_output;
}

