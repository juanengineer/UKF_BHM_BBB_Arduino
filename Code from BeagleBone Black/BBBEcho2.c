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

int main(){
   int file, count;
   int res, strpos;
   unsigned char *ans1;
   char str4[25];
   float val;

   if ((file = open("/dev/ttyO4", O_RDWR | O_NOCTTY | O_NDELAY))<0){
      perror("UART: Failed to open the file.\n");
      return -1;
   }
   struct termios options;               //The termios structure is vital
   tcgetattr(file, &options);            //Sets the parameters associated with file

   // Set up the communications options:
   //   9600 baud, 8-bit, enable receiver, no modem control lines
   options.c_cflag = B9600 | CS8 | CREAD | CLOCAL;
   options.c_iflag = IGNPAR | ICRNL;    //ignore partity errors, CR -> newline
   tcflush(file, TCIFLUSH);             //discard file information not transmitted
   tcsetattr(file, TCSANOW, &options);  //changes occur immmediately


   usleep(100000);                  //give the Arduino a chance to respond

   unsigned char receive[300];      //declare a buffer for receiving data
   if ((count = read(file, (void*)receive, 300))<0){   //receive the data
      perror("Failed to read from the input\n");
      return -1;
   }
   if (count==0) printf("There was no data available to read!\n");
   else {
      printf("The following was read in [%d]: %s\n",count,receive);
      printf("The size of the string is %d: \n",strlen(receive));
      ans1 = strstr(receive,"C1:"); //search for "Voltage" in the receive string
      strpos = ans1 - receive;      //subtract pointers
      printf("strpos = %d\n",strpos);
      str4[0] = receive[strpos+3];
      val = atof(str4);  //Convert string to floating point
      printf("Convert to float and add2: %f\n",val+2);
	
      
   }
   close(file);
   return 0;
}
