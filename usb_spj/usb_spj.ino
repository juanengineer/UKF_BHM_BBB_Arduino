// This function is called once when the program starts. 
void setup() {
   // Choose a baud rate and configuration. 9600 
   // Default is 8-bit with No parity and 1 stop bit
   Serial.begin(9600, SERIAL_8N1);

}

// This function will loop as quickly as possible, forever.
void loop() {
   char charIn;
   /*
   Serial.println("led1on");
   delay(100);
   Serial.println("led1off");
   delay(100);
   */
   int c = 0; //total characters in buffer
   if(Serial.available()>0){          // A byte has been received
      /*
      
      charIn = Serial.read();       // Read the character in from the BBB
      delay(10);
      if (charIn == 'y') Serial.println("u");         // Send the character back to the BBB 
      if (charIn != 'y') Serial.println(charIn);
   }
   
   //else charIn = 'x1';
   //if(Serial.available()<=0) Serial.println("x2"); 
   */
   /*
   c = Serial.available();
   Serial.println(c);
   
   */
   charIn = Serial.read();
   delay(1000);
   if (charIn == 'y') Serial.println("i=59.0000.V=21.0000");
   if (charIn != 'y') Serial.println(charIn);   
  }
  if (Serial.available()<=0){
    //Serial.println("No data sent/received or port not open");
    delay(500);
  }
}
