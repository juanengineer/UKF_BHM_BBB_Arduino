// This function is called once when the program starts. 
void setup() {
   // Choose a baud rate and configuration. 9600 
   // Default is 8-bit with No parity and 1 stop bit
   Serial.begin(9600, SERIAL_8N1);

}

// This function will loop as quickly as possible, forever.
void loop() {
   char charIn;
   if(Serial.available()>0){          // A byte has been received
      charIn = Serial.read();       // Read the character in from the BBB
      if (charIn == 'y') Serial.println("u");         // Send the character back to the BBB 
      if (charIn != 'y') Serial.println(charIn);
      delay(100);
   }
   else charIn = 'x1';
   if(Serial.available()<=0) Serial.println("x2"); 
   delay(100);
}
