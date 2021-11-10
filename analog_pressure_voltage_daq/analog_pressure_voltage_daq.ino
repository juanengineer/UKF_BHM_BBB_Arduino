int analogPin = 5;     // potentiometer wiper (middle terminal) connected to analog pin 3
int pressure_mmHg = 0;
int row = 0;
                       // outside leads to ground and +5V
int val = 0;           // variable to store the value read
//pressure(mmHg) = -0.4827*val + 273.73 

void setup()
{
  Serial.begin(128000);          //  setup serial
  analogReference(INTERNAL);
  Serial.println("CLEARDATA");
  Serial.println("LABEL,Timer,pressure,pulse");
}

void loop()
{
  old_pressure = pressure_mmHg;
  val = analogRead(analogPin);    // read the input pin
  pressure_mmHg = -0.4827*(float)val + 273.73; //truncated integer is OK
 // Serial.println(pressure_mmHg);  //to debug monitor
 
 Serial.print("DATA,TIMER,"); Serial.println(pressure_mmHg);
  row++;
  
/*  
  if (row > 360) 
   {
    row=0;
    Serial.println("ROW,SET,2");
   }
 */  
   delay(10);
  
}

