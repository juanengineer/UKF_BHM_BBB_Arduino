  float a=0;  float P=1; float Q=1; float R=1; int k = 1;
  //float tau = 100*pow(10,3)*220*pow(10,-6);
  float tau = 3300*pow(10,-6)*180;
  const float xa0 = 3.3;
 // const float xa0 = 5;
 // float xa = xa0;
  float xa = 1;
  float x = xa0;
  float x_est;
  float z;
  int row = 0;
  int time_initial;
  float time;
  float t_remain;
  float t_predict;
  
  //const float W0;
  float SQRTM_P;
  float X0;
  float X1;
  float X2;
  //const float W1;
  //const float W2;
  float X0F;
  float X1F;
  float X2F;
  float XF_mean;
  float PX;
  float Z0F;  
  float Z1F;  
  float Z2F;
  float ZF_mean;
  float PZ; 
  float PXZ;
  float K;
  
void setup(){
  Serial.begin(9600);
  Serial.println("CLEARDATA");
 // Serial.println("LABEL,Timer,truth,estimate");
  Serial.println("LABEL,Timer,model,estimate,measurement,time to voltage,time const");
//  pinMode(A0, OUTPUT);
//  digitalWrite(A0,HIGH);
 analogWrite(9,int(3.72/5.0*255)); // Uses PWM to give an average of 3.3V (or value of xa0) to digital pin 9

// Serial.println(int(xa0/5.0*255));
  
  pinMode(A0, INPUT);
  
  delay(1000);
  time_initial = millis();
  time = time_initial;
  
  analogWrite(9,0);            // Turns off digital pin 9 output

}

void loop(){
//  k++;

 // a = xa0*(1-exp(-k/tau));
 // a = xa0*exp(-k*.250/tau);
 time = millis() - time_initial;
 time = (float)time/1000;
 a = xa0*exp(-time/tau);
 x = a + random(0.005);
 
 // z = x + random(0.1);
  z = (float)(analogRead(A0))/1024*5;
//  x_est = xa;
  
  const float W0 = -0.1;
  SQRTM_P = sqrt (2*P/(1-W0));
  X0 = xa;
  X1 = xa + SQRTM_P;
  X2 = xa - SQRTM_P;
  const float W1 = (1-W0)/(2*2);
  const float W2 = W1;
  X0F = a + X0;
  X1F = a + X1 + SQRTM_P;
  X2F = a + X2 - SQRTM_P;
  XF_mean = W0*X0F + W1*X1F + W2*X2F ;
  PX =  W0*(X0F - XF_mean)*(X0F - XF_mean) + W1*(X1F - XF_mean)*(X1F - XF_mean) +  W2*(X2F - XF_mean)*(X2F - XF_mean) + Q;
  Z0F = X0F;  
  Z1F = X1F;  
  Z2F = X2F;
  ZF_mean = W0*Z0F + W1*Z1F + W2*Z2F;
  PZ = W0*(Z0F - ZF_mean)*(Z0F - ZF_mean) + W1*(Z1F - ZF_mean)*(Z1F - ZF_mean) +  W2*(Z2F - ZF_mean)*(Z2F - ZF_mean) + R; 
  PXZ = W0*(X0F - XF_mean)*(Z0F - ZF_mean) + W1*(X1F - XF_mean)*(Z1F - ZF_mean) + W2*(X2F - XF_mean)*(Z2F - ZF_mean);
  K = PXZ/PZ;
  xa = XF_mean + K*(z - ZF_mean);  
  P = PX - K*PZ*K;
  /*
  Serial.print("The truth state is: ");
  Serial.println(x);
  Serial.print("The state estimate is: ");
  Serial.println(xa);
  Serial.print("DATA,TIMER,"); Serial.println(pressure_mmHg);
  */
  //Serial.print("DATA,TIMER,"); Serial.print(x);Serial.print(","); Serial.println(xa);
 
    tau = -time/log(xa/xa0);     
 //   t_predict = -tau*log(0.20);  //Time from max voltage to 20% of max voltage
    t_predict = -tau*log(0.03/5.0);  
 
// Serial.print("DATA,TIMER,"); Serial.print(t_predict);Serial.print(","); Serial.println(z);
  Serial.print("DATA,"); Serial.print(time); Serial.print(","); Serial.print(x);Serial.print(",");
  Serial.print(xa);Serial.print(",");Serial.print(z); Serial.print(",");Serial.print(t_predict); Serial.print(",");
  Serial.println(tau);
  row++;
//  delay(1);
}
