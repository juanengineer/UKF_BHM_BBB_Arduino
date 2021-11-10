  float P=1; float Q=1; float R=1; int k=1;
  //float tau = 100*pow(10,3)*220*pow(10,-6);
 // float tau_model = 3300*pow(10,-6)*180;
 float tau_model = 3300*pow(10,-6)*10000*7;
//  const float xa0 = 3.3;
  const float xa0 = 5.0;
  float a = xa0;
  
 // const float xa0 = 5;
 // float xa = xa0;
  float xa = xa0;
  float x = xa0;
  float x_est;
  float z;
  int row = 0;
  int time_initial;
  float time;
  float t_remain;
  float t_predict;
  float dt = 0.04;
  float time_old;
  
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
  Serial.begin(115200);
 // Serial.println("CLEARDATA");
 // Serial.println("LABEL,Timer,truth,estimate");
 // Serial.println("LABEL,Timer,model,estimate,measurement,time to voltage,time const");
  pinMode(A0, OUTPUT);
  digitalWrite(A0,HIGH);
// analogWrite(9,int(3.72/5.0*255)); // Uses PWM to give an average of 3.3V (or value of xa0) to digital pin 9
// analogWrite(9,int(5.0/5.0*255)); // Uses PWM to give an average of 3.3V (or value of xa0) to digital pin 9
// Serial.println(int(xa0/5.0*255));
  delay(1000);
  
  pinMode(A0, INPUT);
  
 // delay(1000);
  time_initial = millis();
  time = time_initial;
  
 // analogWrite(9,0);            // Turns off digital pin 9 output
 // dt = 0.04;
}

void loop(){
//  k++;


 // a = xa0*(1-exp(-k/tau_model));
 // a = xa0*exp(-k*.250/tau_model);
 //time_old = time;
 time = millis() - time_initial;
 time = (float)time/1000;
 dt = time-time_old;
 //dt = 0.04;
 
 a = dt/tau_model;
 
 a = x/(1+a);
// a = (1/(1+dt/tau_model))*x;
// a = a/(1.0+dt/tau_model);
 //a = (1.0/(1.0+dt/tau_model))*a;
// a = xa0*exp(-time/tau_model);
// a = 1;
// x = a + random(0.005);
   x = a;
 
 // z = x + random(0.1);
  z = (float)(analogRead(A0))/1024*5;
 // x_est = xa;
  
  const float W0 = 0.999;
 // const float W0 = -0.1;
  const float W1 = (1-W0)/(2*2);
  const float W2 = W1;
  
  SQRTM_P = sqrt (2*P/(1-W0));
  X0 = xa;
  X1 = xa + SQRTM_P;
  X2 = xa - SQRTM_P;
  
  //a = (1/(1+dt/tau_model))*X0;
  a = 1 + dt/tau_model;
  a = X0/a;
  X0F = a;
  
  //a = (1/(1+dt/tau_model))*X1;
 // X1F = a + X1 + SQRTM_P;
  a = 1 + dt/tau_model;
  a = X1/a;
  X1F = a;
  
  //a = (1/(1+dt/tau_model))*X2;
  //X2F = a + X2 - SQRTM_P;
  a = 1 + dt/tau_model;
  a = X2/a;
  X2F = a;
  
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
 
 //   tau_model = -time/log(xa/xa0);     
 //   t_predict = -tau*log(0.20);  //Time from max voltage to 20% of max voltage
 //   t_predict = -tau_model*log(0.03/5.0);  
 
// Serial.print("DATA,TIMER,"); Serial.print(t_predict);Serial.print(","); Serial.println(z);
  Serial.print("time,z,"); Serial.print(time); Serial.print(","); Serial.println(z);
  //Serial.print(",");  Serial.print(xa);Serial.print(",");Serial.print(z); Serial.print(",");Serial.print(t_predict); Serial.print(",");
  //Serial.println(tau_model);
//  row++;
  delay(10);
}
