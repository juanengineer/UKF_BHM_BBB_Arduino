  float a=0;  float P=1; float Q=1; float R=1; int k = 1;
  float tau = 100*pow(10,3)*22*pow(10,-6);
  const int xa0 = 21;
  float xa = xa0;
  float x = xa0;
  float x_est;
  float z;
  int row = 0;
  
void setup(){
  Serial.begin(9600);
  Serial.println("CLEARDATA");
  Serial.println("LABEL,Timer,truth,estimate");
}

void loop(){
  k++;
  a = -xa0/tau*pow(2.71828,-k/tau);
  x = a + x + random(0.1);
  z = x + random(0.1);
  x_est = xa;
  
  const float W0 = 0.99;
  float SQRTM_P = sqrt (2*P/(1-W0));
  float X0 = xa;
  float X1 = xa + SQRTM_P;
  float X2 = xa - SQRTM_P;
  float W1 = (1-W0)/(2*2);
  float W2 = W1;
  float X0F = a + X0;
  float X1F = a + X1 + SQRTM_P;
  float X2F = a + X2 - SQRTM_P;
  float XF_mean = W0*X0F + W1*X1F + W2*X2F ;
  float PX =  W0*(X0F - XF_mean)*(X0F - XF_mean) + W1*(X1F - XF_mean)*(X1F - XF_mean) +  W2*(X2F - XF_mean)*(X2F - XF_mean) + Q;
  float Z0F = X0F;  
  float Z1F = X1F;  
  float Z2F = X2F;
  float ZF_mean = W0*Z0F + W1*Z1F + W2*Z2F;
  float PZ = W0*(Z0F - ZF_mean)*(Z0F - ZF_mean) + W1*(Z1F - ZF_mean)*(Z1F - ZF_mean) +  W2*(Z2F - ZF_mean)*(Z2F - ZF_mean) + R; 
  float PXZ = W0*(X0F - XF_mean)*(Z0F - ZF_mean) + W1*(X1F - XF_mean)*(Z1F - ZF_mean) + W2*(X2F - XF_mean)*(Z2F - ZF_mean);
  float K = PXZ/PZ;
  xa = XF_mean + K*(z - ZF_mean);  
  P = PX - K*PZ*K;
  /*
  Serial.print("The truth state is: ");
  Serial.println(x);
  Serial.print("The state estimate is: ");
  Serial.println(xa);
  Serial.print("DATA,TIMER,"); Serial.println(pressure_mmHg);
  */
  Serial.print("DATA,TIMER,"); Serial.print(x);Serial.print(","); Serial.println(xa);
  row++;
  delay(250);
}
