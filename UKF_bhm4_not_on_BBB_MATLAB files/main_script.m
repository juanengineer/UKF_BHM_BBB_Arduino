function [] = main_script()

i_frombattery_model = horzcat(60*ones(1,57),25*ones(1,265),35*ones(1,228),...
    20*ones(1,142),1*ones(1,50),30*ones(1,10),60*ones(1,1),1*ones(1,30),...
    20*ones(1,20),1*ones(1,82),21*ones(1,347),13*ones(1,193),1*ones(1,10));

%i_frombattery_model = i_frombattery_model + 0*randn(1,N);
%i_frombattery_model = horzcat(60*ones(1,5));
%i_frombattery_model = [60 50 40 30 20];


N = length(i_frombattery_model);

%i_meas = i_frombattery_model + 0.1*randn(1,N);
%N = 30;                     %Time steps
%ndim = 4;         %Length of state vector
%ndim = 3;         %Removed Tb from state and meas vectors. Added battery current
ndim = 4;         %Added battery current
%zdim = 1;        %Length of meas vector
zdim = 2;        %Added current sensor value
x = ones(ndim,N); 
%{
 

xa0 = 3.3;
x(1,1) = xa0;
x(2,1) = tau_model;
 
v = 0.05*randn(2,N); 
R = 0.04*eye(2);
%}
%w = 0.05*randn(ndim,N);
w = 0.02*ones(ndim,N);
Q = 1*eye(ndim);
R = 1*eye(zdim);
%R = 0.04*eye(ndim);
dt = 1;
%dt = 0.04;

qmax = 2.88e+4;    
Cmax = 2.85e+4;
Ccb0 = 19.4;
Ccb1 = 1576;
Ccb2 = 41.7;
Ccb3 = -203;
Rs = 2.77e-2;
Cs = 89.3;
Rcp0 = 1.6e-3;
Rcp1 = 8.45;
Rcp2 = -61.9;
Ccp0 = 2689;
Ccp1 = -2285;
Ccp2 = -0.73;
Rp = 1e5;
Jt = 800;
hcp = 19;
hcs = 1;
ha = 0.5;   %heat transfer coeff
Ta = 25;    %ambient temp degC

Param = struct('qmax',qmax,'Cmax',Cmax,'Ccb0',Ccb0,'Ccb1',Ccb1,...
    'Ccb2',Ccb2,'Ccb3',Ccb3,'Rs',Rs,'Cs',Cs,'Rcp0',Rcp0,'Rcp1',Rcp1,'Rcp2',Rcp2,...
    'Ccp0',Ccp0,'Ccp1',Ccp1,'Ccp2',Ccp2,'Rp',Rp,'Jt',Jt,'hcp',hcp,'hcs',hcs,...
    'ha',ha,'Ta',Ta);

%[SOC,x_total,t_EOD, Voltage, Cbout] = generate_measurements(N,dt,i_frombattery);
z = ones(ndim,N);
%[z,t_EOD, SOC_EOD, counter_array] = generate_measurements(N,dt,i_frombattery_model);
[z] = generate_measurements(N,dt,i_frombattery_model);

%z = zeros(2,N);
x_est = ones(ndim,N);
%{
t_to_volt = zeros(1,N);

%P = 1.1*eye(2);
%z_ukf = [1 ; 1];
%}
state_ukf = 1.0*ones(ndim,2); %First column is old, second is new
est_ukf = 1.0*ones(ndim,2); %First column is old estimate
P = 1*eye(ndim);
%z_ukf = [1 ; 1];            %Shape of meas. vector.
z_ukf = [1;
         1];
%x(:,1) = [qmax qmax qmax Ta]'; %Initial state condition
%x(:,1) = [qmax qmax qmax 58.5]'; %Initial state condition
x(:,1) = [28800 0.000534 -81.9523 58.5]'; %Initial state condition
SOC_array = ones(1,N);
SOC_array(1) = 1;


SOC_EOD_array(1) = 0.999;
%Section below has to do with discharge time prediction
%SOC_EOD = 0.28;
%Below is Newton method
%counter = 0;
voltage_knee = 16.7;
%SOC_EOD = SOC_EOD_array(k-1);
SOC_EOD = 0.999;
newton_error = 1;
while newton_error > 0.01 %This is fractional error
    Cb = Ccb0 + Ccb1*SOC_EOD+Ccb2*SOC_EOD^2+Ccb3*SOC_EOD^3;
    Cb_prime = Ccb1 + 2*Ccb2*SOC_EOD + 3*Ccb3*SOC_EOD^2;
    qb_knee = voltage_knee*Cb;
    f_of_SOC = 1-(qmax-qb_knee)/Cmax-SOC_EOD;
    f_prime = -voltage_knee/Cmax*Cb_prime;
    SOC_EOD = SOC_EOD - f_of_SOC/f_prime;
    newton_error = abs((f_of_SOC/f_prime)/SOC_EOD);
 %   counter = counter + 1;
end
Cb;
Cb_prime;
qb_knee;
f_of_SOC;
f_prime;
SOC_EOD;
newton_error;
q_EOD = qmax - Cmax*(1-SOC_EOD)
disp('The q_EOD is above')
t_EOD_sum = 0;

for k = 2:N     %UKF routine
state_ukf(:,1) = x(:,k-1);
w_ukf = w(:,k);
est_ukf(:,1) = x_est(:,k-1);

%i = i_frombattery_model(k);
z_ukf = z(:,k);  



%[P,state_ukf,est_ukf, SOC] = UKF_bhm2(i,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim);

%ukfout = UKF_bhm3(i,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim);

%i_model = 58.5;  %This comes from the battery current load profile that is available
i_model = i_frombattery_model(k);
ukfout = UKF_bhm3(i_model,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim); 

%The battery current i, will be output as an estimated state and updated using current sensor measurements 

P = ukfout.P;
state_ukf = ukfout.state_ukf;
est_ukf = ukfout.est_ukf;
SOC = ukfout.SOC;
i = est_ukf(4,2);   %Updated battery current estimate

%P = P;        %Reminder of P update

%qb = x_est(1,k);
disp('true qb below');
%qb = state_ukf(1,2);

disp('est qb below');
qb = est_ukf(1,2);


x(:,k) = state_ukf(:,2);
x_est(:,k) = est_ukf(:,2);
SOC_array(k)= SOC;
SOC_actual_array(k) = 1 - (qmax-qb)/Cmax;
i_array(1,k) = x_est(4,k);

%t_EOD(k) = (t_EOD(k) + t_EOD(k-1))/2; %take average
t_present = k*dt;

t_EOD(k) = (q_EOD-qb)/(-i) + t_present;  

t_EOD_sum = t_EOD_sum + t_EOD(k);
t_avg(k) = t_EOD_sum/k;

end
plot_results(t_avg,i_array,SOC_actual_array,i_frombattery_model,N,x_total,z,x_est,dt,t_EOD,SOC_EOD, counter_array)
%plot_results(SOC_array,i_frombattery_model,N,x_total,z,x_est,dt,t_EOD,SOC_EOD_array, counter_array)
%plot_results(SOC_array,i_frombattery_model,N,x_total,z,x_est,dt,t_EOD)
%plot_results(SOC_array,i_frombattery_model,N,x_total,z,x_est,dt)
%plot_results(N,x_total,SOC,i_frombattery_model,dt)
%plot_results(N,x,z,x_est,dt,t_to_volt);

function [z] = generate_measurements(N,dt,i_frombattery_model)
%function [z,t_EOD, SOC_EOD, counter_array] = generate_measurements(N,dt,i_frombattery_model)
SOC = ones(1,N); 
state_meas = ones(2,1);
%z(1,1) = xa0;
qmax = 2.88e+4;    
Cmax = 2.85e+4;
Ccb0 = 19.4;
Ccb1 = 1576;
Ccb2 = 41.7;
Ccb3 = -203;
Rs = 2.77e-2;
Cs = 89.3;
Rcp0 = 1.6e-3;
Rcp1 = 8.45;
Rcp2 = -61.9;
Ccp0 = 2689;
Ccp1 = -2285;
Ccp2 = -0.73;
Rp = 1e5;
Jt = 800;
hcp = 19;
hcs = 1;
ha = 0.5;   %heat transfer coeff
Ta = 25;    %ambient temp degC

%x = ones(ndims,1);

%x(:,1) = [qmax qmax qmax Ta]'; %Initial state condition
%x(:,1) = [qmax qmax qmax 58.5]'; %Initial state condition
x(:,1) = [28800 0.000534 -81.9523 58.5]'; %Initial state condition

SOC(1) = 100;

x_total(:,1) = x(:,1);

SOC_EOD_array(1) = 0.999;
counter_array(1) = 0;

%Section below has to do with discharge time prediction
%SOC_EOD = 0.28;
%Below is Newton method
%counter = 0;
voltage_knee = 16.7;
%SOC_EOD = SOC_EOD_array(k-1);
SOC_EOD = 0.999;
newton_error = 1;
while newton_error > 0.01 %This is fractional error
    Cb = Ccb0 + Ccb1*SOC_EOD+Ccb2*SOC_EOD^2+Ccb3*SOC_EOD^3;
    Cb_prime = Ccb1 + 2*Ccb2*SOC_EOD + 3*Ccb3*SOC_EOD^2;
    qb_knee = voltage_knee*Cb;
    f_of_SOC = 1-(qmax-qb_knee)/Cmax-SOC_EOD;
    f_prime = -voltage_knee/Cmax*Cb_prime;
    SOC_EOD = SOC_EOD - f_of_SOC/f_prime;
    newton_error = abs((f_of_SOC/f_prime)/SOC_EOD);
 %   counter = counter + 1;
end
Cb;
Cb_prime;
qb_knee;
f_of_SOC;
f_prime;
SOC_EOD;
newton_error;
q_EOD = qmax - Cmax*(1-SOC_EOD);
%SOC_EOD_array(k) = SOC_EOD;
%counter_array(k) = counter;

for k = 2:N

qb = x(1);
qcp = x(2);
qcs = x(3);
%Tb = x(4);

SOC(k) = 1 - (qmax-qb)/Cmax;
%SOC(k) = SOC(k)*100;
Cb = Ccb0 + Ccb1*SOC(k)+Ccb2*SOC(k)^2+Ccb3*SOC(k)^3;
Ccp = Ccp0 + Ccp1*exp(Ccp2*SOC(k));
Rcp = Rcp0 + Rcp1*exp(Rcp2*SOC(k));
%{    
xB = [qb qcp qCs]';
A = [-1/(Cb*Rp) 1/(Ccp*Rp) 1/(Cs*Rp);
    1/(Cb*Rp) -1/(Ccp*Rp*Rcp) 1/(Cs*Rp);
    1/(Cb*Rp) 1/(Ccp*Rp) 1/(Cs*Rp)];

current_vector = current*ones(3,1);
state_meas = xB + (A*xB + current_vector + v(:,k))*dt; %3 x 1 vector
%}

%qhat = qmax - (qmax*0.99 - Cmax)*(1-SOC(k)); %Causes divergence

%Vb = qhat/Cb;   %Diverges
Vb = qb/Cb;    %Works with this!
Vcs = qcs/Cs;
Vcp = qcp/Ccp;
Vp = Vb - Vcp - Vcs;
%V(k) = Vp;

ip = Vp/Rp;
ib = ip + i_frombattery_model(k);
icp = ib - Vcp/Rcp;

qbdot = -ib;
qcpdot = icp;
qcsnew = ib*Rs*Cs;

%Tbdot = 1/Jt*(hcp*Vcp^2/Rcp + hcs*Vcs^2/Rs - ha*(Tb-Ta));
%Tbdot = 1;
%{
x = [ qb + qbdot*dt;
        qcp + qcpdot*dt;
        qcsnew;
        Tb + Tbdot*dt];
%}

x = [ qb + qbdot*dt;
        qcp + qcpdot*dt;
        qcsnew
        58.5];
  
x_total(:,k) = x; 
qb = x(1);
qcp = x(2);
qcs = x(3);

t_present = k*dt;
qb;
%t_prediction(k) = (q_prediction-qmax)/(-i_frombattery(k));   %Corrected eq!
t_EOD(k) = (q_EOD-qb)/(-i_frombattery_model(k));    



t_EOD(k) = t_EOD(k)+t_present;
%t_EOD(k) = (t_EOD(k) + t_EOD(k-1))/2; %take average
t_EOD_display = t_EOD(k);

z(1,k) = Vp;        %Measured battery voltage
z(2,k) = i_frombattery_model(k)+5*randn(1,1);    %Measured battery current


%z(2,k) = i_frombattery_model(k) + 0.1*randn(1,1);
%z(2,k) = 30;   %Expected temperature
%Voltage(k) = 0;


%z(:,k) = state_meas;
%{
a = (1/(1+dt/tau_model))*z(1,k-1); 
tau_meas = (1/(-1+z(1,k-1)/a))*dt;
%}
%state_meas = [a; tau_meas] + v(:,k);
%z(:,k) = state_meas;
end

end

function [] = plot_results(t_avg,i_array,SOC_actual_array,i_frombattery_model,N,x_total,z,x_est,dt,t_EOD,SOC_EOD, counter_array)
%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure  %Makes a new instance of charts
subplot
k = 1:N;
subplot(2,2,1)
%plot(k*dt,z(2,k),'b.',k*dt,i_array,'b--',k*dt,i_frombattery_model,'r',k*dt,SOC_array*100,'bo',k*dt,z(1,k),'g',k*dt,SOC_EOD*100*ones(1,N),'r--')
plot(k*dt,z(2,k),'b.',k*dt,i_array,'r.',k*dt,SOC_actual_array*100,'bo',k*dt,z(1,k),'g',k*dt,SOC_EOD*100*ones(1,N),'r--')

legend('current measure','batt current estimate','B3 SOC','Vp','SOC_{EOD}')
xlim([0 1435]);
ylim([0 100]);

subplot(2,2,2)
plot(k*dt,t_EOD,'k',k*dt,t_avg(k),'b')
legend('t_{EOD}','t_{avg}')
xlim([0 1435]);
ylim([0 2000]);

subplot(2,2,3)
plot(k*dt,z(1,k),'g')
legend('Vp')
xlim([0 1435]);
ylim([16 21]);

%{
subplot(2,2,4)
plot(k*dt,counter_array,'r*')
legend('iterations')
%ylim([0 5])
%}
end
end