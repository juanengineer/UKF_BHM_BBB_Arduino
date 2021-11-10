function [ukfout] = UKF_bhm3(i,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim)
%#codegen
%function [P,state_ukf,est_ukf, SOC] = UKF_bhm3(i,Param,z_ukf,state_ukf,w_ukf,est_ukf,P,Q,R,dt,ndim,zdim)
%function [t_predict,P,state_ukf,est_ukf] = UKF_bhm(z_ukf,tau_model,state_ukf,w_ukf,est_ukf,P,Q,R,dt,xa0,ndim)

% Function expanded to N dimensional state vector
% Array indices must begin at 1 and not 0
% Sigma points are denoted by the "Chi" variables

%%SOC = ones(1,N); 
%%state_meas = ones(2,1);

qmax = Param.qmax;     
Cmax = Param.Cmax;
Ccb0 = Param.Ccb0;
Ccb1 = Param.Ccb1;
Ccb2 = Param.Ccb2;
Ccb3 = Param.Ccb3;
Rs = Param.Rs;

Cs = Param.Cs;
Rcp0 = Param.Rcp0;
Rcp1 = Param.Rcp1;
Rcp2 = Param.Rcp2;
Ccp0 = Param.Ccp0;
Ccp1 = Param.Ccp1;
Ccp2 = Param.Ccp2;

Rp = Param.Rp;
Jt = Param.Jt;
hcp = Param.hcp;
hcs = Param.hcs;
ha = Param.ha;   %heat transfer coeff
Ta = Param.Ta;    %ambient temp degC


%SOC(1) = 100;

%%x_total(:,1) = x;
qb = state_ukf(1,1);
qcp = state_ukf(2,1);
qcs = state_ukf(3,1);
%Tb = state_ukf(4,1); %Removed Tb from state and meas vectors

SOC = 1 - (qmax-qb)/Cmax;
Cb = Ccb0 + Ccb1*SOC+Ccb2*SOC^2+Ccb3*SOC^3;
Ccp = Ccp0 + Ccp1*exp(Ccp2*SOC);
Rcp = Rcp0 + Rcp1*exp(Rcp2*SOC);

%qhat = qmax - (qmax*0.99 - Cmax)*(1-SOC);
%Vb = qhat/Cb; %Diverges

Vb = qb/Cb;     %Works with this! Only diverges at very end.
Vcs = qcs/Cs;
Vcp = qcp/Ccp;
Vp = Vb - Vcp - Vcs;

ip = Vp/Rp;
ib = ip + i;
icp = ib - Vcp/Rcp;

qbdot = -ib;
qcpdot = icp;
qcsnew = ib*Rs*Cs;
%Tbdot = 1/Jt*(hcp*Vcp^2/Rcp + hcs*Vcs^2/Rs - ha*(Tb-Ta));

x = [ qb + qbdot*dt;
        qcp + qcpdot*dt;
        qcsnew];
%{  
x = [ qb + qbdot*dt;
        qcp + qcpdot*dt;
        qcsnew;
        Tb + Tbdot*dt];    
%}
        
state_model = x + w_ukf;   %voltage and tau
state_ukf(:,2)= state_model;
xa = est_ukf(:,1);
%% Selection of Sigma Points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W0 = -0.8;  % -1 < W0 < 1
%W0 = 0.99;              %Weights for tau are different than for voltage
%W1 = (1-W0)/(2*2);
W = zeros(2*ndim+1,1);
W(1) = 0.99;                            %Sigma point array begins at 1 and not 0
SQRTM_P = real(sqrtm( 2*P/(1-W(1)) ));  %Square root of matrix is converted to real type
Chi = zeros(ndim,2*ndim+1);
Chi(:,1) = xa; % xa is old state estimate feed back from algorithm

for k = 2:ndim+1            %Generate sigma points with +sqrt. Increase indices by 1 since arrays don't start from 0
Chi(:,k) = xa + SQRTM_P(:,k-1); 
W(k) = (1-W(1))/(2*2);

Chi(:,ndim+k) = xa - SQRTM_P(:,k-1); 
W(ndim+k) = (1-W(1))/(2*2);
end
%% Model Forecast Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate Sigma Points through non-linear transform
ChiF = zeros(ndim,2*ndim+1);

for k = 1:2*ndim+1          %Forecast sigma points
ip = Vp/Rp;
ib = ip + i;        
Chi_qcsnew = ib*Rs*Cs;    %%Consider finding qcsdot

ChiF(:,k) = [ Chi(1,k) + qbdot*dt;    %Symbol for sigma pt is chi
        Chi(2,k) + qcpdot*dt;
        Chi_qcsnew];            %%% Not sure about this!!!

end

%% Compute the mean and covariance of forecast
ChiF_mean = zeros(ndim,1);

for k = 1:2*ndim+1          %Calculate mean of sigma points
ChiF_mean = ChiF_mean + W(k)*ChiF(:,k); 
end

PX = zeros(ndim); 
%PX = zeros(ndim,ndim);      
for k = 1:2*ndim+1    %Calculate covariance matrix of sigma pts
PX = PX + W(k)*(ChiF(:,k) - ChiF_mean)*(ChiF(:,k) - ChiF_mean)';
end

PX = PX+Q;
%% Propagate the sigma points through the observation model
%ZF = ChiF;                  %Expected measurements based on forecasted sigma points
ZF = [1/Cb -1/Ccp -1/Param.Cs]*ChiF; %Expected measurements based on forecasted sigma points

%% Compute mean and covariance of transformed observations
ZF_mean = zeros(zdim,1);
%ZF_mean = zeros(ndim,1);

for k = 1:2*ndim+1    %Calculate mean of expected measurements
    ZF_mean = ZF_mean + W(k)*ZF(:,k); 
end

PZ = zeros(zdim);
%PZ = zeros(ndim,ndim);
for k = 1:2*ndim+1    %Covar of expected measurements
PZ = PZ + W(k)*(ZF(:,k) - ZF_mean)*(ZF(:,k) - ZF_mean)';
end
PZ = PZ+R;
%% Compute the cross covariance between XF and Zf
PXZ = zeros(ndim,zdim);
%PXZ = zeros(ndim,ndim);

for k = 1:2*ndim+1    %Cross covariance of sigma points and expected measurements
PXZ = PXZ + W(k)*(ChiF(:,k) - ChiF_mean)*(ZF(:,k) - ZF_mean)'; 
end
%% Data Assimilation Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%K = PXZ*pinv(PZ);
K = PXZ*PZ^(-1);
xa = ChiF_mean + K*(z_ukf - ZF_mean);
P = PX - K*PZ*K';
est_ukf(:,2) = xa;
%t_predict = -xa(2,1)*log(0.5/xa0);
ukfout = struct('P',P,'state_ukf',state_ukf,'est_ukf',est_ukf,'SOC',SOC);

end
