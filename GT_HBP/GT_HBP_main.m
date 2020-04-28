clear variables;
close all;

parameters_HBP;

% get random pooling matrix  
F = get_random_observation_2(M,N,N_G,N_O);

% get random patients
X = get_random_patient(N,K);
   
% true state of pools
Y_0 = (F*X>0);
    
% observation 
Y = get_observation(Y_0,p_TP,p_FP);
    
% variables for HBP
param.beta_a = beta_a;
param.beta_b = beta_b;
param.Y = Y;
param.F = F;
param.p_TP = p_TP;
param.p_FP = p_FP;
param.dmp = dmp;
param.THETA = THETA;
param.BP_STEP_MAX = BP_STEP_MAX;

% calculate infection probability by BP
[prob_infect, rho_dist, conv_check] = GT_HBP(param);

X_MAP = (prob_infect>0.5);
TP = mean((X>0).*X_MAP')/mean(X);
FP = mean((X==0).*X_MAP')/mean(1-X);
rhoh = mean(rho_dist);
