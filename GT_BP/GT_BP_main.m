clear variables;
close all;

% load parameters
parameters;

% generate pooling matrix
F = get_random_observation_2(M,N,N_G,N_O);

% generate true state of patients
X = get_random_patient(N,K);

% true state of pools
Y_0 = (F*X>0);
    
% observation 
Y = get_observation(Y_0,p_TP,p_FP);

% variables for BP
param.Y = Y;
param.F = F;
param.p_TP = p_TP;
param.p_FP = p_FP;
param.rhoh = rhoh;
param.dmp = dmp;
param.THETA = THETA;
param.BP_STEP_MAX = BP_STEP_MAX;

% calculate infection probability by BP
[prob_infect, conv_check] = GT_BP(param);

% MAP estimator
X_MAP = (prob_infect > 0.5);

% true positive rate
TP = mean((X>0).*X_MAP')/mean(X);

% false positive rate
FP = mean((X==0).*X_MAP')/mean(1-X);

    
