clear variables;
close all;

parameters_EM_prevalence;
    
% get random pooling matrix
F = get_random_observation_2(M,N,N_G,N_O);
% generate random patients
X = get_random_patient(N,K);

% true state of pools
Y_0 = (F*X>0);
    
% observation 
Y = get_observation(Y_0,p_TP,p_FP);

% variables for BP
param.Y = Y;
param.F = F;
param.dmp = dmp;
param.THETA = THETA_E;
param.BP_STEP_MAX = BP_STEP_MAX;
    
param.rhoh = rhoh;
param.ph_TP = p_TP;
param.ph_FP = p_FP;

shusoku = 0;
count_rho = 0;

while shusoku == 0
        
    count_rho = count_rho+1;
    rhoh_old = rhoh;
        
    [prob_infect, rhoh, P_neg, conv_check] = GT_BP_EM(param);
        
    shusoku = (abs(rhoh_old-rhoh)<THETA_M);  
    param.rhoh = rhoh;
        
    if(count_rho>rho_STEP_MAX)
        conv_check_rho = 0;
        break;
    end

end
    
