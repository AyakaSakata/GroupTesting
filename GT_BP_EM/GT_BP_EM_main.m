clear variables;
close all;

parameters_EM;
    
% get random pooling matrix
F = get_random_observation_2(M,N,N_G,N_O);
% generate random patients
X = get_random_patient(N,K);

% true state of pools
Y_0 = (F*X>0);
    
% observation 
Y = get_observation(Y_0,p_TP,p_FP);

% Assumed values of parameters
ph_TP = p_TP-0.01;
ph_FP = p_FP+0.01;
    
% variables for BP
param.Y = Y;
param.F = F;
param.p_TP = p_TP;
param.p_FP = p_FP;
param.rhoh = rhoh;
param.dmp = dmp;
param.THETA = THETA_E;
param.BP_STEP_MAX = BP_STEP_MAX;

param.ph_TP = ph_TP;
param.ph_FP = ph_FP;
param.rhoh = rhoh;

shusoku = 0;
count_M = 0;
while shusoku == 0
        
    count_M = count_M+1;
        
    ph_TP_old = ph_TP;
    ph_FP_old = ph_FP;
    
    [prob_infect, rhoh, P_neg, conv_check] = GT_BP_EM(param);
        
    TPFP = zeros(2,1);
    TPFP(1) = ph_TP;
    TPFP(2) = ph_FP;
        
    U = ph_TP*Y+(1-ph_TP)*(1-Y);
    W = ph_FP*Y+(1-ph_FP)*(1-Y);
    Z_mu = U.*(1-P_neg)+W.*P_neg;
        
    fg = zeros(2,1);
    fg(1) = mean((Y-(1-Y)).*(1-P_neg)./Z_mu);
    fg(2) = mean((Y-(1-Y)).*P_neg./Z_mu);
        
    tmp = (Y-(1-Y)).^2;
    G = zeros(2,2);
    G(1,1) = -mean((tmp.*(1-P_neg).^2)./Z_mu);
    G(2,1) = -mean((tmp.*(1-P_neg).*P_neg)./Z_mu);
    G(1,2) = G(2,1);
    G(2,2) = -mean((tmp.*P_neg.^2)./Z_mu);

    TPFP_new = TPFP-dmp*rand(1)*inv(G)*fg;
    ph_TP = TPFP_new(1);
    ph_FP = TPFP_new(2);

    shusoku = (abs(ph_TP_old-ph_TP)<THETA_M) && (abs(ph_FP_old-ph_FP)<THETA_M);
    
    param.ph_TP = ph_TP;
    param.ph_FP = ph_FP;
    param.rhoh = rhoh;
        
    if(count_M>M_STEP_MAX || sum(isnan(G(:))) > 1)
        conv_check_uv = 0;
        break;
    end

end
    
