clear variables;
close all;

parameters_boot;

% generate pooling matrix
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
param.p_TP = p_TP;
param.p_FP = p_FP;
param.rhoh = rhoh;
param.dmp = dmp;
param.THETA = THETA;
param.BP_STEP_MAX = BP_STEP_MAX;

% calculate infection probability by BP
[prob_infect, conv_check] = GT_BP(param);

X_MAP = (prob_infect > 0.5);
TP_MAP = mean((X>0).*X_MAP')/mean(X);
FP_MAP = mean((X==0).*X_MAP')/mean(1-X);

T = log(prob_infect./(1-prob_infect));
    
% start bootstrap
prob_infect_b = zeros(N,boot_sample);
conv_check_b = zeros(boot_sample,1);
b_size = conv_check_b;
Xh_b = zeros(N,boot_sample);
phat = zeros(N,2);

for b = 1:boot_sample

    boot_sample = randi(M,M,1);
    boot_sample = sort(boot_sample);

    for mu = 1:M
        tmp = find(boot_sample == mu);
        if size(tmp,1) > 1
            tmp(1) = [];
            boot_sample(tmp) = [];
        end
    end
    b_size(b) = size(boot_sample,1);
    F_b = F(boot_sample,:);
    Y_b = Y(boot_sample);

    param.Y = Y_b;
    param.F = F_b;

    [prob_infect_b(:,b), conv_check_b(b)] = GT_BP(param);
end
T_b = log(prob_infect_b./(1-prob_infect_b));
sigma_b = sqrt(var(T_b,[],2));
 
X_left = T'-sigma_b*1.96;
X_right = T'+sigma_b*1.96;
    
X_boot = (X_right>0);
    
TP_boot = mean((X>0).*X_boot)/mean(X);
FP_boot = mean((X==0).*X_boot)/mean(1-X);
