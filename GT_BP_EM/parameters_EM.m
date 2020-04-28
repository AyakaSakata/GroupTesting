%% parameters used in GT_BP

% damping factor
dmp = 0.1;
% threshold for convergence of E-step
THETA_E = 1e-5;
% threshold for convergence of M-step
THETA_M = 1e-3;
% maximum length of BP step
BP_STEP_MAX = 1e+4;
% maximum length of M-step
M_STEP_MAX = 1e+2;

% number of patients
N = 1000;
% number of pools
M = 500;
% group size
N_G = 10;
% overlap size
N_O = M*N_G/N;

% true positive probability in the test
p_TP = 0.95;
% false positive probability in the test
p_FP = 0.1;

% True prevalence
rho = 0.1;
% Assumed prevalence
rhoh = rho+0.01;
K = ceil(N*rho);
