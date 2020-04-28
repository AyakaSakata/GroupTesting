%% parameters used in GT_BP

% damping factor
dmp = 0.1;
% threshold for convergence
THETA = 1e-5;
% maximum length of BP step
BP_STEP_MAX = 1e+4;

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
rho = 0.05;
% Assumed prevalence
rhoh = rho;
K = ceil(N*rho);

% number of bootstrap samples
boot_sample = 100;

