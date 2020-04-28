function [prob_infect, rhoh_dist, conv_check] = GT_HBP(param)

    BP_STEP_MAX = param.BP_STEP_MAX;
    dmp = param.dmp;
    THETA = param.THETA;
    conv_check = 1;
    
    Y = param.Y;
    F = param.F;
    
    N = size(F,2);
    M = size(Y,1);
    
    d_rho = 0.001;
    rho_tmp = d_rho:d_rho:1;
    rho_size = size(rho_tmp,2);
    rho_mat = rho_tmp'*ones(1,N);
    
    beta_a = param.beta_a;
    beta_b = param.beta_b;
    
    phi = (rho_tmp.^(beta_a-1)).*((1-rho_tmp).^(beta_b-1))/beta(beta_a,beta_b)*d_rho;
    
    p_TP = param.p_TP;
    p_FP = param.p_FP;
    
    theta = rand(N,M).*F';
    theta_t = rand(M,N).*F;
    
    pi = rand(N,1);
    rho_t = M/N*rand(N,1);

    shusoku = 0;
    bp_step = 0;
    
    while shusoku == 0
    
        bp_step = bp_step + 1;
        
        theta_old = theta;
        thetat_old = theta_t;

        pi_old = pi;
        rhot_old = rho_t;

        rhoh_dist = rho_t;
        rhoh_mat = rhoh_dist*ones(1,M);

        neg_mat = (1-theta);
        neg_vec = prod(neg_mat);
        P_neg = neg_vec'*ones(1,N);
        P_neg_div = 1-theta';
        P_neg = F.*(P_neg./P_neg_div);

        Y_mat = Y*ones(1,N);
        Zt_p = p_TP*Y_mat+(1-p_TP)*(1-Y_mat);
        Zt_n = (p_TP*Y_mat+(1-p_TP)*(1-Y_mat)).*(1-P_neg)...
            +(p_FP*Y_mat+(1-p_FP)*(1-Y_mat)).*P_neg;
        Zt = Zt_p+Zt_n;
        theta_t = F.*((p_TP*Y_mat+(1-p_TP)*(1-Y_mat))./Zt);
    
        shusoku = (sum(abs(theta_t(:)-thetat_old(:))<THETA)==M*N);
    
        theta_t = dmp*(theta_t-thetat_old)+thetat_old;

        thetat_for_prod = theta_t + (1-F);
        prod_thetat = prod(thetat_for_prod);
        prod_thetat_mat = (prod_thetat'*ones(1,M))./thetat_for_prod';

        neg_mat_h = (1-theta_t);
        Ph_neg = prod(neg_mat_h);
        Ph_neg_mat = (Ph_neg'*ones(1,M))./(1-theta_t');

        Z_p = rhoh_mat.*prod_thetat_mat;
        Z_n = (1-rhoh_mat).*Ph_neg_mat;
        theta = F'.*(Z_p./(Z_p+Z_n));
    
        shusoku = shusoku & (sum(abs(theta(:)-theta_old(:))<THETA)==M*N);
    
        theta = dmp*(theta-theta_old)+theta;
        
        pi = (prod_thetat./(prod_thetat+Ph_neg))';

        shusoku = shusoku & (sum(abs(pi-pi_old)<THETA) == N);        
        pi = dmp*(pi-pi_old)+pi;
        
        pi_mat = ones(rho_size,1)*pi';
        pi_den_tmp_mat = 2*(rho_mat.*pi_mat+(1-rho_mat).*(1-pi_mat));
        pi_den_tmp = prod(pi_den_tmp_mat,2);
        pi_den_mat = pi_den_tmp*ones(1,N)./pi_den_tmp_mat;
        pi_den = phi*pi_den_mat;
        
        rho_t = ((phi*(rho_mat.*pi_den_mat))./pi_den)';

        shusoku = shusoku & (sum(abs(rho_t-rhot_old)<THETA) == N);
        rho_t = dmp*(rho_t-rhot_old)+rho_t;

        if bp_step > BP_STEP_MAX
            conv_check = 0;
            break;
        end
    end
    
    rhoh_dist = rho_t;
    
    thetat_for_prod = theta_t + (1-F);
    prod_thetat = prod(thetat_for_prod);

    neg_mat_h = (1-theta_t);
    Ph_neg = prod(neg_mat_h);

    Z_p = rhoh_dist'.*prod_thetat;
    Z_n = (1-rhoh_dist').*Ph_neg;
    prob_infect = Z_p./(Z_p+Z_n);
    
end

