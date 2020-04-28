function [prob_infect, conv_check] = GT_BP(param)

    % "theta_t" corresponds to theta_tilde in our paper. 

    BP_STEP_MAX = param.BP_STEP_MAX;
    dmp = param.dmp;
    THETA = param.THETA;
    conv_check = 1;
        
    Y = param.Y;
    F = param.F;
    
    N = size(F,2);
    M = size(Y,1);
    
    rhoh = param.rhoh;
    
    p_TP = param.p_TP;
    p_FP = param.p_FP;
    
    theta = rand(N,M).*F';
    theta_t = rand(M,N).*F;

    shusoku = 0;
    bp_step = 0;
    
    while shusoku == 0
    
        bp_step = bp_step + 1;
        
        theta_old = theta;
        theta_t_old = theta_t;

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
    
        % convergence check
        shusoku = (sum(abs(theta_t(:)-theta_t_old(:))<THETA)==M*N);
    
        theta_t = dmp*(theta_t-theta_t_old)+theta_t_old;

        thetat_for_prod = theta_t + (1-F);
        prod_thetat = prod(thetat_for_prod);
        prod_thetat_mat = (prod_thetat'*ones(1,M))./thetat_for_prod';

        neg_mat_h = (1-theta_t);
        thetat_neg = prod(neg_mat_h);
        thetat_neg_mat = (thetat_neg'*ones(1,M))./(1-theta_t');

        Z_p = rhoh*prod_thetat_mat;
        Z_n = (1-rhoh)*thetat_neg_mat;
        theta = F'.*(Z_p./(Z_p+Z_n));
    
        % convergence check
        shusoku = shusoku & (sum(abs(theta(:)-theta_old(:))<THETA)==M*N);
    
        theta = dmp*(theta-theta_old)+theta;
    
        if bp_step > BP_STEP_MAX
            conv_check = 0;
            break;
        end
    end
    
    thetat_for_prod = theta_t + (1-F);
    prod_thetat = prod(thetat_for_prod);

    neg_mat_h = (1-theta_t);
    thetat_neg = prod(neg_mat_h);

    Z_p = rhoh*prod_thetat;
    Z_n = (1-rhoh)*thetat_neg;
    prob_infect = Z_p./(Z_p+Z_n);
    
end

