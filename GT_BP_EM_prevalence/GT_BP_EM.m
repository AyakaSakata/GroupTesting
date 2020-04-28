function [prob_infect, rhoh, P_neg, conv_check] = GT_BP_EM(param)

    BP_STEP_MAX = param.BP_STEP_MAX;
    dmp = param.dmp;
    THETA = param.THETA;
    conv_check = 1;
    
    Y = param.Y;
    F = param.F;
    
    N = size(F,2);
    M = size(Y,1);
    
    rhoh = param.rhoh;
    ph_TP = param.ph_TP;
    ph_FP = param.ph_FP;

    theta = rand(N,M).*F';
    theta_t = rand(M,N).*F;
    
    Y_mat = Y*ones(1,N);
    U_mat = ph_TP*Y_mat+(1-ph_TP)*(1-Y_mat);
    W_mat = ph_FP*Y_mat+(1-ph_FP)*(1-Y_mat);

    shusoku = 0;
    bp_step = 0;
    
    while shusoku == 0
    
        bp_step = bp_step + 1;
        
        theta_old = theta;
        thetat_old = theta_t;
        
        neg_mat = (1-theta);
        neg_vec = prod(neg_mat);
        P_neg = neg_vec'*ones(1,N);
        P_neg_div = 1-theta';
        P_neg = F.*(P_neg./P_neg_div);

        Zt_p = U_mat;
        Zt_n = U_mat.*(1-P_neg)+W_mat.*P_neg;
        Zt = Zt_p+Zt_n;
        theta_t = F.*(U_mat./Zt);
    
        shusoku = (sum(abs(theta_t(:)-thetat_old(:))<THETA)==M*N);
    
        theta_t = dmp*(theta_t-thetat_old)+thetat_old;

        thetat_for_prod = theta_t + (1-F);
        prod_thetat = prod(thetat_for_prod);
        prod_thetat_mat = (prod_thetat'*ones(1,M))./thetat_for_prod';

        neg_mat_t = (1-theta_t);
        Pt_neg = prod(neg_mat_t);
        Pt_neg_mat = (Pt_neg'*ones(1,M))./(1-theta_t');

        Z_p = rhoh*prod_thetat_mat;
        Z_n = (1-rhoh)*Pt_neg_mat;
        theta = F'.*(Z_p./(Z_p+Z_n));
    
        shusoku = shusoku & (sum(abs(theta(:)-theta_old(:))<THETA)==M*N);
    
        theta = dmp*(theta-theta_old)+theta;
        
        neg_mat = 1-theta;
        P_neg = prod(neg_mat)';
        
        if (bp_step > BP_STEP_MAX) || (sum(isnan([ph_TP,ph_FP]))>0)
            conv_check = 0;
            break;
        end
        
    end
    
    % update of rhoh
    thetat_for_prod = theta_t + (1-F);
    prod_thetat = prod(thetat_for_prod);
    
    neg_mat_t = (1-theta_t);
    Pt_neg = prod(neg_mat_t);

    Z_p = rhoh*prod_thetat;
    Z_n = (1-rhoh)*Pt_neg;
    prob_infect = Z_p./(Z_p+Z_n);
    
    rhoh = mean(prob_infect);
    
end

