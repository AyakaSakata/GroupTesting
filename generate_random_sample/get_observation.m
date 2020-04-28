function [Y] = get_observation(Y_0,p_TP,p_FP)

    pos_ind = find(Y_0==1);
    pos_size = size(pos_ind,1);
    neg_ind = find(Y_0==0);
    neg_size = size(neg_ind,1);

    pos_mut = (rand(pos_size,1)<1-p_TP);
    neg_mut = (rand(neg_size,1)<p_FP);

    Y = Y_0;
    Y(pos_ind) = Y_0(pos_ind)-pos_mut;
    Y(neg_ind) = Y_0(neg_ind)+neg_mut;
end

