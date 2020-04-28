function [X] = get_random_patient(N,K)

    X = zeros(N,1);
    tmp = randperm(N);
    infect = tmp(1:K);
    X(infect) = 1;
    
end

