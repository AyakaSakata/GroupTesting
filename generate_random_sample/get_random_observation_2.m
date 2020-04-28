function [ F ] = get_random_observation_2(M,N,N_G,N_O)

    max_loop = 10000;
    
    labels_length = N_G*M;
    
    test_labels = [];
    for mu = 1: M
        tmp = mu*ones(1,N_G);
        test_labels = [test_labels,tmp];
    end
    var_labels = [];
    for i = 1:N
        tmp = i*ones(1,N_O);
        var_labels = [var_labels,tmp];
    end

    restart = 1;
    while restart == 1
        restart = 0;
        arrange = randperm(labels_length);
        F = zeros(M,N);
    
        for i = 1:labels_length
            F(test_labels(arrange(i)),var_labels(i)) = 1;
        end

        check_var = (sum((sum(F)==N_O)) == N);
        check_test = (sum((sum(F,2)==N_G)) == M);
        check = check_var && check_test;

        while check == 0
            test_less = find(sum(F,2) < N_G);
            n_test_less = N_G - sum(F(test_less,:),2);
            var_less = find(sum(F) < N_O);
            n_var_less = N_O - sum(F(:,var_less));
        
            test_less_size = size(test_less,1);
            test_less_labels = [];
            for i = 1:test_less_size
                tmp = test_less(i)*ones(1,n_test_less(i));
                test_less_labels = [test_less_labels,tmp];
            end
    
            var_less_size = size(var_less,2);
            var_less_labels = [];
            for i = 1:var_less_size
                tmp = var_less(i)*ones(1,n_var_less(i));
                var_less_labels = [var_less_labels,tmp];
            end
    
            less_length = sum(n_var_less);
    
            double_check = 1;
            loop_count = 0;
            while double_check >= 1
            
                loop_count = loop_count + 1;
                arrange = randperm(less_length);
    
                double_check = 0;
                for i = 1:less_length
                    double_check = double_check + F(test_less_labels(arrange(i)),var_less_labels(i));
                end
                if loop_count > max_loop
                    restart = 1;
                    break;
                end
            end
            
            if restart == 1
                break;
            end
    
            for i = 1:less_length
                F(test_less_labels(arrange(i)),var_less_labels(i)) = 1;
            end
    
            check_var = (sum((sum(F) == N_O)) == N);
            check_test = (sum((sum(F,2) == N_G)) == M);
            check = check_var && check_test;

        end
    end
end

