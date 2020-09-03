%ALEXA DATA tested on Robert N's neg binomial dist 
load testing_number_of_locs_alexa_2.mat

testing_number_of_locs_alexa_negbin_2 = testing_number_of_locs_alexa_2; 

for i=1:27 
        Delta = testing_number_of_locs_alexa_2(i).Delta/10; 
        Y_train = testing_number_of_locs_alexa_2(i).Y_train;
        Y_test = testing_number_of_locs_alexa_2(i).Y_test;
        N_F = size(testing_number_of_locs_alexa_2(i).Y_test,2); 
        
        %3 STATE 
        zhat = [[10 50] 20 Delta/10 1e-6 0.05];
        lb = [0 0 0 0 0 0];
        ub = [Inf Inf Inf Delta 1 1];
        A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
        myfun = @(z) -eval_loglik_m_wp2(Y_train,[z(1:2) 0 z(3:5)],Delta,[z(6) 1-sum(z(6)) 0],0); %Negative log-likelihood function to be minimised.
        %myfun = @(z) -no_locs_pmf_negbin_all(Y_train,[z(1:2) 0 z(3:5)],Delta,[z(6) 1-sum(z(6)) 0],0); %Negative log-likelihood function to be minimised.
        options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-6,'TolX',1e-6);
        z = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
        z_mle = z; 

        lambda = z(1:2); 
        mu = [0 z(3)]; 
        delta = z(4); 
        alpha = z(5); 
        raw_params = [lambda mu delta alpha]; 
        nu_X = [z(6) 1-sum(z(6)) 0]; 
%         [BP0, BP1] = emission_delta_endstate_3state(lambda, mu, delta, Delta);
%         B0 = (1-alpha).*BP0; 
%         B1 = BP1 + alpha*BP0; 
        testing_number_of_locs_alexa_negbin_2(i).raw_params = raw_params;
        testing_number_of_locs_alexa_negbin_2(i).nu_X = nu_X;


        n_obs_in_frame_test = sum(Y_test,1);
        testing_number_of_locs_alexa_negbin_2(i).n_obs_in_frame = n_obs_in_frame_test;

        [probs, expected,variance] = no_locs_pmf_negbin(lambda,mu,delta,N_F,Delta); 
        khat = ceil(sum(n_obs_in_frame_test)/expected); 
        kmin = max(ceil((sum(n_obs_in_frame_test))/N_F), floor(khat - 3*sqrt(variance))); 
        kmax = ceil(khat + 3*sqrt(variance));
        p_no_locs = zeros(1,kmax+1); 
        prior = 1:kmax+1; 

        for j=kmin:kmax
            vec_probs = abs(ifft(fft([probs zeros(1,j*N_F+1-length(probs))]).^j));
            probv = vec_probs(sum(n_obs_in_frame_test)+1);
            p_no_locs(j) = probv;
            j
        end 
        clear vec_probs 

        p_no_locs_blinks_alexa = p_no_locs/sum(p_no_locs);
        testing_number_of_locs_alexa_negbin_2(i).p_no_locs_blinks_alexa = p_no_locs_blinks_alexa;
        [N_val,testing_number_of_locs_alexa_negbin_2(i).cum_mode] = max(p_no_locs_blinks_alexa); 
        %testing_number_of_locs_alexa_negbin_2(i).cum_mode = cum_mode; 
        [testing_number_of_locs_alexa_negbin_2(i).cum_CI, testing_number_of_locs_alexa_negbin_2(i).post_coverage] = hpd_interval_discrete(p_no_locs_blinks_alexa,testing_number_of_locs_alexa_negbin_2(i).cum_mode,.05);
%         
%         v = cumsum(p_no_locs_blinks_alexa) < 0.975 & cumsum(p_no_locs_blinks_alexa) > 0.025; %95% of the variance
%         ks = prior(diff(v)~=0);

%         if N_true >= ks(1) && N_true <= ks(2)
%             in_CI = 1; 
            save -v7.3 testing_number_of_locs_alexa_negbin_2.mat testing_number_of_locs_alexa_negbin_2
%         end 
%     end 
end
