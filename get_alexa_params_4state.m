%ALEXA DATA
load testing_number_of_locs_alexa_2.mat
testing_number_of_locs_alexa_4state = testing_number_of_locs_alexa_2; 
clear testing_number_of_locs_alexa_2

for i=2:27 
Delta = testing_number_of_locs_alexa_4state(i).Delta; 
Y_train = testing_number_of_locs_alexa_4state(i).Y_train; 
Y_test = testing_number_of_locs_alexa_4state(i).Y_test; 
N_F = size(Y_train,2);
    
%4 STATE 
%zhat = [lambda_001 lambda_01 lambda_011 lambda_10 mu_0 mu_01 mu_1 delta
%alpha]
zhat = [[10 10 1 30] 20 Delta/10 1e-6 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Delta 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y_train,[z(1:4) 0 0 z(5:7)],Delta, ...
    [z(8:9) 1-sum(z(8:9)) 0],1); 
%Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType',...
    'central','UseParallel',1,'TolFun',1e-8,'TolX',1e-8);
z = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
z_mle = z; 

lambda = z(1:4); 
mu = [0 0 z(5)]; 
delta = z(6); 
alpha = z(7); 
raw_params = [lambda mu delta alpha]; 
nu_X = [z(8:9) 1-sum(z(8:9)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(lambda, mu, delta, Delta);
B0 = (1-alpha).*BP0; 
B1 = BP1 + alpha*BP0; 
testing_number_of_locs_alexa_4state(i).raw_params = raw_params;
testing_number_of_locs_alexa_4state(i).nu_X = nu_X;
testing_number_of_locs_alexa_4state(i).B0 = B0;
testing_number_of_locs_alexa_4state(i).B1 = B1;

n_obs_in_frame_test = sum(Y_test,1);
testing_number_of_locs_alexa_4state(i).n_obs_in_frame = n_obs_in_frame_test;

[probs, expected,variance] = no_locs_pmf(nu_X,B0,B1,N_F); 
khat = ceil(sum(n_obs_in_frame_test)/expected); 
kmin = max(ceil((sum(n_obs_in_frame_test))/N_F), ...
    floor(khat - 3*sqrt(variance))); 
kmax = ceil(khat + 3*sqrt(variance));
p_no_locs = zeros(1,kmax+1); 
prior = 1:kmax+1; 

parfor j=kmin:kmax
    vec_probs = abs(ifft(fft([probs zeros(1,j*N_F+1-length(probs))]).^j));
    probv = vec_probs(sum(n_obs_in_frame_test)+1);
    p_no_locs(j) = probv;
end 
clear vec_probs 

p_no_locs_blinks_alexa = p_no_locs/sum(p_no_locs);
testing_number_of_locs_alexa_4state(i).p_no_locs_blinks_alexa = ...
    p_no_locs_blinks_alexa;
[N_val,testing_number_of_locs_alexa_4state(i).cum_mode] = ...
    max(p_no_locs_blinks_alexa); 
%testing_number_of_locs_alexa_negbin(i).cum_mode = cum_mode; 
[testing_number_of_locs_alexa_4state(i).cum_CI, ...
    testing_number_of_locs_alexa_4state(i).post_coverage] = ...
    hpd_interval_discrete(p_no_locs_blinks_alexa,...
    testing_number_of_locs_alexa_4state(i).cum_mode,.05);
end
save -v7.3 testing_number_of_locs_alexa_4state.mat ...
    testing_number_of_locs_alexa_4state