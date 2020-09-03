function [nu_X,lambda,mu, delta, alpha, B0,B1] = PSHMM_model_select(Y,Delta)

nus = cell(16,3);
lambdas = cell(16,3);
mus = cell(16,3); 
deltas = zeros(16,3);
alphas = zeros(16,3); 
B0s = cell(16,3); 
B1s = cell(16,3); 
BICs = zeros(16,3); 

%First try an m=0 model 
%3 STATE 
lambdas_start = ExampleFitDwellTimes_mstate(Y, Delta,0,1);

%No photo-bleaching 
%zhat = [lambda_01 lambda_10 mu_0 mu_01 delta alpha]
sprintf("Testing model M^0_{NULL}")
zhat = [lambdas_start(1:2) Delta/10 1e-10 0.2];
lb = [0 0 0 0 0];
ub = [Inf Inf Delta 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:2) 0 0 z(3:4)],Delta,[z(5) 1-z(5) 0],0); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_3_null, fval_3_null] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{1,1} = [z_3_null(end) 1-sum(z_3_null(end)) 0]; 
[BP0, BP1] = emission_delta_endstate_3state(z_3_null(1:2), [0 0], z_3_null(3), Delta);
B0s{1,1} = (1-z_3_null(4)).*BP0; 
B1s{1,1} = BP1 + z_3_null(4)*BP0;
lambdas{1,1} = z_3_null(1:2); 
mus{1,1} = [0 0]; 
deltas(1,1) = z_3_null(3);
alphas(1,1) = z_3_null(4); 
BICs(1,1) = 2*fval_3_null + length(z_3_null)*log(numel(Y));

%Photo-bleaching from the Off state 0
sprintf("Testing model M^0_{0}")
zhat = [lambdas{1,1} lambdas_start(end) Delta/10 1e-10 0.2];
lb = [0 0 0 0 0 0];
ub = [Inf Inf Inf Delta 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:2) z(3) 0 z(4:5)],Delta,[z(6) 1-z(6) 0],0); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_3_0, fval_3_0] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{2,1} = [z_3_0(end) 1-sum(z_3_0(end)) 0]; 
[BP0, BP1] = emission_delta_endstate_3state(z_3_0(1:2), [z_3_0(3) 0], z_3_0(4), Delta);
B0s{2,1} = (1-z_3_0(5)).*BP0; 
B1s{2,1} = BP1 + z_3_0(5)*BP0;
lambdas{2,1} = z_3_0(1:2); 
mus{2,1} = [z_3_0(3) 0]; 
deltas(2,1) = z_3_0(4);
alphas(2,1) = z_3_0(5); 
BICs(2,1) = 2*fval_3_0 + length(z_3_0)*log(numel(Y));

%Photo-bleaching from the On state 1
sprintf("Testing model M^0_{1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:2) 0 z(3) z(4:5)],Delta,[z(6) 1-z(6) 0],0); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_3_1, fval_3_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{3,1} = [z_3_1(end) 1-sum(z_3_1(end)) 0]; 
[BP0, BP1] = emission_delta_endstate_3state(z_3_1(1:2),[0 z_3_1(3)], z_3_1(4), Delta);
B0s{3,1} = (1-z_3_1(5)).*BP0; 
B1s{3,1} = BP1 + z_3_1(5)*BP0;
lambdas{3,1} = z_3_1(1:2); 
mus{3,1} = [0 z_3_1(3)]; 
deltas(3,1) = z_3_1(4);
alphas(3,1) = z_3_1(5);
BICs(3,1) = 2*fval_3_1 + length(z_3_1)*log(numel(Y));

%Photo-bleaching from both On and Off states 0,1
sprintf("Testing model M^0_{0,1}")
zhat = [lambdas{1,1} lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2];
lb = [0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Delta 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,z(1:6),Delta,[z(7) 1-z(7) 0],0); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_3_0_1, fval_3_0_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{4,1} = [z_3_0_1(end) 1-sum(z_3_0_1(end)) 0]; 
[BP0, BP1] = emission_delta_endstate_3state(z_3_0_1(1:2), z_3_0_1(3:4), z_3_0_1(5), Delta);
B0s{4,1} = (1-z_3_0_1(6)).*BP0; 
B1s{4,1} = BP1 + z_3_0_1(6)*BP0;
lambdas{4,1} = z_3_0_1(1:2); 
mus{4,1} = z_3_0_1(3:4); 
deltas(4,1) = z_3_0_1(5);
alphas(4,1) = z_3_0_1(6);
BICs(4,1) = 2*fval_3_0_1 + length(z_3_0_1)*log(numel(Y));

%4 STATE 
lambdas_start = ExampleFitDwellTimes_mstate(Y, Delta,1,1);
%zhat = [lambda_001 lambda_01 lambda_011 lambda_10 mu_0 mu_01 mu_1 delta
%alpha]

%Photo-bleaching from the Off state 0
sprintf("Testing model M^1_{NULL}")
zhat = [lambdas{1,1}(1) lambdas_start(2:3) lambdas{1,1}(end) Delta/10 1e-10 0.2 0.1];
lb = [0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Delta 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) 0 0 0 z(5:6)],Delta,[z(7:8) 1-z(7)-z(8) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_null, fval_4_null] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{1,2} = [z_4_null(end-1) z_4_null(end) 1-sum(z_4_null(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_null(1:4), [0 0 0], z_4_null(5), Delta);
B0s{1,2} = (1-z_4_null(6)).*BP0; 
B1s{1,2} = BP1 + z_4_null(6)*BP0;
lambdas{1,2} = z_4_null(1:4); 
mus{1,2} = [0 0 0]; 
deltas(1,2) = z_4_null(5);
alphas(1,2) = z_4_null(6);
BICs(1,2) = 2*fval_4_null + length(z_4_null)*log(numel(Y));

%Photo-bleaching from the first Off state 0
sprintf("Testing model M^1_{0}")
zhat = [lambdas{1,2}/2 lambdas_start(end) Delta/10 1e-10 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Delta 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) z(5) 0 0 z(6:7)],Delta,[z(8:9) 1-z(8)-z(9) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_0, fval_4_0] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{2,2} = [z_4_0(end-1) z_4_0(end) 1-sum(z_4_0(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_0(1:4), [z_4_0(5) 0 0], z_4_0(6), Delta);
B0s{2,2} = (1-z_4_0(7)).*BP0; 
B1s{2,2} = BP1 + z_4_0(7)*BP0;
lambdas{2,2} = z_4_0(1:4); 
mus{2,2} = [z_4_0(5) 0 0]; 
deltas(2,2) = z_4_0(6);
alphas(2,2) = z_4_0(7);
BICs(2,2) = 2*fval_4_0 + length(z_4_0)*log(numel(Y));

%Photo-bleaching from the second Off state 0_1
sprintf("Testing model M^1_{0_1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) 0 z(5) 0 z(6:7)],Delta,[z(8:9) 1-z(8)-z(9) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_01, fval_4_01] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{3,2} = [z_4_01(end-1) z_4_01(end) 1-sum(z_4_01(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_01(1:4), [0 z_4_01(5) 0], z_4_01(6), Delta);
B0s{3,2} = (1-z_4_01(7)).*BP0; 
B1s{3,2} = BP1 + z_4_01(7)*BP0;
lambdas{3,2} = z_4_01(1:4); 
mus{3,2} = [0 z_4_01(5) 0]; 
deltas(3,2) = z_4_01(6);
alphas(3,2) = z_4_01(7);
BICs(3,2) = 2*fval_4_01 + length(z_4_01)*log(numel(Y));

%Photo-bleaching from the On state 1
sprintf("Testing model M^1_{1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) 0 0 z(5) z(6:7)],Delta,[z(8:9) 1-z(8)-z(9) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_1, fval_4_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{4,2} = [z_4_1(end-1) z_4_1(end) 1-sum(z_4_1(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_1(1:4), [0 0 z_4_1(5)], z_4_1(6), Delta);
B0s{4,2} = (1-z_4_1(7)).*BP0; 
B1s{4,2} = BP1 + z_4_1(7)*BP0;
lambdas{4,2} = z_4_1(1:4); 
mus{4,2} = [0 0 z_4_1(5)]; 
deltas(4,2) = z_4_1(6);
alphas(4,2) = z_4_1(7);
BICs(4,2) = 2*fval_4_1 + length(z_4_1)*log(numel(Y));

%Photo-bleaching from the Off states 0 and 0_1
sprintf("Testing model M^1_{0,0_1}")
zhat = [lambdas{1,2}/2 lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2 0.1];
lb = [0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Delta 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) z(5) z(6) 0 z(7:8)],Delta,[z(9:10) 1-z(9)-z(10) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_0_01, fval_4_0_01] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{5,2} = [z_4_0_01(end-1) z_4_0_01(end) 1-sum(z_4_0_01(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_0_01(1:4), [z_4_0_01(5:6) 0], z_4_0_01(7), Delta);
B0s{5,2} = (1-z_4_0_01(8)).*BP0; 
B1s{5,2} = BP1 + z_4_0_01(8)*BP0;
lambdas{5,2} = z_4_0_01(1:4); 
mus{5,2} = [z_4_0_01(5:6) 0]; 
deltas(5,2) = z_4_0_01(7);
alphas(5,2) = z_4_0_01(8);
BICs(5,2) = 2*fval_4_0_01 + length(z_4_0_01)*log(numel(Y));

%Photo-bleaching from the Off and On states 0 and 1
sprintf("Testing model M^1_{0,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) z(5) 0 z(6) z(7:8)],Delta,[z(9:10) 1-z(9)-z(10) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_0_1, fval_4_0_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{6,2} = [z_4_0_1(end-1) z_4_0_1(end) 1-sum(z_4_0_1(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_0_1(1:4), [z_4_0_1(5) 0 z_4_0_1(6)], z_4_0_1(7), Delta);
B0s{6,2} = (1-z_4_0_1(8)).*BP0; 
B1s{6,2} = BP1 + z_4_0_1(8)*BP0;
lambdas{6,2} = z_4_0_1(1:4); 
mus{6,2} = [z_4_0_1(5) 0 z_4_0_1(6)]; 
deltas(6,2) = z_4_0_1(7);
alphas(6,2) = z_4_0_1(8);
BICs(6,2) = 2*fval_4_0_1 + length(z_4_0_1)*log(numel(Y));

%Photo-bleaching from the Off and On states 0_1 and 1
sprintf("Testing model M^1_{0_1,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:4) 0 z(5:6) z(7:8)],Delta,[z(9:10) 1-z(9)-z(10) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_01_1, fval_4_01_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{7,2} = [z_4_01_1(end-1) z_4_01_1(end) 1-sum(z_4_01_1(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_01_1(1:4), [0 z_4_01_1(5:6)], z_4_01_1(7), Delta);
B0s{7,2} = (1-z_4_01_1(8)).*BP0; 
B1s{7,2} = BP1 + z_4_01_1(8)*BP0;
lambdas{7,2} = z_4_01_1(1:4); 
mus{7,2} = [0 z_4_01_1(5:6)]; 
deltas(7,2) = z_4_01_1(7);
alphas(7,2) = z_4_01_1(8);
BICs(7,2) = 2*fval_4_01_1 + length(z_4_01_1)*log(numel(Y));

%Photo-bleaching from all Off and On states 0, 0_1 and 1
sprintf("Testing model M^1_{0,0_1,1}")
zhat = [lambdas{1,2}/2 lambdas_start(end) lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Inf Delta/2 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,z(1:9),Delta,[z(10:11) 1-z(10)-z(11) 0],1); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_4_0_01_1, fval_4_0_01_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{8,2} = [z_4_0_01_1(end-1) z_4_0_01_1(end) 1-sum(z_4_0_01_1(end-1:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_4state(z_4_0_01_1(1:4), z_4_0_01_1(5:7), z_4_0_01_1(8), Delta);
B0s{8,2} = (1-z_4_0_01_1(9)).*BP0; 
B1s{8,2} = BP1 + z_4_0_01_1(9)*BP0;
lambdas{8,2} = z_4_0_01_1(1:4); 
mus{8,2} = z_4_0_01_1(5:7); 
deltas(8,2) = z_4_0_01_1(8);
alphas(8,2) = z_4_0_01_1(9);
BICs(8,2) = 2*fval_4_0_01_1 + length(z_4_0_01_1)*log(numel(Y));

%5 STATE 
lambdas_start = ExampleFitDwellTimes_mstate(Y, Delta,2,1);
%zhat = [lambda_001 lambda_01 lambda_0102 lambda_011 lambda_021 lambda_10
%mu_0 mu_01 mu_1 delta
%alpha]

%No photo-bleaching 
sprintf("Testing model M^2_{NULL}")
zhat = [lambdas{1,2}(1:2) lambdas_start(3) lambdas{1,2}(3) lambdas_start(5) lambdas{1,2}(4) Delta/10 1e-10 0.2 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Delta 1 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 0 0 0 z(7:8)],Delta,[z(9:11) 1-sum(z(9:11)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_null, fval_5_null] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{1,3} = [z_5_null(end-2:end) 1-sum(z_5_null(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_null(1:6), [0 0 0 0], z_5_null(7), Delta);
B0s{1,3} = (1-z_5_null(8)).*BP0; 
B1s{1,3} = BP1 + z_5_null(8)*BP0;
lambdas{1,3} = z_5_null(1:6); 
mus{1,3} = [0 0 0 0]; 
deltas(1,3) = z_5_null(7);
alphas(1,3) = z_5_null(8);
BICs(1,3) = 2*fval_5_null + length(z_5_null)*log(numel(Y));

%Photo-bleaching from the first Off state 0 
sprintf("Testing model M^2_{0}")
zhat = [lambdas{1,3} lambdas_start(end) Delta/10 1e-10 0.2 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Inf Delta 1 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) 0 0 0 z(8:9)],Delta,[z(10:12) 1-sum(z(10:12)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0, fval_5_0] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{2,3} = [z_5_0(end-2:end) 1-sum(z_5_0(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0(1:6), [z_5_0(7) 0 0 0], z_5_0(8), Delta);
B0s{2,3} = (1-z_5_0(9)).*BP0; 
B1s{2,3} = BP1 + z_5_0(9)*BP0;
lambdas{2,3} = z_5_0(1:6); 
mus{2,3} = [z_5_0(7) 0 0 0]; 
deltas(2,3) = z_5_0(8);
alphas(2,3) = z_5_0(9);
BICs(2,3) = 2*fval_5_0 + length(z_5_0)*log(numel(Y));

%Photo-bleaching from the second Off state 0_1 
sprintf("Testing model M^2_{0_1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 z(7) 0 0 z(8:9)],Delta,[z(10:12) 1-sum(z(10:12)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_01, fval_5_01] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{3,3} = [z_5_01(end-2:end) 1-sum(z_5_01(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_01(1:6), [0 z_5_01(7) 0 0], z_5_01(8), Delta);
B0s{3,3} = (1-z_5_01(9)).*BP0; 
B1s{3,3} = BP1 + z_5_01(9)*BP0;
lambdas{3,3} = z_5_01(1:6); 
mus{3,3} = [0 z_5_01(7) 0 0]; 
deltas(3,3) = z_5_01(8);
alphas(3,3) = z_5_01(9);
BICs(3,3) = 2*fval_5_01 + length(z_5_01)*log(numel(Y));

%Photo-bleaching from the third Off state 0_2 
sprintf("Testing model M^2_{0_2}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 0 z(7) 0 z(8:9)],Delta,[z(10:12) 1-sum(z(10:12)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_02, fval_5_02] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{4,3} = [z_5_02(end-2:end) 1-sum(z_5_02(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_02(1:6), [0 0 z_5_02(7) 0], z_5_02(8), Delta);
B0s{4,3} = (1-z_5_02(9)).*BP0; 
B1s{4,3} = BP1 + z_5_02(9)*BP0;
lambdas{4,3} = z_5_02(1:6); 
mus{4,3} = [0 0 z_5_02(7) 0]; 
deltas(4,3) = z_5_02(8);
alphas(4,3) = z_5_02(9);
BICs(4,3) = 2*fval_5_02 + length(z_5_02)*log(numel(Y));

%Photo-bleaching from the On state 1
sprintf("Testing model M^2_{1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 0 0 z(7) z(8:9)],Delta,[z(10:12) 1-sum(z(10:12)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_1, fval_5_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{5,3} = [z_5_1(end-2:end) 1-sum(z_5_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_1(1:6), [0 0 0 z_5_1(7)], z_5_1(8), Delta);
B0s{5,3} = (1-z_5_1(9)).*BP0; 
B1s{5,3} = BP1 + z_5_1(9)*BP0;
lambdas{5,3} = z_5_1(1:6); 
mus{5,3} = [0 0 0 z_5_1(7)]; 
deltas(5,3) = z_5_1(8);
alphas(5,3) = z_5_1(9);
BICs(5,3) = 2*fval_5_1 + length(z_5_1)*log(numel(Y));

%Photo-bleaching from the Off states 0,0_1
sprintf("Testing model M^2_{0,0_1}")
zhat = [lambdas{1,3} lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Inf Inf Delta 1 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) z(8) 0 0 z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_01, fval_5_0_01] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{6,3} = [z_5_0_01(end-2:end) 1-sum(z_5_0_01(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_01(1:6), [z_5_0_01(7:8) 0 0], z_5_0_01(9), Delta);
B0s{6,3} = (1-z_5_0_01(10)).*BP0; 
B1s{6,3} = BP1 + z_5_0_01(10)*BP0;
lambdas{6,3} = z_5_0_01(1:6); 
mus{6,3} = [z_5_0_01(7:8) 0 0]; 
deltas(6,3) = z_5_0_01(9);
alphas(6,3) = z_5_0_01(10);
BICs(6,3) = 2*fval_5_0_01 + length(z_5_0_01)*log(numel(Y));

%Photo-bleaching from the Off states 0, 0_2
sprintf("Testing model M^2_{0,0_2}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) 0 z(8) 0 z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_02, fval_5_0_02] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{7,3} = [z_5_0_02(end-2:end) 1-sum(z_5_0_02(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_02(1:6), [z_5_0_02(7) 0 z_5_0_02(8) 0], z_5_0_02(9), Delta);
B0s{7,3} = (1-z_5_0_02(10)).*BP0; 
B1s{7,3} = BP1 + z_5_0_02(10)*BP0;
lambdas{7,3} = z_5_0_02(1:6); 
mus{7,3} = [z_5_0_02(7) 0 z_5_0_02(8) 0]; 
deltas(7,3) = z_5_0_02(9);
alphas(7,3) = z_5_0_02(10);
BICs(7,3) = 2*fval_5_0_02 + length(z_5_0_02)*log(numel(Y));

%Photo-bleaching from the Off and On states 0,1
sprintf("Testing model M^2_{0,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) 0 0 z(8) z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_1, fval_5_0_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{8,3} = [z_5_0_1(end-2:end) 1-sum(z_5_0_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_1(1:6), [z_5_0_1(7) 0 0 z_5_0_1(8)], z_5_0_1(9), Delta);
B0s{8,3} = (1-z_5_0_1(10)).*BP0; 
B1s{8,3} = BP1 + z_5_0_1(10)*BP0;
lambdas{8,3} = z_5_0_1(1:6); 
mus{8,3} = [z_5_0_1(7) 0 0 z_5_0_1(8)]; 
deltas(8,3) = z_5_0_1(9);
alphas(8,3) = z_5_0_1(10);
BICs(8,3) = 2*fval_5_0_1 + length(z_5_0_1)*log(numel(Y));

%Photo-bleaching from the Off states 0_1,0_2
sprintf("Testing model M^2_{0_1,0_2}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 z(7) z(8) 0 z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_01_02, fval_5_01_02] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{9,3} = [z_5_01_02(end-2:end) 1-sum(z_5_01_02(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_01_02(1:6), [0 z_5_01_02(7:8) 0], z_5_01_02(9), Delta);
B0s{9,3} = (1-z_5_01_02(10)).*BP0; 
B1s{9,3} = BP1 + z_5_01_02(10)*BP0;
lambdas{9,3} = z_5_01_02(1:6); 
mus{9,3} = [0 z_5_01_02(7:8) 0]; 
deltas(9,3) = z_5_01_02(9);
alphas(9,3) = z_5_01_02(10);
BICs(9,3) = 2*fval_5_01_02 + length(z_5_01_02)*log(numel(Y));

%Photo-bleaching from the Off and On states 0_1, 1
sprintf("Testing model M^2_{0_1,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 z(7) 0 z(8) z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_01_1, fval_5_01_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{10,3} = [z_5_01_1(end-2:end) 1-sum(z_5_01_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_01_1(1:6), [0 z_5_01_1(7) 0 z_5_01_1(8)], z_5_01_1(9), Delta);
B0s{10,3} = (1-z_5_01_1(10)).*BP0; 
B1s{10,3} = BP1 + z_5_01_1(10)*BP0;
lambdas{10,3} = z_5_01_1(1:6); 
mus{10,3} = [0 z_5_01_1(7) 0 z_5_01_1(8)]; 
deltas(10,3) = z_5_01_1(9);
alphas(10,3) = z_5_01_1(10);
BICs(10,3) = 2*fval_5_01_1 + length(z_5_01_1)*log(numel(Y));

%Photo-bleaching from the Off and On states 0_2, 1
sprintf("Testing model M^2_{0_2,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 0 z(7) z(8) z(9:10)],Delta,[z(11:13) 1-sum(z(11:13)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_02_1, fval_5_02_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{11,3} = [z_5_02_1(end-2:end) 1-sum(z_5_02_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_02_1(1:6), [0 0 z_5_02_1(7:8)], z_5_02_1(9), Delta);
B0s{11,3} = (1-z_5_02_1(10)).*BP0; 
B1s{11,3} = BP1 + z_5_02_1(10)*BP0;
lambdas{11,3} = z_5_02_1(1:6); 
mus{11,3} = [0 0 z_5_02_1(7:8)]; 
deltas(11,3) = z_5_02_1(9);
alphas(11,3) = z_5_02_1(10);
BICs(11,3) = 2*fval_5_02_1 + length(z_5_02_1)*log(numel(Y));

%Photo-bleaching from the Off states 0, 0_1,0_2
sprintf("Testing model M^2_{0,0_1,0_2}")
zhat = [lambdas{1,3} lambdas_start(end) lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Delta 1 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) z(8) z(9) 0 z(10:11)],Delta,[z(12:14) 1-sum(z(12:14)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_01_02, fval_5_0_01_02] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{12,3} = [z_5_0_01_02(end-2:end) 1-sum(z_5_0_01_02(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_01_02(1:6), [z_5_0_01_02(7:9) 0], z_5_0_01_02(10), Delta);
B0s{12,3} = (1-z_5_0_01_02(11)).*BP0; 
B1s{12,3} = BP1 + z_5_0_01_02(11)*BP0;
lambdas{12,3} = z_5_0_01_02(1:6); 
mus{12,3} = [z_5_0_01_02(7:9) 0]; 
deltas(12,3) = z_5_0_01_02(10);
alphas(12,3) = z_5_0_01_02(11);
BICs(12,3) = 2*fval_5_0_01_02 + length(z_5_0_01_02)*log(numel(Y));

%Photo-bleaching from the Off and On states 0,0_1,1
sprintf("Testing model M^2_{0,0_1,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) z(8) 0 z(9) z(10:11)],Delta,[z(12:14) 1-sum(z(12:14)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_01_1, fval_5_0_01_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{13,3} = [z_5_0_01_1(end-2:end) 1-sum(z_5_0_01_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_01_1(1:6), [z_5_0_01_1(7:8) 0 z_5_0_01_1(9)], z_5_0_01_1(10), Delta);
B0s{13,3} = (1-z_5_0_01_1(11)).*BP0; 
B1s{13,3} = BP1 + z_5_0_01_1(11)*BP0;
lambdas{13,3} = z_5_0_01_1(1:6); 
mus{13,3} = [z_5_0_01_1(7:8) 0 z_5_0_01_1(9)]; 
deltas(13,3) = z_5_0_01_1(10);
alphas(13,3) = z_5_0_01_1(11);
BICs(13,3) =  2*fval_5_0_01_1 + length(z_5_0_01_1)*log(numel(Y));

%Photo-bleaching from the Off and On states 0,0_2,1
sprintf("Testing model M^2_{0,0_2,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) z(7) 0 z(8) z(9) z(10:11)],Delta,[z(12:14) 1-sum(z(12:14)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_02_1, fval_5_0_02_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{14,3} = [z_5_0_02_1(end-2:end) 1-sum(z_5_0_02_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_02_1(1:6), [z_5_0_02_1(7) 0 z_5_0_02_1(8:9)], z_5_0_02_1(10), Delta);
B0s{14,3} = (1-z_5_0_02_1(11)).*BP0; 
B1s{14,3} = BP1 + z_5_0_02_1(11)*BP0;
lambdas{14,3} = z_5_0_02_1(1:6); 
mus{14,3} = [z_5_0_02_1(7) 0 z_5_0_02_1(8:9)]; 
deltas(14,3) = z_5_0_02_1(10);
alphas(14,3) = z_5_0_02_1(11);
BICs(14,3) =  2*fval_5_0_02_1 + length(z_5_0_02_1)*log(numel(Y));

%Photo-bleaching from the Off and On states 0_1, 0_2,1
sprintf("Testing model M^2_{0_1,0_2,1}")
myfun = @(z) -eval_loglik_m_wp2(Y,[z(1:6) 0 z(7) z(8) z(9) z(10:11)],Delta,[z(12:14) 1-sum(z(12:14)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_01_02_1, fval_5_01_02_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{15,3} = [z_5_01_02_1(end-2:end) 1-sum(z_5_01_02_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_01_02_1(1:6), [0 z_5_01_02_1(7:9)], z_5_01_02_1(10), Delta);
B0s{15,3} = (1-z_5_01_02_1(11)).*BP0; 
B1s{15,3} = BP1 + z_5_01_02_1(11)*BP0;
lambdas{15,3} = z_5_01_02_1(1:6); 
mus{15,3} = [0 z_5_01_02_1(7:9)]; 
deltas(15,3) = z_5_01_02_1(10);
alphas(15,3) = z_5_01_02_1(11);
BICs(15,3) = 2*fval_5_01_02_1 + length(z_5_01_02_1)*log(numel(Y));

%Photo-bleaching from all Off and On states 0,0_1, 0_2,1
sprintf("Testing model M^2_{0,0_1,0_2,1}")
zhat = [lambdas{1,3} lambdas_start(end) lambdas_start(end) lambdas_start(end) lambdas_start(end) Delta/10 1e-10 0.2 0.1 0.1];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Delta 1 1 1 1];
A = []; b=[]; Aeq=[]; beq=[]; nonlcol = [];
myfun = @(z) -eval_loglik_m_wp2(Y,z(1:12),Delta,[z(13:15) 1-sum(z(13:15)) 0],2); %Negative log-likelihood function to be minimised.
options = optimset('MaxIter',2e4,'MaxFunEvals',2e4,'FinDiffType','central','UseParallel',0,'TolFun',1e-10,'TolX',1e-10);
[z_5_0_01_02_1, fval_5_0_01_02_1] = fmincon(myfun,zhat,A,b,Aeq,beq,lb,ub,nonlcol,options);
nus{16,3} = [z_5_0_01_02_1(end-2:end) 1-sum(z_5_0_01_02_1(end-2:end)) 0]; 
[BP0, BP1] = emission_delta_endstate_5state(z_5_0_01_02_1(1:6), z_5_0_01_02_1(7:10), z_5_0_01_02_1(11), Delta);
B0s{16,3} = (1-z_5_0_01_02_1(12)).*BP0; 
B1s{16,3} = BP1 + z_5_0_01_02_1(12)*BP0;
lambdas{16,3} = z_5_0_01_02_1(1:6); 
mus{16,3} = z_5_0_01_02_1(7:10); 
deltas(16,3) = z_5_0_01_02_1(11);
alphas(16,3) = z_5_0_01_02_1(12);
BICs(16,3) = 2*fval_5_0_01_02_1 + length(z_5_0_01_02_1)*log(numel(Y))

sprintf("Calculating BICs")
BICs(BICs==0) = Inf; 
minimum=min(min(BICs));
[x,y]=find(BICs==minimum);

nu_X = nus{x,y};
lambda = lambdas{x,y};
mu = mus{x,y}; 
delta = deltas(x,y); 
alpha = alphas(x,y); 
B0 = B0s{x,y}; 
B1 = B1s{x,y}; 
sprintf("Outputting selected model")
end 

