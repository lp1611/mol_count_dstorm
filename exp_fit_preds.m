function [rates, fval] = exp_fit_preds(off_times,lambda,m)
%This function predicts transition rates using the exponential fitting log-likelihood, when m=1,2. 
if m==2
    zhat = log([lambda/2 lambda/3 lambda/3 lambda/5 lambda/3]); %transformation to real line 
    myfun = @(z) -exp_fit_likelihood([z(1) z(2) z(3) z(4) z(5)],off_times,2); 
    options = optimset('MaxIter',2e4,'MaxFunEvals',2e4);
    [z fval] = fminsearch(myfun,zhat,options);
    rates = exp(z); %inverse transformation
elseif m==1 
    zhat = log([lambda lambda lambda 0 0]); 
    myfun = @(z) -exp_fit_likelihood([z(1) z(2) z(3) -Inf -Inf],off_times,1); 
    options = optimset('MaxIter',2e4,'MaxFunEvals',2e4);
    [z fval] = fminsearch(myfun,zhat,options);
    rates = exp(z); 
end 
end 
