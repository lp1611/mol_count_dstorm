function [Q_0, Q_1] = emission_delta_endstate_5state(lambda, mu, delta, Delta) 
%This function computes emission matrices when m=2 (3 dark states). It uses the function qij_5state for inverse Laplace transforms and distributions of waiting times. 
%NOTE: We only compute matrices up until n=2 jumps within a time frame, there will be numerical underflow if the number of frames sampled per second is low. 
%m+1 zero states and 1 on state. 
%lambda = (lambda_001 lambda_01 lambda_0102 lambda_011 lambda_0203 lambda_021 ... lambda_0m1 lambda_10
%mu = [mu_0 mu_01 mu_02 ... mu_0m mu_1]
m = length(mu); %m states without the bleached state  

G = zeros(m,m); 
G(1,2) = lambda(1); %l_001
G(1,end) = lambda(2); %l_01
G(end,1) = lambda(end); %l_10

lam_zerom = lambda(3:2:end-1);
lam_ones = [lambda(4:2:end-1) lambda(end-1)];
for i=1:(m-2) 
    G(1+i,end) = lam_ones(i);
    G(1+i,i+2) = lam_zerom(i);
end 


for i=1:m 
    G(i,i) = -(sum(G(i,:))+mu(i));
end 

sigma = -diag(G);
sigma = sigma(1:end-1); 

Q_0 = zeros(4,5);
Q_1 = zeros(4,5);
Q0_0 = zeros(4,4); 
Q0_1 = Q0_0; 

conv = 1; 
n = 0; 
while conv > 1e-3 && n<=4
    %sprintf("Here")
    [Q, endstate] = qij_5state(n,lambda,mu,Delta);
    
    if n==2
        Q = abs(Q); 
    end 

    Q0 = Q(:,[1:end-2 end]); %ends not at one
    Q1 = Q(:,end-1); %ends at 1 
    
    p0 = cdf('Gamma',delta,n,(-G(m,m))^-1)/cdf('Gamma',Delta,n,(-G(m,m))^-1); %starts not at zero 
    p1 = cdf('Gamma',delta,n+1,(-G(m,m))^-1)/cdf('Gamma',Delta,n+1,(-G(m,m))^-1); %starts at one
    
    Q0_0(end,:) = p1*Q0(end,:); 
    Q0_0(1:end-1,:) = p0*Q0(1:end-1,:); 
    Q0_1(end,:) = (1-p1)*Q0(end,:); 
    Q0_1(1:end-1,:) = (1-p0)*Q0(1:end-1,:); 
    
    if n==0 
        Q1_0 = zeros(m,1); 
        Q1_1 = Q1;
    elseif n==1 
        probs_in_state = 1:m-1; 
        for i=1:m-1
            probs_in_state(i) = (exp(-(-G(i,i))*(Delta-delta))-exp(-(-G(i,i))*Delta))/(1-exp(-(-G(i,i))*Delta));
        end        
        Q1_0 = [probs_in_state probs_in_state(1)]'.*Q1;
        Q1_1 = (ones(1,m) - [probs_in_state probs_in_state(1)])'.*Q1;
    elseif n==2
        probs = zeros(m,m-1); 
        probs(1,1) = cdf('Gamma',Delta-delta,2,(-G(1,1))^-1)/cdf('Gamma',Delta,2,(-G(1,1))^-1);
        if abs(G(1,1) - G(2,2)) < 1e-4 %To reduce numerical overflow. 
            probs(1,2) = cdf('Gamma',Delta-delta,2,(-G(1,1))^-1)/cdf('Gamma',Delta,2,(-G(1,1))^-1);
        else %sum of 2 exponentials (closed form). 
            probs(1,2) = 1 - (-G(2,2)*(1-exp(-G(1,1)*(Delta-delta)))+G(1,1)*(1-exp(-G(2,2)*(Delta-delta))))/(-G(2,2)*(1-exp(-G(1,1)*(Delta)))+G(1,1)*(1-exp(-G(2,2)*(Delta))));
        end 
        probs(2,2) = probs(1,2);
        if abs(G(2,2) - G(3,3)) < 1e-4 
            probs(2,3) = cdf('Gamma',Delta-delta,2,(-G(2,2))^-1)/cdf('Gamma',Delta,2,(-G(2,2))^-1);
        else %sum of 2 exponentials (closed form). 
            probs(2,3) = 1 - (-G(2,2)*(1-exp(-G(3,3)*(Delta-delta)))+G(3,3)*(1-exp(-G(2,2)*(Delta-delta))))/(-G(2,2)*(1-exp(-G(3,3)*(Delta)))+G(3,3)*(1-exp(-G(2,2)*(Delta))));
        end

        if abs(G(1,1) - G(3,3)) < 1e-4
            probs(3,3) = cdf('Gamma',Delta-delta,2,(-G(1,1))^-1)/cdf('Gamma',Delta,2,(-G(1,1))^-1);
        else %sum of 2 exponentials (closed form). 
            probs(3,3) = 1 - (-G(3,3)*(1-exp(-G(1,1)*(Delta-delta)))+G(1,1)*(1-exp(-G(3,3)*(Delta-delta))))/(-G(3,3)*(1-exp(-G(1,1)*(Delta)))+G(1,1)*(1-exp(-G(3,3)*(Delta))));
        end 
        probs(4,:) = probs(1,:); %starts at 1

        Q1_1 = sum(probs.*endstate,2);
        Q1_0 = sum((1-probs).*endstate,2);
    elseif n==3 %THIS SHOULD BE CORRECT 
        probs = zeros(m,m-1); 
        probs(1,1) = cdf('Gamma',Delta-delta,3,(-G(1,1))^-1)/cdf('Gamma',Delta,3,(-G(1,1))^-1);
        probs(1,2) = Sum_of_gamma_cdf(sigma, [2 1 0]', Delta-delta)/Sum_of_gamma_cdf(sigma, [2 1 0]', Delta);
        probs(1,3) = Sum_of_gamma_cdf(sigma, [1 1 1]', Delta-delta)/Sum_of_gamma_cdf(sigma, [1 1 1]', Delta);
        probs(2,1) = probs(1,2); 
        probs(2,3) = Sum_of_gamma_cdf(sigma, [1 2 0]', Delta-delta)/Sum_of_gamma_cdf(sigma, [1 2 0]', Delta);
        probs(2,2) = probs(1,3); 
        probs(3,2) = Sum_of_gamma_cdf(sigma, [2 0 1]', Delta-delta)/Sum_of_gamma_cdf(sigma, [2 0 1]', Delta);
        probs(3,3) = probs(1,3); 
        probs(4,:) = probs(1,:); 
        
        Q1_1 = sum((probs).*endstate,2);
        Q1_0 = sum((1-probs).*endstate,2);
    elseif n==4
        probs = zeros(m,m);
        probs(1,1) = cdf('Gamma',Delta-delta,4,(-G(1,1))^-1)/cdf('Gamma',Delta,4,(-G(1,1))^-1);
        probs(1,2) = Sum_of_gamma_cdf(sigma, [3 1 0]', Delta-delta)/Sum_of_gamma_cdf(sigma, [2 1 0]', Delta);
        probs(1,3) = Sum_of_gamma_cdf(sigma, [2 1 1]', Delta-delta)/Sum_of_gamma_cdf(sigma, [2 1 1]', Delta);
        probs(1,4) = Sum_of_gamma_cdf(sigma, [2 2 0]', Delta-delta)/Sum_of_gamma_cdf(sigma, [2 2 0]', Delta);
        probs(2,1) = Sum_of_gamma_cdf(sigma, [1 2 1]', Delta-delta)/Sum_of_gamma_cdf(sigma, [1 2 1]', Delta);
        probs(2,2) = probs(1,2); 
        probs(2,3) = probs(1,3); 
        probs(2,4) = probs(1,4); 
        probs(3,1) = Sum_of_gamma_cdf(sigma, [3 0 1]', Delta-delta)/Sum_of_gamma_cdf(sigma, [3 0 1]', Delta);
        probs(3,2) = Sum_of_gamma_cdf(sigma, [1 1 2]', Delta-delta)/Sum_of_gamma_cdf(sigma, [1 1 2]', Delta);
        probs(3,3) = probs(1,3); 
        probs(4,:) = probs(1,:); 
        
        Q1_1 = sum((probs).*endstate,2);
        Q1_0 = sum((1-probs).*endstate,2);
        
    end 
    
    Q_0_ = Q_0 + [Q0_0(:,1:m-1) Q1_0 Q0_0(:,end)]; 
    Q_1_ = Q_1 + [Q0_1(:,1:m-1) Q1_1 Q0_1(:,end)]; 
        
    conv_0 = norm(Q_0_ - Q_0,'fro');
    conv_1 = norm(Q_1_ - Q_1,'fro');
    conv = max(conv_0, conv_1); 
    Q_0 = Q_0_; 
    Q_1 = Q_1_; 
    n = n+1; 
end 
%sprintf("n is %d", n)

Q_0(end+1,:) = zeros(1,m+1);
Q_0(end,end) = 1; %emission matrix when 0 is observed. 
Q_1(end+1,:) = zeros(1,m+1); %emission matrix when 1 is observed. 

G = [-sum(lambda(1:2))-mu(1) lambda(1) 0 lambda(2) mu(1); 0 -lambda(3)-lambda(4)-mu(2) lambda(3) lambda(4) mu(2); 0 0 -lambda(5)-mu(3) lambda(5) mu(3); lambda(6) 0 0 -lambda(6)-mu(end) mu(end) ; 0 0 0 0 0];
Q_0 = expm(G*Delta)-Q_1;

end 


    





