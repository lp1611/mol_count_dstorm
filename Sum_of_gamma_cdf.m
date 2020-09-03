function p = Sum_of_gamma_cdf(lambda, n, x) 
%lambda is a vector of lambdas (s_0, s_01, s_02, ... s_0m) 
%n is a vector of the number of time pieces in each state... INTEGERS

lambda = 1./lambda; 
lambda = lambda(~~n); 
n = n(~~n); 

% sum = 0; 
% for j=1:m 
%     for r=1:n(j)
%         bjr = b_jr(lambda, n,j,r); % need to do this function  
%         sum = sum + (-lambda(j))^r * bjr * gammainc(x/lambda(j),r,'upper'); 
%     end 
% end 
% 
% p = 1 - sum*prod((-lambda).^n); 

beta = min(lambda); 
C = prod((beta./lambda).^n);
gamma = zeros(1,100+1); 

for k=1:100
    gamma(k) = sum(n.*(1-beta./lambda).^k./k);
end 
rho = sum(n); 

delta = ones(1,100+1);  
for k=1:100
    for i=1:k+1 
    delta(k+1) = delta(k+1) + (i*gamma(i)*delta(k+2-i))/(k+1);
    end 
end 

p(1)=0;
conv = 1; 
i=0; 
while (conv > 1e-6) 
    p(i+2) = p(i+1) + C*delta(i+1)*(beta^-1)*gammainc(x/beta,rho+i,'lower');
    conv = abs(p(i+2)-p(i+1));
    i=i+1; 
end 
p = p(end); 
end 


    
    
