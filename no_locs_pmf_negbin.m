function [p, expected, variance] = no_locs_pmf_negbin(lambda,mu,delta,N_F,Delta)

t = N_F*Delta; 
k0 = lambda(1)+mu(1); 
k1 = lambda(2)+mu(2); 

Ploc = exp(-k1*delta); %was 1 before. 
% Ploc = 1; 
ksw = k0*k1/(k0+k1); 
kbl = (lambda(1)*mu(2) + lambda(2)*mu(1) + mu(1)*mu(2))/(k0+k1);  
b = kbl/ksw; %This is b
kfrac = b/(Ploc + b*(1-Ploc));
r = Ploc*ksw*t;
%kfrac = lambda(1)*mu(end)/(lambda(2)+mu(end)); 

p = 0:N_F; 
p(1) = exp(-r); 
for i=1:N_F
    p(i+1) = ((1-kfrac)^i * poisspdf(i,r)) + ... %(i*log(1-kfrac) + i*log(r) + t*ksw - logfactorial(i))
        (kfrac*(1-kfrac)^(i-1) * (poisscdf(i-1,r,'upper'))); 
end 
expected = (0:N_F).*p * ones(N_F+1,1);
expected2 = ((0:N_F).^2).*p*ones(N_F+1,1); 
variance = expected2 - expected^2; 

end 
