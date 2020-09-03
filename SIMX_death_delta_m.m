function [X, X1,states,Y] = SIMX_death_delta_m(lambda2,mu2,Delta,delta,lx,nu2,FPR)
%INPUTS 
%SIMULATE A REALISATION FROM THE PSHMM for m=0,1,2. 
%lambda2 - a 1 by 2*(m+1) vector of transition rates in the order lambda_001, lambda_01, lambda_0102, lambda_011, ..., lambda_0m-10m, lambda_0m-11, lambda_0m1.
%mu2 - a 1 by m+2 vector of absorption/bleached rates in the order mu_0, mu_01, .., mu_0m, mu_1. 
%Delta - 1/Number of frames sampled per second.
%delta - the noise parameter (0<= delta<= Delta) reliable for false negative observations and setting the lower threshold of photon counts (in time) used to determine whether a localisation has been produced in a frame. 
%lx - Total number of frames the experiment has run for. 
%nu2 - the initial probability mass of hidden states in the order nu =
%[nu_0, nu_01, ..., nu_0m, nu_1 nu_2] where nu_2 is usually set to 0. 
%FPR - false positive rate (0<FPR<1). 

%OUTPUTS 
%X - a 1 by lx realisation from the PSHMM with input paramaters for a
%single flurorophore. 

%%GNERATOR MATRIX
m = length(mu2); %m states without the bleached state 
mu2 = [mu2(1) mu2(end) mu2(2:end-1)]; 
nu2 = [nu2(1) nu2(end-1) nu2(2:end-2) nu2(end)]; 

G = zeros(m,m); 
if m==2 
    G(1,2) = lambda2(1); 
    G(2,1) = lambda2(2); 
else 
G(1,2) = lambda2(2); %l_001
G(1,3) = lambda2(1); 
G(2,1) = lambda2(end); %l_10
G(end,2) = lambda2(end-1); %l_(m+1)1
lambda2 = lambda2(3:end-2); 
lam_zerom = lambda2(1:2:end);
lam_ones = lambda2(2:2:end);
end 

if m>3
    for i=1:(m-3) 
        G(2+i,2) = lam_ones(i);
        G(2+i,i+2+1) = lam_zerom(i);
    end 
end 

for i=1:m 
    G(i,i) = -(sum(G(i,:))+mu2(i));
end 

m0 = m-2; %There are m0=m-2 more dark states 
G1 = G; 
G1(m+1,m+1) = 0;
G1(:,m+1) = [mu2 0]';
J = G1; 

for i=1:m
    J(i,:) = -G1(i,:)./G1(i,i);
    J(i,i) = 0; 
end 
J(m+1,m+1) = 1; %%JUMP CHAIN MATRIX

tstop = Delta*lx; 
%nu = [lambda(2),lambda(1)]./(sum(lambda)-lambda(3)); %INITIAL PROBABILITY 
Y(1) = datasample([0,1,2:m0+1,m0+2],1,'Weights',nu2); %initialising the system 
hold = -diag(G); %rate of holding time in state 0,1,0_1,..0_m0 resp

%m0+3 is now the bleached state. 

T(1) = 0; 
n=1; 
while T(n) <= tstop
    if Y(n)~=(m0+2)
        time = -log(rand(1))/hold(Y(n)+1);
        Y(n+1) = datasample([0,1,2:m0+1,m0+2],1,'Weights',J(Y(n)+1,:));
        T(n+1) = T(n) + time;
    else 
        T(n+1)=tstop+1;
    end 
n=n+1;
end 

Y1 = Y; 
if Y(end)~=m0+2 
    Y=Y(1:n-1); 
else 
    Y(end)=0;
end 
%Y1 = Y; 
Y(Y~=1)=0; 

X1 = discritiseTransitionTimes(T, Y, Delta, lx);
X = X1; 
X(X1<=delta)=0; 
X(X1>delta)=1; 
X(X==0) = datasample([0,1],length(X(X==0)),'Weights',[1-FPR FPR]);

states = 1:lx; 
% states(1) = Y1(1); 

T = T(2:end); 
for i=1:lx
    %time i*Delta 
    p = i*Delta < T; 
    f_p = find(p); 
    if sum(p == 0) == 0 
        states(i) = states(1); 
    else
        states(i) = Y1(f_p(1)); 
    end 
end 

states = [Y1(1) states]; 
end 







    



