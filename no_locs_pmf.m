function [probs, exp, variance] = no_locs_pmf(nu,B0,B1,NF) 

probs = analytic_prs_3(nu,B0,B1,NF);

exp = (0:NF).*probs * ones(NF+1,1);
exp2 = ((0:NF).^2).*probs*ones(NF+1,1); 
variance = exp2 - exp^2; 

end 
