function [I, p] = hpd_interval_discrete(p_no_locs_blinks_alexa,N_mode,alpha)

i=1; 
p = p_no_locs_blinks_alexa(N_mode);
I = N_mode; 

while p < 1-alpha
    vec = [p_no_locs_blinks_alexa(N_mode-i),p_no_locs_blinks_alexa(N_mode+i)];
    kvec = [N_mode-i,N_mode+i]; 
    [min_int, I_min] = min(vec);
    [max_int, I_max] = max(vec);
    
    p = p + min_int; 
    I = [I kvec(I_min)]; 
    if p >= 1-alpha 
        break;
    else 
        p = p + max_int; 
        I = [I kvec(I_max)]; 
        i = i+1; 
    end 
end 

I = [min(I) max(I)]; 

end 
