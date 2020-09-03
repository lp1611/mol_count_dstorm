function value = nchoosek_ln(a,b)
value = exp(log(factorial(a)) - log(factorial(b))-log(factorial(a-b))); 
end 