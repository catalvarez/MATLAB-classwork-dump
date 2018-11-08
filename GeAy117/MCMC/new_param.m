function trial = new_param(param,stepsize,index)
% generates a new parameter set 
 % draw normally distributed random number for the perturbed parameter
    trial = param;
    trial(index) = randn(1).*stepsize(index)+param(index);
end