function [new_param,logL_new] = get_new_sample(A_param,logL,P,rho,sigma)
% Takes parameters, likelihoods, and data and generates a new nested
% sampling point
% All input values should be matrices of all active points

% initial step sizes, from HW5_1
rho0_sigma = 0.0011;
K_prime_sigma = 0.6776;
K0_sigma = 5.1415;
stepsize = [K_prime_sigma; K0_sigma; rho0_sigma];

% variables for number of steps, steps accepted, and steps rejected
nsteps = 40;
nAccepts = 0;
nRejects = 0;

% get lowest likelihood value
[log_Lstar,worst] = min(logL);

% random pick of active point, not the lowest likelihood point
i = ceil(rand(1)*length(A_param(1,:)));
while i == worst
    i = ceil(rand(1)*length(A_param(1,:)));
end

% parameter values from the random active point
new_param = A_param(:,i);

for k = 1:nsteps
    % draw normally distributed random number for each parameter
    d_param = randn(3,1).*stepsize;
    % trial set of parameters, current guess+random number
    try_param = new_param+d_param;
    % compute likelihood of trial parameters
    try_logL = get_log_likelihood(try_param(3,1),try_param(2,1),try_param(1,1),P,rho,sigma);
    % sort trial parameter as either accepted or rejected point
    if try_logL >= log_Lstar
        new_param = try_param;
        nAccepts = nAccepts+1;
    else
        nRejects = nRejects+1;
    end
    % keeps the acceptance ratio near 50%
    if nAccepts >= nRejects
        stepsize = stepsize.*exp(1./nAccepts);
    elseif nRejects >= nAccepts
        stepsize = stepsize./exp(1.0/nRejects);
    end
end

% Likelihood value for new point
logL_new = get_log_likelihood(new_param(3,1),new_param(2,1),new_param(1,1),P,rho,sigma);

end