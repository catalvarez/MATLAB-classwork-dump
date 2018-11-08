function f_lnL = get_log_likelihood(rho0,K0,K_prime,P,rho,sigma)
% takes parameters and data and solves for the natural log of the
% likelihood of the data given the model

f_lnL = zeros(length(K_prime),1);
for j = 1:length(K_prime)
    rho_m = get_model(K_prime(j),K0(j),rho0(j),P);
    f_lnL(j,1) = -0.5.*sum((rho-rho_m).^2./(sigma.^2));
end

end