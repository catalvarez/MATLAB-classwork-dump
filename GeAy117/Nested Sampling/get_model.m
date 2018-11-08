function rho_m = get_model(K_prime,K0,rho0,P)
% Takes parameters and solves the EOS model for density
rho_m = rho0.*(1+K_prime.*P./K0).^(1./K_prime);

end