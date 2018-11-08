sigma1 = 2;
A1 = 41.4;
sigma2 = 3;
A2 = 46.9;
mu = 35:0.1:55;

P_A_mu = 1./(sigma2.*sqrt(2*pi)).*exp(-.5*((A2-mu)./sigma2).^2) ... 
    .*1./(sigma1.*sqrt(2*pi)).*exp(-.5*((A1-mu)./sigma1).^2);

% sum(P_A_mu)
plot(mu,P_A_mu)
title('The Answer to LUE')
xlabel('Mean value')
ylabel('Probability of A given a mean value')
% xlim([35 55])
% ylim([0 0.21])
