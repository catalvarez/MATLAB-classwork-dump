sigma = 2;
A = 41.4;
mu = 34:0.1:48;

P_A_mu = 1./(sigma.*sqrt(2*pi)).*exp(-.5*((A-mu)./sigma).^2);

sum(P_A_mu)
plot(mu,P_A_mu)
title('The Answer to LUE')
xlabel('Mean value')
ylabel('Probability of A given a mean value')
xlim([35 55])
ylim([0 0.21])
