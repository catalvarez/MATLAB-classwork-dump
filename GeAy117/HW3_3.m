f = 0:0.01:1;

P = factorial(21)./(factorial(3).*factorial(18)).*(0.6*f).^2.*f.*(1-0.6*f).^8.*(1-f).^10;

plot(f,P)
title('Stars with Binary Companions: Sensitivity')
xlabel('Possible fraction with companions, f')
ylabel('Probability')