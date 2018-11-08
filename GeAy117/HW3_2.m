f = 0:0.01:1;

P = factorial(21)./(factorial(3).*factorial(18)).*f.^3.*(1-f).^18;

plot(f,P)
title('Stars with Binary Companions: Basic')
xlabel('Possible fraction with companions, f')
ylabel('Probability')