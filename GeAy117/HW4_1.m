P = 124000;
x = sin(2.*pi./P.*t);
alpha = sum(2./sigma.^2.*x.^2);
beta = sum(2./sigma.^2);
gamma = sum(2./sigma.^2.*x);
p = sum(2./sigma.^2.*x.*Y);
q = sum(2./sigma.^2.*Y);
N = length(Y);

m0 = (beta.*p-gamma.*q)/(alpha.*beta-gamma.^2);
c0 = (alpha.*q-gamma.*p)/(alpha.*beta-gamma.^2);

covar = 2./(alpha.*beta-gamma.^2).*[beta -gamma; -gamma alpha];
msigma = sqrt(covar(1,1));
csigma = sqrt(covar(2,2));

tline = 0:3720:372000;
line = m0.*sin(2.*pi./P.*tline)+c0;
errorbar(t,Y,sigma,'o')
hold on
plot(tline,line,'g')
title('Method One: Analytical Best Fit')
xlabel('time (years)')
ylabel('obliquity (degrees)')
