qs = 0.07;
F = 6.2*10^-12;
Ts = 0+273;
cp = 10^3;
rho = 3000;
k = 3;
kappa = 10^-6;

z = 0:500;
T = -F/(2*cp*kappa*rho).*z.^2-qs/k.*z+Ts;
plot(T,z)