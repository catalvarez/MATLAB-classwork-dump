k = 10^-6:10^-6:10^-3;
L = 10^7;
l = 10^6;
tao = 0.02:0.01:0.2;
K = 10^3;
rho = 1025;
f = -10^(-4);
s = tao./(rho*f*K);
z = 0:100:4000;
z = z*-1;
psi = zeros(length(k),length(tao));

for i = 1:length(k)
   for j = 1:length(tao)
       for n = 1:length(z)
           %n = 1;
           depth(n) = (k(i)*L*(4*exp(z(n)/(s(j)*l))+1/2)/(4*s(j)*l*exp(z(n)/(s(j)*l))+z(n)));
       end
       psi(i,j) = max(depth);
   end
end

contour(tao,k,psi)
xlabel('tao')
ylabel('k')
%set(gca,'yscale','log')