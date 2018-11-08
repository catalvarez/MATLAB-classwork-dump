function f = StatsHW5(a)
% input a = rho0,K_prime,K0

mat = dlmread('ps5_data2.txt');
f = 0.5.*sum((mat(:,2)-a(1).*(1+a(2).*mat(:,1)./a(3)).^(1./a(2))).^2./mat(:,3).^2);

end