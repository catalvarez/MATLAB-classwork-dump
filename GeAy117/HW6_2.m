function f = HW6_2(a)
% input a = rho0,K_prime,K0

mat = dlmread('ps6_rv_data.txt');
f = sum((mat(:,2)-a(1).*cos(2*pi/6.5.*mat(:,1)+a(2))+a(3)).^2./mat(:,3).^2);

end